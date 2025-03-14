# This file is part of geoslurp.
# geoslurp is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.

# geoslurp is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with geoslurp; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

# Author Roelof Rietbroek (r.rietbroek@utwente.nl), 2021
from geoslurp.dataset.pandasbase import PandasBase
from geoslurp.config.slurplogger import slurplogger
from geoslurp.datapull.http import Uri as http
from geoslurp.datapull import UriFile
import os
from datetime import datetime
import pandas as pd

schema='shxarray'

class GRACEfilter(PandasBase):
    """Class for registering SH filters (downloads from github) """
    schema=schema
    version=(0,0,0)
    def __init__(self,dbconn):
        super().__init__(dbconn)
        self.pdfile=os.path.join(self.cacheDir(),'inventory_upd.csv') 
    def pull(self):
        """Pulls the dataset from github and unpacks it in the cache directory"""
        #download the inventory file 
        lastchanged=datetime(2021,11,5)
        inventory="https://github.com/strawpants/GRACE-filter/raw/master/inventory.xlsx"
        uri,upd=http(inventory,lastmod=lastchanged).download(self.cacheDir(),check=True)
        pdinvent=pd.read_excel(uri.url,engine="openpyxl")
        #download all the files
        ddir=self.dataDir()
        for idx,row in pdinvent.iterrows():
            ffile,upd=http(row.uri,lastmod=lastchanged).download(ddir,check=True)
            #update file with newly downloaded file
            pdinvent.at[idx,'uri']=self.conf.generalize_path(ffile.url)

        #write updated excel file
        pdinvent.to_csv(os.path.join(self.pdfile))

