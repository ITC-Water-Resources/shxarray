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

# Author Roelof Rietbroek (roelof@geod.uni-bonn.de), 2025
from geoslurp.dataset import DataSet
from geoslurp.datapull.http import Uri as http
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column,Integer,String,MetaData
from sqlalchemy.dialects.postgresql import TIMESTAMP,JSONB
from geoslurp.datapull.github import Crawler as ghCrawler
from geoslurp.datapull.github import GithubFilter as ghfilter
from geoslurp.config.slurplogger import slurplogger
from geoslurp.types.json import DataArrayJSONType
from geoslurp.datapull import UriFile
from geoslurp.datapull import findFiles
import numpy as np
import xarray as xr
import gzip
import re
import os

def lloveMetaExtractor(uri):
    """extract some metainfo from the load Lovenumber file"""
    #extract maximum degree from file and heuristically derive loadtype from the filename)
    if re.search("body",uri.url):
        ltype="body"
    else:
        ltype="surface"
    
    nmax=0
    reentry=re.compile('^ *[0-9]')
    hn=[]
    ln=[]
    kn=[]
    deg=[]
    slurplogger().info(f"Processing {uri.url}")
    descr=""
    ref=None
    with gzip.open(uri.url,'rt') as fid:
        for line in fid:
            if reentry.search(line):
                linespl=line.split()
                n=int(linespl[0])
                if n == 1:
                    #look for CF degree 1 coefficients only
                    ref="CF"
                    if linespl[4] != ref:
                        #only use the degree 1 numbers of the chosen reference system
                        continue
                deg.append(n)
                hln=[float(el.replace('D','E')) for el in linespl[1:4]]
                #possibly replace infinity values with NaN
                hln=[None if np.isinf(el) else el for el in hln]

                hn.append(hln[0])
                ln.append(hln[1])
                kn.append(hln[2])
            else:
                #append comment to description
                descr+=line

    #create an xarray dataset
    dslove=xr.Dataset(data_vars=dict(kn=(["degree"],kn),hn=(["degree"],hn),ln=(["degree"],ln)),coords=dict(degree=(["degree"],deg)))
    
    #extract the maximum degree
    nmax=dslove.degree.max().data.item()
    meta={"name":os.path.basename(uri.url).replace(".love.gz",""),"lastupdate":uri.lastmod,
            "descr":descr,"loadtype":ltype,"nmax":nmax,"ref":ref,"data":dslove}


    return meta

schema='shxarray'
LLoveTBase=declarative_base(metadata=MetaData(schema=schema))

class LLoveTable(LLoveTBase):
    """Defines the Load Love number PostgreSQL table"""
    __tablename__='llove'
    id=Column(Integer,primary_key=True)
    name=Column(String,unique=True)
    loadtype=Column(String)
    ref=Column(String)
    lastupdate=Column(TIMESTAMP)
    nmax=Column(Integer)
    descr=Column(String)
    data=Column(DataArrayJSONType)

class LLove(DataSet):
    """Class for registering load love numbers (downloads from github) """
    schema=schema
    version=(1,0,0)
    table=LLoveTable
    def __init__(self,dbconn):
        super().__init__(dbconn)
    
    def pull(self):
        """Pulls the dataset from github and unpacks it in the cache directory"""
        #crawls the github repository for Load Love numbers

        reponame="strawpants/snrei"
        commitsha="e61d3e2a9eb328d48f56f5aa73fa2aaba60f1d5c"
    
        try:
            cred=self.conf.authCred("github",['oauthtoken'])
            token=cred.oauthtoken
        except:
            token=None
        ghcrawler=ghCrawler(reponame,commitsha=commitsha,
                           filter=ghfilter({"type":"blob","path":"\.love"}),
                           followfilt=ghfilter({"type":"tree","path":"Love"}),
                           oauthtoken=token)
        
        #download all datasets
        ghcrawler.parallelDownload(self.cacheDir(),check=True,maxconn=3,gzip=True)

    def register(self):
        slurplogger().info("Building file list..")
        files=[UriFile(file) for file in findFiles(self.cacheDir(),'.*love',self._dbinvent.lastupdate)]

        if len(files) == 0:
            slurplogger().info("LLove: No new files found since last update")
            return
        
        self.truncateTable()
        #loop over files
        for uri in files:
            self.addEntry(lloveMetaExtractor(uri))
        self.updateInvent()









