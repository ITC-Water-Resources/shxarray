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
# License along with Frommle; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

# Author Roelof Rietbroek (roelof@geod.uni-bonn.de), 2025

from geoslurp.dataset.dataSetBase import DataSet
from geoslurp.datapull.ftp import Crawler as ftpCrawler
from geoslurp.config.slurplogger import slurplog
from geoslurp.datapull import findFiles
from glob import glob
import gzip
import yaml
from geoslurp.datapull import UriFile
from io  import StringIO
import os
from datetime import datetime
from shxarray.geoslurp.gravity import GravitySHTBase
import re

schema="shxarray"

def graceMetaExtractor(uri):
    """Extract meta information from a GRACE file"""

    #some dirty search rand replace hacks to fix faulty yaml header in the grace/fo data
    hdrpatches=[(re.compile("0000-00-00T00:00:00"),"1970-01-01T00:00:00"),
            (re.compile(r"Dahle et al\. \(2019\)\:"),"Dahle et al. (2019),"),
            (re.compile(r"Dobslaw et al\. \(2019\)\:"),"Dobslaw et al. (2019),")]
    patchedLines=0

    buf=StringIO()
    with gzip.open(uri.url,'rt') as fid:
        slurplog.info("Extracting info from %s"%(uri.url))
        for ln in fid:
            if '# End of YAML header' in ln:
                #parse the yaml header
                hdr=yaml.safe_load(buf.getvalue())["header"]
                break
            else:
                # if re.search("Dahle",ln):
                    # import pdb;pdb.set_trace()
                #see if the line needs patching
                for reg,repl in hdrpatches:
                    ln,nr=re.subn(reg,repl,ln,count=1)
                    patchedLines+=nr
                #hack replace 0000-00-00 dates because yaml can't parse them
                buf.write(ln)
        if patchedLines > 0:
            #we want to fix the header and patch the input file
            buf.write(ln) #write end of YAML file
            #dunp the remainder of the file in the stringio buffer
            buf.write(fid.read())

    
    if patchedLines > 0:
        slurplog.info("Patching faulty yaml header in file %s"%uri.url)
        with gzip.open(uri.url,'wt') as fidout:
            fidout.write(buf.getvalue())

    nonstand=hdr["non-standard_attributes"]


    meta={"nmax":hdr["dimensions"]["degree"],
          "omax":hdr["dimensions"]["order"],
          "tstart":hdr["global_attributes"]["time_coverage_start"],
          "tend":hdr["global_attributes"]["time_coverage_end"],
          "lastupdate":uri.lastmod,
          "format":nonstand["format_id"]["short_name"],
          "gm":nonstand["earth_gravity_param"]["value"],
          "re":nonstand["mean_equator_radius"]["value"],
          "uri":uri.url,
          "type":nonstand["product_id"][0:3],
          "data":{"description":hdr["global_attributes"]["title"]}
          }

    #add tide system
    try:
        tmp=nonstand["permanent_tide_flag"]
        if re.search('inclusive',tmp):
            meta["tidesystem"]="zero-tide"
        elif re.search('exclusive'):
            meta["tidesystem"]="tide-free"
    except:
        pass

    return meta

class GRACEL2Base(DataSet):
    """Derived type representing GRACE spherical harmonic coefficients on the podaac server"""
    release=None
    center=None
    updated=None
    schema=schema
    stripuri=True
    def __init__(self,dbconn):
        super().__init__(dbconn)
        #initialize postgreslq table
        GravitySHTBase.metadata.create_all(self.db.dbeng, checkfirst=True)

    def pull(self):
        url="ftp://isdcftp.gfz.de/"+self.mission+"/Level-2/"+self.center+"/"+self.release+"/"
        ftp=ftpCrawler(url,pattern='G.*gz$')

        self.updated=ftp.parallelDownload(self.dataDir(),maxconn=3,check=True)

    def register(self):

        #create a list of files which need to be (re)registered
        if self.updated:
            files=self.updated
        else:
            files=[UriFile(file) for file in findFiles(self.dataDir(),r'G.*\.gz$',self._dbinvent.lastupdate)]
        filesnew=self.retainnewUris(files)
        
        #loop over the newer files
        for uri in filesnew:
            meta=graceMetaExtractor(uri)
            self.addEntry(meta)
        
        self.updateInvent()


def GRACEL2ClassFactory(clsName):
    """Dynamically construct GRACE Level 2 dataset classes"""
    mission,center,release=clsName.split("_")
    clsname_ad=clsName.replace("-","_").replace(".","")
    table=type(clsname_ad+"Table", (GravitySHTBase,), {})
    return type(clsname_ad, (GRACEL2Base,), {"release": release, "center":center,"table":table,"mission":mission})

# setup GRACE datasets
def GRACEDsets(conf):
    out=[]
    

    combos=[("grace","CSR","RL06"),
            ("grace","GFZ","RL06"),
            ("grace","JPL","RL06"),
            ("grace-fo","CSR","RL06.3"),
            ("grace-fo","GFZ","RL06.3"),
            ("grace-fo","JPL","RL06.3")]

    for mission,center,release in combos:
        clsName=mission+"_"+center+"_"+release
        out.append(GRACEL2ClassFactory(clsName))
   
    return out





