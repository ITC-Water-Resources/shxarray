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

from sqlalchemy.ext.declarative import declared_attr, as_declarative
from sqlalchemy import MetaData
from sqlalchemy import Column,Integer,String, Boolean,Float
from sqlalchemy.dialects.postgresql import TIMESTAMP
from geoslurp.dataset import DataSet
from geoslurp.datapull.ftp import Crawler as ftpCrawler
from geoslurp.datapull.uri import findFiles
from geoslurp.config.slurplogger import slurplogger
from geoslurp.datapull import UriFile
from datetime import datetime,timedelta
from shxarray.geoslurp.gravity import GravitySHTBase,GravitySHinDBTBase
from shxarray.geoslurp.icgem import icgemMetaExtractor
import re
import os
from geoslurp.db.settings import getCreateDir
import gzip as gz

schema='shxarray'

@as_declarative(metadata=MetaData(schema=schema))
class GravityNormals(object):
    @declared_attr
    def __tablename__(cls):
        #strip of the 'Table' from the class name
        return cls.__name__[:-5].lower()
    id = Column(Integer, primary_key=True)
    lastupdate=Column(TIMESTAMP)
    tstart=Column(TIMESTAMP,index=True)
    tend=Column(TIMESTAMP,index=True)
    time=Column(TIMESTAMP,index=True)
    nmax=Column(Integer)
    format=Column(String)
    uri=Column(String)

def snxdate2datetime(snxd):
    yr,doy,sec=[int(x) for x in snxd.split(":")]
    if yr > 50:
        yr+=1900
    else:
        yr+=2000

    return datetime(yr,1,1)+timedelta(days=doy-1,seconds=sec)

def NormalsMetaExtractor(uri):

    #read some stuff from the sinex files
    with gz.open(uri.url,'rt') as fid:
        slurplogger().info("Extracting info from %s"%(uri.url))
        ln=fid.readline()
        if ln[0:5] != '%=SNX':
            raise IOError("This doesn't seem like a valid SINEX file")
        lnspl=ln.split()
        tstart=snxdate2datetime(lnspl[5])
        tend=snxdate2datetime(lnspl[6])
        tcent=tstart+(tend-tstart)/2

    meta={"uri":uri.url,"lastupdate":uri.lastmod,"tend":tend,"tstart":tstart,"time":tcent}
    return meta


def enhanceMeta(meta):
    """Extract addtional timestamps from the TU graz filename data"""
    #from the TU GRAZ files one can extract the time period information
    yyyymm_match=re.match(".*([0-9]{4})-([0-9]{2})\.gfc\.gz$",meta["uri"])
    if yyyymm_match:
        yr=int(yyyymm_match.group(1))
        mn=int(yyyymm_match.group(2))
        meta["tstart"]=datetime(yr,mn,1)
        if mn ==12:
            meta['tend']=datetime(yr+1,1,1)
        else:
            meta["tend"]=datetime(yr,mn+1,1)
    else:
        yyyymmdd_match=re.match(".*([0-9]{4})-([0-9]{2})-([0-9]{2})\.gfc\.gz$",meta["uri"])
        if yyyymmdd_match:
            yr=int(yyyymmdd_match.group(1))
            mn=int(yyyymmdd_match.group(2))
            dd=int(yyyymmdd_match.group(3))
            meta['tstart']=datetime(yr,mn,dd)
            meta['tend']=datetime(yr,mn,dd,23,59)
    
    meta["time"]=meta["tstart"]+(meta["tend"]-meta["tstart"])/2
    if re.search("_background",meta["uri"]):
        meta["origin"]="CE"
    else:
        meta["origin"]="CM"
    return meta

class TUGRAZGRACEL2Base(DataSet):
    """Derived type representing GRACE spherical harmonic coefficients from the TU GRAZ"""
    table=None
    updated=None
    schema=schema
    release=''
    subdirs=''
    stripuri=True
    indb=False
    def __init__(self,dbconn):
        super().__init__(dbconn)
        if not self._dbinvent.datadir:
            self._dbinvent.datadir=getCreateDir(os.path.join(self.conf.getDataDir(self.schema),self.release,self.subdirs))

    def pull(self):
        url=os.path.join("ftp://ftp.tugraz.at/outgoing/ITSG/GRACE/",self.release,self.subdirs)
        ftp=ftpCrawler(url,pattern='.*.gfc',followpattern='([0-9]{4})')

        self.updated=ftp.parallelDownload(self._dbinvent.datadir,check=True,gzip=True, maxconn=3)

    def register(self):

        #create a list of files which need to be (re)registered
        if self.updated:
            files=self.updated
        else:
            files=[UriFile(file) for file in findFiles(self._dbinvent.datadir,'.*gfc.gz',since=self._dbinvent.lastupdate)]
        # breakpoint()
        newfiles=self.retainnewUris(files)
        #loop over files
        for uri in newfiles:
            slurplogger().info("extracting meta info from %s"%(uri.url))
            meta=icgemMetaExtractor(uri)
            meta=enhanceMeta(meta)
            self.addEntry(meta)

        self.updateInvent()


class TUGRAZGRACEL2NormalBase(DataSet):
    """Derived type representing GRACE Normal equation systems from the TU GRAZ"""
    table=None
    schema=schema
    release=''
    subdirs=''
    updated=False
    stripuri=True
    def __init__(self,dbconn):
        super().__init__(dbconn)

    def pull(self):
        url=os.path.join("ftp://ftp.tugraz.at/outgoing/ITSG/GRACE/",self.release,self.subdirs)
        ftp=ftpCrawler(url,pattern='.*.snx.gz',followpattern='([0-9]{4})')

        self.updated=ftp.parallelDownload(self.dataDir(),check=True, maxconn=5)

    def register(self):

        #create a list of files which need to be (re)registered
        if self.updated:
            files=self.updated
        else:
            files=[UriFile(file) for file in findFiles(self.dataDir(),'.*snx.gz',since=self._dbinvent.lastupdate)]

        newfiles=self.retainnewUris(files)
        #loop over files
        for uri in newfiles:
            slurplogger().info("extracting meta info from %s"%(uri.url))
            meta=NormalsMetaExtractor(uri)
            self.addEntry(meta)

        self.updateInvent()

def TUGRAZGRACEL2ClassFactory(release,subdirs,indb=False):
    """Dynamically construct GRACE Level 2 dataset classes for TU GRAZ"""
    base,gtype=subdirs.split('/')
    clsName="_".join([release,gtype])
    if indb:
        table=type(clsName.replace('-',"_") +"Table", (GravitySHinDBTBase,), {})
    else:
        table=type(clsName.replace('-',"_") +"Table", (GravitySHTBase,), {})
    return type(clsName, (TUGRAZGRACEL2Base,), {"release": release, "table":table,"subdirs":subdirs,"indb":indb})

def TUGRAZNormalsClassFactory(release,subdirs):
    """Dynamically construct GRACE Level 2 dataset classes for TU GRAZ"""
    base,_,gtype=subdirs.split('/')
    clsName="_".join([release,"normals",gtype])
    table=type(clsName.replace('-',"_") +"Table", (GravityNormals,), {})
    return type(clsName, (TUGRAZGRACEL2NormalBase,), {"release": release, "table":table,"subdirs":subdirs})

# setup GRACE datasets
def TUGRAZGRACEDsets(conf):
    out=[]
    release='ITSG-Grace2018'
    for subdirs in ["daily_kalman/daily_background","daily_kalman/daily_n40","monthly/monthly_background","monthly/monthly_n60","monthly/monthly_n96","monthly/monthly_n120"]:
        out.append(TUGRAZGRACEL2ClassFactory(release,subdirs))

    #and normal equation systems
    out.append(TUGRAZNormalsClassFactory(release,'monthly/normals_SINEX/monthly_n96')) 

    #also add GRACE-FO
    release="ITSG-Grace_operational"
    for subdirs in ["monthly/monthly_signalSeparation","monthly/monthly_background","monthly/monthly_n60","monthly/monthly_n96","monthly/monthly_n120"]:
        out.append(TUGRAZGRACEL2ClassFactory(release,subdirs))


    #and normal equation systems
    out.append(TUGRAZNormalsClassFactory(release,'monthly/normals_SINEX/monthly_n96')) 

    return out

