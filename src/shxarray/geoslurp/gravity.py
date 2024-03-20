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

# Author Roelof Rietbroek (roelof@geod.uni-bonn.de), 2018

from sqlalchemy.ext.declarative import declared_attr, as_declarative
from sqlalchemy import MetaData
from sqlalchemy import Column,Integer,String, Boolean,Float
from sqlalchemy.dialects.postgresql import TIMESTAMP, JSONB
from geoslurp.types.json import DataArrayJSONType


schema="shxarray"

#define a declarative baseclass for spherical harmonics gravity  data
@as_declarative(metadata=MetaData(schema=schema))
class GravitySHTBase(object):
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
    omax=Column(Integer)
    gm=Column(Float)
    re=Column(Float)
    tidesystem=Column(String)
    origin=Column(String)
    format=Column(String)
    type=Column(String)
    uri=Column(String)
    data=Column(JSONB)

@as_declarative(metadata=MetaData(schema=schema))
class GravitySHinDBTBase(object):
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
    omax=Column(Integer)
    gm=Column(Float)
    re=Column(Float)
    tidesystem=Column(String)
    origin=Column(String)
    format=Column(String)
    type=Column(String)
    data=Column(DataArrayJSONType)
    __table_args__ = {'extend_existing': True}


class Trig():
    """Enum to distinguish between a trigonometric cosine and sine coefficient"""
    c = 0
    s = 1

class JSONSHArchive():
    """JSON Archive which stores SH data, with sigmas and possibly a covariance
    Note this mimics the Archive interface of frommle without actually requiring its import"""
    def __init__(self,nmax=None,datadict=None):
        
        if nmax:
            #create from maximum degree
            self.data_={"attr":{},"vars":{"shg":[]}}
            shg=[]
            for n in range(nmax+1):
                for m in range(n+1):
                    shg.append((n,m,Trig.c))
                    if m> 0:
                        shg.append((n,m,Trig.s))
            self.data_["vars"]["shg"]=shg
            self.data_["attr"]["nmax"]=nmax
        elif not nmax and datadict:
            self.data_=datadic
        else:
            raise RunTimeError("Can only construct a JSONSHArchive from either nmax or a datadict")

    def __getitem__(self,key):
        """retrieves a named variable, and lazily creates allowed variables when requested"""
        if not key in self.data_["vars"]:
            if key in ["cnm","sigcnm"]:
                self.data_["vars"][key]=[0]*len(self.data_["vars"]["shg"])
            elif key == "covcnm":
                nl=len(self.data_["vars"]["shg"])
                self.data_["vars"][key]=[[0]*nl]*nl
        return self.data_["vars"][key]
    
    def __setitem__(self,key,item):
        self.data_["vars"][key]=item


    def idx(self,nmt):
        """returns the index of the n,m,t tuple"""
        return self.data_["vars"]["shg"].index(nmt)

    @property
    def attr(self):
        """get the stored global attributes of the file"""
        return self.data_["attr"]

    @attr.setter
    def attr(self,attrdict):
        """sets the stored global attributes of the file"""
        self.data_["attr"]=attrdict
    
    @property    
    def dict(self):
        return self.data_
    
