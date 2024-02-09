# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

from shxarray.core.logging import logger
from shxarray.core.admin import defaultcache
import xarray as xr
import os
import requests
import json

def download_snrei(dbfilename):
    if os.path.exists(dbfilename):
        logger.info(f"{dbfilename} already exists, no need to download)")
        return
    else:
        logger.info(f" Downloading {dbfilename}")

    url="https://github.com/strawpants/snrei/raw/master/Love/geoslurp_dump_llove.sql"
    req=requests.get(url)
    with open(dbfilename,'wb') as fid:
        fid.write(req.content)


def change_frame(dsLove,frame):
    ds1=dsLove.sel(n=1)
    if frame == "CF":
        alpha=((ds1.hn+2*ds1.ln)/3).item()
    elif frame == "CM":
        alpha=(1+ds1.kn).item()
    elif frame == "CE":
        alpha=ds1.kn.item()
    elif frame == "CH":
        alpha=ds1.hn.item()
    elif frame == "CL":
        alpha=ds1.ln.item()
    else:
        raise RuntimeError(f"Unknown frame {frame} selected")

    return dsLove.where(dsLove.n != 1,dsLove-alpha)



def snrei_load(model,dbfile_or_con,frame="CF",nmax=None):
    """Loads (load) Love numbers"""
    if dbfile_or_con is None:
        #create a default filename
        dbfile_or_con=os.path.join(defaultcache('Love'),"geoslurp_dump_llove.sql")
    
    if str(dbfile_or_con).endswith(".sql"):
        download_snrei(dbfile_or_con)
        import sqlite3
        qry=f"SELECT data from llove WHERE name = '{model}'"
        dbcon = sqlite3.connect(dbfile_or_con)
        res=json.loads(dbcon.execute(qry).fetchone()[0])
    else:
        #assume the input is already a open (sqlalchemy-like ) database connection
        qry=f"SELECT data from earthmodels.llove WHERE name = '{model}'"
        res=dbfile_or_con.execute(qry).first()[0]

    dsout=xr.Dataset.from_dict(res).rename({"degree":"n"})
    if nmax is not None:
        dsout=dsout.sel(n=slice(0,nmax))
    #possibly convert to a different isomorphic reference frame
    dsout=change_frame(dsout,frame)
    return dsout
    



class SnreiFactory:
    @staticmethod
    def load(model="PREM",dbfile_or_con=None,frame="CF",nmax=None):
        return snrei_load(model=model,dbfile_or_con=dbfile_or_con,frame=frame,nmax=nmax)



