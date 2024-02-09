# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2024
#

from shxarray.core.admin import defaultcache
import requests
import pandas as pd
from shxarray.io.binv_legacy import readBINV
from shxarray.kernels.anisokernel import AnisoKernel
import os



def load_catalogue():
    catalogfile=os.path.join(defaultcache('grace-filter'),'inventory.xlsx')
    if not os.path.exists(catalogfile):
        #dowload from github
        url="https://github.com/strawpants/GRACE-filter/raw/master/inventory.xlsx"
        r=requests.get(url)
        with open(catalogfile,'wb') as fid:
            fid.write(r.content)

    dfcat=pd.read_excel(catalogfile,index_col='name')
    return dfcat


def load_ddk(ddkversion,trans=False,nmax=-1,truncate=True):

    dfcat=load_catalogue()
    url=dfcat.loc[ddkversion].uri
    ddkfile=os.path.join(defaultcache('grace-filter'),os.path.basename(url))
    
    if not os.path.exists(ddkfile):
        #dowload from github
        r=requests.get(url)
        with open(ddkfile,'wb') as fid:
            fid.write(r.content)
    df=readBINV(ddkfile,trans,nmax)
    
    return AnisoKernel(df,name=ddkversion,truncate=truncate)

        



