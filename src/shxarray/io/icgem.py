# This file is part of frommle2.
# frommle2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.

# frommle2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with Frommle; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

# Author Roelof Rietbroek (r.rietbroek@utwente.nl), 2021

import gzip
import xarray as xr
import re
import sys
import numpy as np
from shxarray.core.sh_indexing import SHindexBase
from shxarray.core.logging import shxlogger 
from datetime import datetime,timedelta
from shxarray.core.cf import get_cfatts
import pandas as pd
from functools import partial

def get_gfc(lnspl,errors,nmaxstop=sys.maxsize):
    #parse a gfc line (ascii)
    n=int(lnspl[1])
    m=int(lnspl[2])
    
    if n > nmaxstop:
        return None,None
    
    c=dict(ityp=lnspl[0],n=n,m=m,cnm=float(lnspl[3]))

    if m != 0:
        s=dict(ityp=lnspl[0],n=n,m=-m,cnm=float(lnspl[4]))
    else:
        s=None
    
    if errors == 1:
        c['sigcnm']=float(lnspl[5])
        if m != 0:
            s['sigcnm']=float(lnspl[6])
    
    return c,s



def get_gfct(lnspl,errors,nmaxstop=sys.maxsize):
    #parse a gfc line (ascii)
    n=int(lnspl[1])
    m=int(lnspl[2])
    if n > nmaxstop:
        return None,None

    time=datetime.strptime(lnspl[5+errors*2],'%Y%m%d')
    c=dict(ityp=lnspl[0],n=n,m=m,cnm=float(lnspl[3]),t0=time)

    if m != 0:
        s=dict(ityp=lnspl[0],n=n,m=-m,cnm=float(lnspl[4]),t0=time)
    else:
        s=None
    
    if errors == 1:
        c['sigcnm']=float(lnspl[5])
        if m != 0:
            s['sigcnm']=float(lnspl[6])
    
    return c,s


def get_trig(lnspl,errors):
    #parse a acos/asin line (ascii)
    
    n=int(lnspl[1])
    m=int(lnspl[2])
    period=float(lnspl[5+errors*2])
    c=dict(ityp=lnspl[0],n=n,m=m,cnm=float(lnspl[3]),period_yr=period)

    if m != 0:
        s=dict(ityp=lnspl[0],n=n,m=-m,cnm=float(lnspl[4]),period_yr=period)
    else:
        s=None
    
    if errors == 1:
        c['sigcnm']=float(lnspl[5])
        if m != 0:
            s['sigcnm']=float(lnspl[6])
    
    return c,s

def readIcgem(fileobj,nmaxstop=sys.maxsize):
    
    if nmaxstop is None:
        nmaxstop=sys.maxsize
    
    needsClosing=False
    if type(fileobj) == str:
        needsClosing=True
        if fileobj.endswith('.gz'):
            fileobj=gzip.open(fileobj,'rt')
        else:
            fileobj=open(fileobj,'rt')

    hassigma=False
    #first read the icgem header 
    hdr={}
    for ln in fileobj:
        if 'begin_of_head' in ln:
            continue
        if 'end_of_head' in ln:
            break
        
        spl=ln.split()
        if len(spl) == 2:
            #insert name value pairs in the hdr dict
            hdr[spl[0]]=spl[1]
        elif len(spl) > 4 and spl[0] == 'key':
            if "sigma" in spl:
                hassigma=True
        shxlogger.warning("product_type not specified in icgem header, assuming gravity_field")
    #extract relevant parameters from the header
    attr={}
    try:
        nmaxsupp=int(hdr["max_degree"])
        attr["nmaxfile"]=nmaxsupp
        if nmaxsupp < nmaxstop:
            attr["nmax"]=nmaxsupp
        else:
            attr["nmax"]=nmaxstop
        nmax=attr["nmax"]

        if nmax > nmaxsupp:
            shxlogger.warning("Nmax ({nmax}) requested larger than supported, higher degree coefficients will be set to zero")


        if 'format' in hdr:
            attr["format"]=hdr['format']
        else:
            attr["format"]="icgem1.0"
        
        if "norm" in hdr:
            attr["norm"]=hdr["norm"]
        
        attr["gm"]=float(hdr["earth_gravity_constant"].replace("D","E"))
        attr["re"]=float(hdr["radius"])
        attr["modelname"]=hdr["modelname"]
    except KeyError:
    #some values may not be present but that is ok
        pass
    
    try:
        if hdr['product_type'] != 'gravity_field':
            raise ValueError(f"Only gravity_field product_type is supported, not {hdr['product_type']}")
    except KeyError:
        shxlogger.warning("product_type not specified in icgem header, assuming gravity_field")

    #Non standard HACK to try to retrieve the epoch from the modelname (GRAZ monthly solutions only)
    if "modelname" in hdr:
        try:
            time=[datetime.strptime(hdr['modelname'][-7:],"%Y-%m")+timedelta(days=14.5)]
        except ValueError:
            time=None
    else:
        time=None
   


    if attr["format"] != "icgem1.0":
        raise ValueError(f"Only icgem1.0 format is supported, not {attr['format']}")

    parser={}
    try:
        if hdr['errors'] == 'no':
            errors=0
        elif hdr['errors'] == 'formal' or hdr['errors'] == 'calibrated':
            errors=1
        else:
            raise ValueError(f"Cannot handle error specification {hdr['error']} in icgem header")
    except KeyError: 
        shxlogger.warning("errors specification not found in icgem header, assuming no errors")
        errors=0 


    parser['gfc']=partial(get_gfc,errors=errors,nmaxstop=nmaxstop)
    parser['gfct']=partial(get_gfct,errors=errors,nmaxstop=nmaxstop)
    parser['trnd']=partial(get_gfc, errors=errors)
    parser['asin']=partial(get_trig, errors=errors)
    parser['acos']=partial(get_trig, errors=errors)

    rowdicts=[]
    for ln in fileobj:
        lnspl=ln.replace("D","E").split()
        ky=lnspl[0]
        try:
            c,s=parser[ky](lnspl)
            if c is None and s is None:
                #stop parsing
                break
            rowdicts.append(c)
            if s is not None:
                rowdicts.append(s)

        except KeyError:
            #ok (not supported key)
            shxlogger.warning(f"key not supported {ky}, ignoring")


    if needsClosing:
        fileobj.close()
    
    #stuff everything in a dataframe
    df=pd.DataFrame.from_dict(rowdicts)
    
    #convert ot xarray
    ds=df.to_xarray().set_coords(['n','m','ityp']).set_xindex('ityp')
    ds=ds.drop_vars('index').rename(dict(index='nm'))
    
    itypes=df.ityp.unique()
    if len(itypes) == 1:
        if itypes[0] !=  'gfc':
            raise ValueError(f"Only {itypes[0]} coefficients found in the file, don't know how to handle")

        dsout=ds.drop_vars(['ityp'])#.sh.build_nmindex()
        
        if time is not None:
            dsout['cnm']=dsout.cnm.expand_dims('time')
            if errors == 1:
                dsout['sigcnm']=dsout.sigcnm.expand_dims('time')
            dsout=dsout.assign_coords(time=time)
        dsout=dsout.sh.build_nmindex()
    else:
        dsout=None 
        for ityp,dsgrp in ds.groupby('ityp'):
            
            if ityp in ['acos','asin']:
                dsgrp=dsgrp.drop_vars(['ityp','t0'])
            elif ityp == 'trnd':
                dsgrp=dsgrp.drop_vars(['ityp','t0','period_yr'],errors='ignore')
            elif ityp == 'gfct':
                dsgrp=dsgrp.drop_vars(['ityp','period_yr'],errors='ignore')
            else:
                dsgrp=dsgrp.drop_vars(['ityp','t0','period_yr'],errors='ignore')
            
            #add attributes
            dsgrp['cnm'].attrs.update(get_cfatts("stokes"))

            if errors == 1:
                dsgrp['sigcnm'].attrs.update(get_cfatts("stokes stdv"))
                renamedict=(dict(cnm=f'cnm_{ityp}',sigcnm=f'sigcnm_{ityp}'))
            else:
                renamedict=(dict(cnm=f'cnm_{ityp}'))
            dsgrp=dsgrp.rename(renamedict).sh.build_nmindex()
            if dsout is None:
                dsout=dsgrp
            else:
                #add the data to the output dataset
                if dsout.sizes['nm'] != dsgrp.sizes['nm']:
                    #add the data but it has it's own nm index
                    dsgrp=dsgrp.rename(dict(nm=f'nm_{ityp}',n=f'n_{ityp}',m=f'm_{ityp}'))

            
            dsout=dsout.merge(dsgrp)
    
    dsout.attrs.update(attr)

    return dsout




# def get_icgem_at_time(dsicgem,time,nmax=None):
    # """
    # #extract a set of Stokes coefficients at a given time from an ICGEM dataset with time variable components
    # """

    # import pdb;pdb.set_trace()
    
    # if nmax is None:
        # nmax=dsicgem.sh.nmax
    # if "
    # dsout=dsicgem.sh.slice_nm(nmax=nmax).sh.build_nmindex()
