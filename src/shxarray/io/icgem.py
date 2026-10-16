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

def parse_dt(data):
        """
            parse time epoch strings into datetime objects
            Note: strptime with %Y%m%d.%H%M doesn't work as there are files out there which use 60 for the minutes
        """
        dtag=data.split(".")
        if len(dtag) == 1 :
            dt=datetime.strptime(dtag[0],'%Y%m%d')
        else: 
            dtag,hmtag=dtag
            dt=datetime.strptime(dtag,'%Y%m%d')
            dt+=timedelta(hours=int(hmtag[0:2]),minutes=int(hmtag[2:4]))


        return dt

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

    time=parse_dt(lnspl[5+errors*2])
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

def get_gfct2(lnspl,errors,nmaxstop=sys.maxsize):
    #parse a gfc line (ascii)
    n=int(lnspl[1])
    m=int(lnspl[2])
    if n > nmaxstop:
        return None,None
    tstart=parse_dt(lnspl[5+errors*2])
    tend=parse_dt(lnspl[6+errors*2])
    c=dict(ityp=lnspl[0],n=n,m=m,cnm=float(lnspl[3]),t0=tstart,t1=tend)

    if m != 0:
        s=dict(ityp=lnspl[0],n=n,m=-m,cnm=float(lnspl[4]),t0=tstart,t1=tend)
    else:
        s=None
    
    if errors == 1:
        c['sigcnm']=float(lnspl[5])
        if m != 0:
            s['sigcnm']=float(lnspl[6])
    return c,s



def get_trig(lnspl,errors,nmaxstop=sys.maxsize):
    #parse a acos/asin line (ascii)
    
    n=int(lnspl[1])
    m=int(lnspl[2])
    if n > nmaxstop:
        return None,None
    
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

def get_trig2(lnspl,errors,nmaxstop=sys.maxsize):
    #parse a acos/asin line (ascii)
    n=int(lnspl[1])
    m=int(lnspl[2])
    if n > nmaxstop:
        return None,None
    
    tstart=parse_dt(lnspl[5+errors*2])
    tend=parse_dt(lnspl[6+errors*2])

    period=float(lnspl[7+errors*2])
    c=dict(ityp=lnspl[0],n=n,m=m,cnm=float(lnspl[3]),period_yr=period,t0=tstart,t1=tend)

    if m != 0:
        s=dict(ityp=lnspl[0],n=n,m=-m,cnm=float(lnspl[4]),period_yr=period,t0=tstart,t1=tend)
    else:
        s=None
    
    if errors == 1:
        c['sigcnm']=float(lnspl[5])
        if m != 0:
            s['sigcnm']=float(lnspl[6])
    
    return c,s

def reorganize_ds_simple(ds,time,errors):
    if ds.ityp[0] !=  'gfc':
        raise ValueError(f"Only {ds.ityp[0]} coefficients found in the file, don't know how to handle")

    dsout=ds.drop_vars(['ityp'])#.sh.build_nmindex()
    
    if time is not None:
        dsout['cnm']=dsout.cnm.expand_dims('time')
        if errors == 1:
            dsout['sigcnm']=dsout.sigcnm.expand_dims('time')
        dsout=dsout.assign_coords(time=time)
    dsout=dsout.sh.build_nmindex()
    return dsout

def reorganize_ds1(ds,errors):
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
    return dsout


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
    inheader=False
    for ln in fileobj:
        if 'begin_of_head' in ln:
            inheader=True
            continue
        if 'end_of_head' in ln:
            break
        
        if not inheader:
            continue

        spl=ln.split()
        if len(spl) == 2:
            #insert name value pairs in the hdr dict
            hdr[spl[0]]=spl[1]
        elif len(spl) > 4 and spl[0] == 'key':
            if "sigma" in spl:
                hassigma=True
        elif len(spl) > 2 and spl[0] == "errors":
            #Exception Clause as some data providers do not comply with adding a single parameter in the header
            hdr[spl[0]]=spl[1]
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
        #assume gravity field
        hdr['product_type'] = 'gravity_field'
        # shxlogger.warning("product_type not specified in icgem header, assuming gravity_field")

    #Non standard HACK to try to retrieve the epoch from the modelname (GRAZ monthly solutions only)
    if "modelname" in hdr:
        try:
            time=[datetime.strptime(hdr['modelname'][-7:],"%Y-%m")+timedelta(days=14.5)]
        except ValueError:
            time=None
    else:
        time=None
   


    if attr["format"] not in ["icgem1.0","icgem2.0"]:
        raise ValueError(f"Only icgem1.0/2.0 format is supported, not {attr['format']}")

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
    if attr["format"] == "icgem1.0":
        parser['gfct']=partial(get_gfct,errors=errors,nmaxstop=nmaxstop)
        parser['asin']=partial(get_trig, errors=errors,nmaxstop=nmaxstop)
        parser['acos']=partial(get_trig, errors=errors,nmaxstop=nmaxstop)
        parser['trnd']=partial(get_gfc, errors=errors,nmaxstop=nmaxstop)
    elif attr["format"] == "icgem2.0":
        parser['gfct']=partial(get_gfct2,errors=errors,nmaxstop=nmaxstop)
        parser['asin']=partial(get_trig2, errors=errors,nmaxstop=nmaxstop)
        parser['acos']=partial(get_trig2, errors=errors,nmaxstop=nmaxstop)
        parser['trnd']=partial(get_gfct2, errors=errors,nmaxstop=nmaxstop)


    rowdicts=[]
    for ln in fileobj:
        lnspl=ln.replace("D","E").split()
        ky=lnspl[0]
        try:
            c,s=parser[ky](lnspl)
            if c is None and s is None:
                continue

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
    
    #convert to xarray
    ds=df.to_xarray().set_coords(['n','m','ityp']).set_xindex('ityp')
    ds=ds.drop_vars('index').rename(dict(index='nm'))
    
    itypes=df.ityp.unique()
    if len(itypes) == 1:
        dsout=reorganize_ds_simple(ds,time,errors)
    elif attr["format"] in ["icgem1.0","icgem2.0"]:
        # do not attempt to reorganize icgem2.0 data in blocks, as it is too fragmented
        dsout=ds
    else:
        raise RuntimeError("cannot reorganize icgem data")
    
    dsout.attrs.update(attr)

    return dsout

def sample_tv_grav(dsicgem,time,nmax=None, excludetypes=None):
    """
    Extract a set of Stokes coefficients at a given time from an ICGEM dataset with time variable components and intervals
    """
    if 'ityp' not in dsicgem:
        raise RuntimeError("The input dataset does not seem to originate from an icgem2.0 file")
    nmaxsupp=dsicgem.n.max().item()

    if nmax is None:
        nmax=nmaxsupp
    elif nmax > nmaxsupp:

        shxlogger.warning(f"Nmax ({nmax}) requested larger than supported, truncating to {nmaxsup}")
        nmax=nmaxsupp
    if not np.isscalar(time):
        #assume time is array like and call this function iteratively
        dsout=xr.DataArray.sh.zeros(nmax=nmax,auxcoords={"time":time}).to_dataset()
        dsout.attrs=dsicgem.attrs
        for t in time:
            # print(f"working on {t}")
            dssub=sample_tv_grav(dsicgem,time=t,nmax=nmax,excludetypes=excludetypes)
            dsout.cnm.loc[t,dssub.nm]=dssub.cnm[0,:]

        return dsout

    #dataarray to collect all the cnm coefficients
    #initiate with timeconstant gfc coefficients
    dacollect=xr.DataArray.sh.zeros(nmax)
    dacollect.name='cnm'
    dagfc=dsicgem.cnm.sel(ityp='gfc').drop_vars('ityp')
    if dagfc.size == 1:
        #single entry only
        dacollect.loc[dict(nm=(dagfc.n.item(),dagfc.m.item()))]=dagfc.item()
    elif dagfc.size > 1:
        dagfc=dagfc.sh.build_nmindex().sh.truncate(nmax) 
        dacollect.loc[dagfc.nm]=dagfc
    
    # get an index vector for the time variable components
    if dsicgem.attrs['format'] == 'icgem2.0':
        idx=(dsicgem.t0 <= time) & (dsicgem.t1 > time)
        dssub=dsicgem.sel(nm=idx)

    else:
        #select all coefficients as they are valid for all intervals
        dssub=dsicgem
        #add t0 entries
        t0=[t0 for t0 in np.unique(dsicgem.t0) if ~np.isnat(t0)]
        if len(t0) != 1:
            raise ValueError("Cannot figure out unique t0")

        #overwrite values
        dssub.t0.loc[:]=t0[0]


    days_per_year=365.25 
    ohm_an=2*np.pi

    for ityp, dsgrp in dssub.groupby('ityp'):
        
        if excludetypes is not None and ityp in excludetypes:
            continue
        
        if ityp == 'acos':
            t0=xr.DataArray([np.datetime64(datetime(yr,1,1)) for yr in dsgrp.t0.dt.year.data],dims=['nm'])
            dt_yr= (time-t0)/np.timedelta64(int(days_per_year*24),'h')
            fac=np.cos(ohm_an*dt_yr/dsgrp.period_yr)
            # breakpoint()
            dsub=(dsgrp.cnm*fac).drop_vars('ityp').sh.build_nmindex().groupby('nm').sum()
        elif ityp == 'asin':
            t0=xr.DataArray([np.datetime64(datetime(yr,1,1)) for yr in dsgrp.t0.dt.year.data],dims=['nm'])
            dt_yr= (time-t0)/np.timedelta64(int(days_per_year*24),'h')
            # dt_yr= (time-t0).dt.total_seconds()/86400/days_per_year
            fac=np.sin(ohm_an*dt_yr/dsgrp.period_yr)
            dsub=(dsgrp.cnm*fac).drop_vars('ityp').sh.build_nmindex().groupby('nm').sum()
        elif ityp == 'trnd':
            t0=xr.DataArray([np.datetime64(datetime(yr,1,1)) for yr in dsgrp.t0.dt.year.data],dims=['nm'])
            # fac=(time-t0).dt.total_seconds()/86400/days_per_year
            dt_yr= (time-t0)/np.timedelta64(int(days_per_year*24),'h')
            dsub=(dsgrp.cnm*dt_yr).drop_vars('ityp').sh.build_nmindex()
        elif ityp == 'gfct':
            dsub=(dsgrp.cnm).drop_vars('ityp').sh.build_nmindex()
        dacollect.loc[dsub.nm]+=dsub

    dacollect.attrs.update(get_cfatts("stokes stdv"))
    
    dsout=dacollect.to_dataset()
    dsout.attrs=dsicgem.attrs
    dsout=dsout.expand_dims(time=[time])
    return dsout
