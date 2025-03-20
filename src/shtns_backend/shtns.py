# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2025

 
import shtns 
import numpy as np
import xarray as xr

from shxarray.core.cf import find_lon,find_lat,get_cfatts,get_cfglobal,change_central_longitude
from shxarray.core.shcomputebase import SHComputeBackendBase

def get_shtns_nm_order(nmax,nmin=0):
    #cosine coefficients (related to real part)
    nm_c=[(n,m) for m in range(nmax+1) for n in range(m,nmax+1)]
    #sine coefficients (related to imaginary part)
    nm_s=[(n,-m) for m in range(nmax+1) for n in range(m,nmax+1)]
   

    scalem0=np.sqrt(4*np.pi)
    scalemaux=scalem0/np.sqrt(2)
    conv=[scalem0 if m == 0 else scalemaux for _,m in nm_c]

    if nmin >0 :
        #restrict indices to n >=nmin and create an boolean index vector
        #also create a mask for the coefficients with n < nmin
        mask=[(n>=nmin) for m in range(nmax+1) for n in range(m,nmax+1)]
        nm_c=np.fromiter(nm_c, dtype='i,i')
        nm_s=np.fromiter(nm_s, dtype='i,i')
        nm_c=nm_c[mask]
        nm_s=nm_s[mask]
    else:
        #this mask will just select everything in the shtns SH vector
        mask=slice(None)


    return nm_c,nm_s,conv,mask


def lonlat_contigous_flags(dain):
    loninfo=find_lon(dain.coords)
    latinfo=find_lat(dain.coords)

    cont_stride_inbytes=dain.data.dtype.itemsize
    strides=dain.data.strides
    ilat=dain.dims.index(latinfo.var.name)
    ilon=dain.dims.index(loninfo.var.name)
    latstride=strides[ilat]
    lonstride=strides[ilon]
    
    
    if latstride == cont_stride_inbytes and lonstride == latstride*latinfo.var.size:
        return shtns.SHT_THETA_CONTIGUOUS,ilon,ilat
    elif lonstride == cont_stride_inbytes and latstride == lonstride*loninfo.var.size:
        return shtns.SHT_PHI_CONTIGUOUS,ilon,ilat
    else:
        return -1,ilon,ilat #no contiguous lon or lat

def south_pole_first(latinfo):
    if latinfo.direction == 'ascending':
        #note: different order of co-latitudes (south to north)
        grdflags=shtns.SHT_SOUTH_POLE_FIRST
    else:
        grdflags=0

    return grdflags

def infer_gtype(loninfo,latinfo):
    gtype=None
    grdflags=0
    if 'shxarray_gtype' in loninfo.var.attrs and 'shxarray_gtype' in latinfo.var.attrs:
        gtype=loninfo.var.attrs['shxarray_gtype']
    elif latinfo.var.dims[0] == loninfo.var.dims[0]:
        #ok set to point grid type
        gtype='point'
    elif latinfo.step  is not None and loninfo.step is not None:
        if loninfo.step != latinfo.step:
            raise RuntimeError("SHTns Backend: longitude and latitude do not having the same spacing")

        if not (latinfo.max-latinfo.min+latinfo.step) >= 180:
            raise RuntimeError("SHTns Backend: for regular grids, latitude does not cover 180 degrees")
        if not(loninfo.min == 0 and (loninfo.max-loninfo.min+loninfo.step) >=360):
            raise RuntimeError("SHTns Backend: for regular grids, longitude must start at degree 0 up to 360 (not included)")
        gtype='regular_lon0'

    if gtype == "point":
        grdflags=-1
    elif gtype == "regular_lon0":
        grdflags=shtns.sht_reg_fast
    elif gtype == "gauss":
        grdflags=shtns.sht_gauss
    else:
        raise RuntimeError("SHTns Backend: Only 'gauss','regular_lon0' and 'point' grid types are supported")
    return gtype,grdflags


def get_compatible_grid(dslonlat):
    loninfo=find_lon(dslonlat.coords)
    latinfo=find_lat(dslonlat.coords)

    grdflags=south_pole_first(latinfo)
    
    gtype,tmpflags=infer_gtype(loninfo,latinfo)
    grdflags|=tmpflags

    
    dslonlat=xr.Dataset(coords=dict(lon=loninfo.var,lat=latinfo.var))
    dslonlat.lat.attrs['shxarray_gtype']=gtype
    dslonlat.lon.attrs['shxarray_gtype']=gtype
    dslonlat.lat.attrs['shtns_gridflags']=grdflags
    dslonlat.lon.attrs['shtns_gridflags']=grdflags
    return dslonlat

class Analysis_shtns:
    shcreateflags=shtns.sht_orthonormal+shtns.SHT_NO_CS_PHASE
    def __init__(self,nmax):
        self.nmax=nmax
        self.shtns=shtns.sht(self.nmax, self.nmax, 1, norm=self.shcreateflags)

    def __call__(self,dain):
        loninfo=find_lon(dain.coords)
        latinfo=find_lat(dain.coords)

        grdflags=south_pole_first(latinfo)
    
        gtype,tmpflags=infer_gtype(loninfo,latinfo)
        if gtype == "point":
            raise RuntimeError("SHTns Backend: Analysis is not supported for point grids")
        
        grdflags|=tmpflags
        auxdims=[dim for dim in dain.dims if dim not in [loninfo.var.dims[0],latinfo.var.dims[0]]]

        if len(auxdims) > 1:
            dain=dain.stack(auxdim=auxdims)
            auxdim='auxdim'
        elif len(auxdims) == 0:
            auxdim="auxsingle"
            #expand for consistency with the below code
            dain=dain.expand_dims(auxdim)
        else:
            auxdim=auxdims[0]
        
        naux=dain.sizes[auxdim]

        if dain.dims.index(auxdim) == len(dain.dims)-1:
            dain=dain.T #transpose so that the auxiliary dimension will be the first
        elif dain.dims.index(auxdim) != 0:
            raise RuntimeError("auxiliary dimensions must be either first or last")
        
        cont_flag,ilon,ilat=lonlat_contigous_flags(dain)

        if cont_flag == -1:
            #We need to copy the data and can't work with the original data views
            usecopy=True
            if ilon < ilat:
                cont_flag=shtns.SHT_THETA_CONTIGUOUS
            else:
                cont_flag=shtns.SHT_PHI_CONTIGUOUS
        else: 
            usecopy=False

        
        grdflags|=cont_flag

        # #create an spherical harmonic output datarray of zeros
        daout=xr.DataArray.sh.zeros(self.nmax,auxcoords={auxdim:dain.coords[auxdim]})
        
        # #setup SHTns 
        # # sh = shtns.sht(nmax, nmax,1,shtns.sht_orthonormal)
        nlon=loninfo.var.size
        nlat=latinfo.var.size
        self.shtns.set_grid(nlat,nlon,grdflags)

        nm_c,nm_s,conv,mask=get_shtns_nm_order(self.nmax)
        cnm_comp=self.shtns.spec_array()
        for i in range(naux):
            if usecopy:
                grddata=dain.data[i,:,:].copy(order='C')
            else:
                grddata=dain.data[i,:,:]
            self.shtns.spat_to_SH(grddata,cnm_comp)
            cnm_comp/=conv
            #scale and asssign values to real spherical harmonics
            daout[i].loc[nm_s]=-np.imag(cnm_comp)
            #note: order matters, m=0 coefficients will be overwritten with correct ones below
            daout[i].loc[nm_c]=np.real(cnm_comp)
        

        #unstack the auxiliary dimensions if needed 
        if auxdim == 'auxdim':
            daout=daout.unstack(auxdim)
        elif auxdim == 'auxsingle':
            daout=daout.squeeze(auxdim,drop=True)
        return daout

class Synthesis_shtns:
    shcreateflags=shtns.sht_orthonormal+shtns.SHT_NO_CS_PHASE
    def __init__(self,dslonlat):
        #check for grid compatibility
        self._dsobj=get_compatible_grid(dslonlat)
        self._sh=None
        self.nmax=-1
    @property 
    def shtns(self):
        
        if self.nmax < 0:
            raise RuntimeError("SHTns Backend: nmax should be set before calling the synthesis method")

        if self._sh is None:
            self._sh=shtns.sht(self.nmax, self.nmax, norm=self.shcreateflags)
        
        return self._sh 

    def __call__(self,dain):
        self.nmax=dain.sh.nmax
        self.nmin=dain.sh.nmin
        nm_c,nm_s,conv,mask=get_shtns_nm_order(self.nmax,self.nmin)

        grdflags=self._dsobj.lat.attrs['shtns_gridflags']
        auxdims=[dim for dim in dain.dims if dim not in ['nm']]

        if len(auxdims) > 1:
            dain=dain.stack(auxdim=auxdims)
            auxdim='auxdim'
        elif len(auxdims) == 0:
            auxdim="auxsingle"
            #expand for consistency with the below code
            dain=dain.expand_dims(auxdim)
        else:
            auxdim=auxdims[0]
        
        naux=dain.sizes[auxdim]
        if dain.dims.index(auxdim) == len(dain.dims)-1:
            dain=dain.T #transpose so that the auxiliary dimension will be the first
        elif dain.dims.index(auxdim) != 0:
            raise RuntimeError("auxiliary dimension must be either first or last")

        if grdflags >= 0:
            nlon=self._dsobj.sizes['lon'] 
            nlat=self._dsobj.sizes['lat']
            #contiguous flags actually do not matter here
            #grdflags!=shtns.SHT_PHI_CONTIGUOUS
            #grdflags!=shtns.SHT_THETA_CONTIGUOUS
            self.shtns.set_grid(nlat,nlon,grdflags)


            #create a array holding the complex SH coefficients
            cnmcom=self.shtns.spec_array()
        
            #create an spherical harmonic output datarray of zeros
            coordsout={auxdim:dain.coords[auxdim],'lat':self._dsobj.coords['lat'],'lon':self._dsobj.coords['lon']}
            outshape=[naux,nlat,nlon]
            daout=xr.DataArray(np.zeros(outshape,order='C'),coords=coordsout,dims=[auxdim,'lat','lon'])
        
            for i in range(naux):
                #scale and assign values to real spherical harmonics
                cnmcom[mask]=dain[i].loc[nm_c].data-1j*dain[i].loc[nm_s].data
                cnmcom*=conv
                self.shtns.SH_to_spat(cnmcom,daout.data[i,:,:])
        else:
            nlonlat=self._dsobj.sizes['nlonlat']
            #create a array holding the complex SH coefficients
            cnmcom=self.shtns.spec_array()
            #create an spherical harmonic output datarray of zeros
            coordsout={auxdim:dain.coords[auxdim],'lat':("nlonlat",self._dsobj.coords['lat'].data),'lon':("nlonlat",self._dsobj.coords['lon'].data)}
            
            outshape=[naux,nlonlat]
            daout=xr.DataArray(np.zeros(outshape,order='C'),coords=coordsout,dims=[auxdim,'nlonlat'])
            for i in range(naux):
                #scale and assign values to real spherical harmonics
                cnmcom[mask]=dain[i].loc[nm_c].data-1j*dain[i].loc[nm_s].data
                cnmcom*=conv
                j=0
                for cost,phi in zip(np.sin(np.deg2rad(self._dsobj.lat.data)),np.deg2rad(self._dsobj.lon.data)):
                    daout.data[i,j]=self.shtns.SH_to_point(cnmcom,cost,phi)
                    j+=1
        
        if auxdim == 'auxdim':
            daout=daout.unstack(auxdim)
        elif auxdim == 'auxsingle':
            daout=daout.squeeze(auxdim,drop=True)
        daout.name=dain.name
        return daout 



class SHTnsBackend(SHComputeBackendBase):
    _credits="Used backend: SHTns High performance Spherical Harmonic Transform for numerical simulations (https://nschaeff.bitbucket.io/shtns)"
    
    def __init__(self):
        super().__init__()
        try:
            import shtns
        except ImportError:
            raise ImportError("SHTns Backend: shtns package is not installed. Please install it using 'pip install shtns")

    def synthesis(self, dain,dslonlat):
        syn=Synthesis_shtns(dslonlat)
        
        daout=syn(dain)
        
        daout.name=dain.name
        daout.attrs.update(get_cfglobal())
        daout.attrs['comments']=self._credits
        daout.attrs['history']="Synthesis operation"
        return daout
    
    def analysis(self,dain,nmax,**kwargs):
          

        ana=Analysis_shtns(nmax)

        daout=ana(dain)
        
        daout.name=dain.name
        daout.attrs.update(get_cfglobal())
        daout.attrs['comments']=self._credits
        daout.attrs['history']="Analysis operation"
        return daout

    
    def lonlat_grid(self,nmax=None,lon=None,lat=None,gtype="gauss"):
        """
        Create/check a lon-lat grid compatible with the SHTns backend
        Parameters
        ----------
        nmax : int, optional
            Maximum degree of the spherical harmonic expansion.
        lon : array-like, optional
            longitude coordinates
        lat : array-like, optional
            latitude coordinates
        type : str, optional default "gauss"
            Type of grid to be created. 'gauss','regular_lon0' or 'point' is supported.
        Returns
        -------
            An Dataset with lon, lat coordinates variables compatible with the SHTns backend

        """
        if gtype is None:
            #force to regular_lon0`
            gtype="regular_lon0"


        shtns_type=shtns.sht_quick_init
        if gtype == "gauss":
            shtns_type=shtns.sht_gauss
        elif gtype == "regular_lon0":
            shtns_type=shtns.sht_reg_fast
        elif gtype == "point":
            shtns_type=-1
        else:
            raise ValueError("SHTns Backend: Only 'gauss','regular' and 'point' grid types are supported")

        if nmax is not None:
            if lon is not None or lat is not None:
                raise ValueError("SHTns Backend: Either nmax or lon and lat should be specified")
            sh = shtns.sht(nmax, nmax, norm=shtns.sht_orthonormal+shtns.SHT_NO_CS_PHASE)

            nlat, nphi = sh.set_grid(flags=shtns_type)
            lon=np.linspace(0,360,nphi,endpoint=False)
            lat=np.rad2deg(np.asin(sh.cos_theta))
            dslonlat=xr.Dataset(coords=dict(lon=lon,lat=lat))
        elif lon is not None and lat is not None:
            if gtype == "point":
                if len(lon) != len(lat):
                    raise ValueError("SHTns Backend: lon and lat should have the same length fotr 'point' grids")
                dslonlat=xr.Dataset(coords=dict(lon=("nlonlat",lon),lat=("nlonlat",lat)))
            else:
                #we need to check whether the proposed grid is compatible with SHTns
                #rough estimate of nmax needed
                dres=min(np.median(np.diff(lon)),np.median(np.diff(lat)))
                nmax=int(180/dres)
                sh = shtns.sht(nmax, nmax, norm=shtns.sht_orthonormal+shtns.SHT_NO_CS_PHASE)
                sh.set_grid(len(lat),len(lon),shtns_type)
                dslonlat=xr.Dataset(coords=dict(lon=lon,lat=lat))

        dslonlat.lon.attrs.update(get_cfatts("longitude"))
        dslonlat.lat.attrs.update(get_cfatts("latitude"))
        if dslonlat.lon.min().item() < 0 and dslonlat.lon.max().item() <= 180:
            if gtype != 'point':
                dslonlat=change_central_longitude(dslonlat,central_longitude=180,resort= (gtype != 'point'))

        #add specific attributes for shtns
        dslonlat.lon.attrs['shtns_type']=shtns_type
        dslonlat.lat.attrs['shtns_type']=shtns_type
        dslonlat.lon.attrs["shxarray_gtype"]=gtype
        dslonlat.lat.attrs["shxarray_gtype"]=gtype
        dslonlat.attrs.update(get_cfglobal())
        dslonlat.attrs['comments']=self._credits
        return dslonlat




