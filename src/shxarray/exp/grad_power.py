# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2025
#

from shxarray.core.logging import shxlogger
from shxarray.shlib import Pnm
from shxarray.core.sh_indexing import SHindexBase
import numpy as np
import xarray as xr
from scipy.sparse import block_diag,diags
from sparse import as_coo
from shxarray.earth.constants import a_earth
from shxarray.core.cf import get_cfatts,get_cfglobal 
def pnmbar_overlap_order(m1:int,m2:int,nmax:int,dx=1e-4):
    """
    Computes all non-zero combinations of Legendre overlap integrals for a given combination of orders and maximum degree
This routine uses a simple numerical Euler integration to approximate the overlap integral.
    I(n1, m1, n2,m2) = ∫ P_n1^m1(x) * P_n2^m2(x) dx from -1 to 1

    Parameters
    ----------
    m1 : int
        Order of the first associated Legendre function
        
    m2 : int
        Order of the second associated Legendre function
        
    nmax : int
        Maximum degree of the associated Legendre functions
        
    dx : float, optional
        Step size for numerical integrationiover the x axis, by default 1e-4
        
    
    Returns
    -------
    xr.DataArray 
        a checkerboard array spanning the valid degrees for the given orders m1 and m2. 

    """
    
    try:
        assert nmax >= m1 and nmax >= m2, "Maximum degree must be larger than the orders"
    except AssertionError as e:
        breakpoint()
    
    pnmbase=Pnm(nmax)
    
    #index relate to the requested order
    in1=[i for i,nm in enumerate(pnmbase.index()) if nm[1] == m1]
    in2=[i for i,nm in enumerate(pnmbase.index()) if nm[1] == m2]
    nm=np.array(pnmbase.index())

    #short cut diagonal matrix (when needed)
    if m1 == m2:
        diagval=2.0 if m1==0 else 4.0
        valout=np.diag(diagval*np.ones(len(in1)))

    else:

        x=np.arange(-1+dx/2,1,dx)
        pnmout1=np.zeros([len(in1),len(x)])
        pnmout2=np.zeros([len(in2),len(x)])
        for i,xx in enumerate(x):
            pnmtmp=pnmbase(xx)
            pnmout1[:,i]=pnmtmp[in1]
            pnmout2[:,i]=pnmtmp[in2]


        #simple Euler integration of all cross terms for simplicities sake
        valout=dx*pnmout1@pnmout2.T
        
        #note the constraints (see Mavromatis 1999) cause a checkerboard pattern where thje first entry is always non-zero
        # explicitly set to zeros
        valout[0::2,1::2]=0.0
        valout[1::2,0::2]=0.0


    #return the result as a xarray dataarray

    nm1=SHindexBase.mi_fromarrays(nm[in1].T,names=['n1','m1'])
    nm2=SHindexBase.mi_fromarrays(nm[in2].T,names=['n2','m2'])
    daout=xr.DataArray(valout,dims=['nm1','nm2'],coords=dict(nm1=(['nm1'],nm1),nm2=(['nm2'],nm2)))
    return daout
    

def getab(m,nmax):
    """
        Compute normalization adjustment factors
    """
    
    n=np.arange(m,nmax+1)
    sz=len(n)


    if m == 0:
        fac_a=1
    else:
        fac_a=2

    fac_b=fac_a/2


    if m != 1:
        fac_a/=2

    fac_2n=(2*n+1)/(2*n-1)
    #allocate space
    A=np.sqrt(fac_a*fac_2n*(n+m)*(n+m-1))
    B=np.sqrt(fac_b*fac_2n*(n-m)*(n-m-1))

    return A,B



def gradient_power_matrix(nmax,sat_radius=None,kn=None,direction='east_west'):
    """
        Compute a spherical harmonic matrix, Phi, which relates the power of an east-west gradient to a set of unknown Stokes Coefficients cnm:
        
        int (dphi/dlambda)^2 d omega = cnm' Phi cnm

    """

    #check for valid directions
    validdir=['east_west','north_south','horizontal']
    if direction not in validdir:
        raise RuntimeError(f"Invalid direction {direction}, must be one of {validdir}")



    #overall scale to apply (note that we leave out GM/A to be consistent with normalized unitless Stokes Coefficients.
    scale=np.pi/4

    #setup a degree dependent weighting
    if sat_radius is not None:
        rn=np.power(a_earth/sat_radius,np.arange(1,nmax+2))
    else:
        rn=np.ones([nmax+1])

    if kn is not None:
        #additionally apply degree dependent weighting
        if not hasattr(kn,'__getitem__'):
            raise RuntimeError("kn must be array-like")
            
        rn=np.array([rn[n]*kn[n] for n in np.arange(0,nmax+1)])


    if direction == 'north_south':
        north_south=True
    else:
        north_south=False

    if direction in ['north_south','horizontal']:
        #In order to compute north-south variation we need to first compute the diagonal matrix which represents the power of the horizontal gradient ( grad_h V)^2
        hor_grad_diag=np.diag([4*np.pi*n*(n+1)*np.power(rn[n],2) for n in range(0,nmax+1)])


    if direction == 'horizontal':
        # Shortcut : just return the diagonal horizontal gradient power matrix
        #generate the entire (diagonal) matrix for horizontal gradient power
        nm=SHindexBase.nm_mi(nmax,sort='m++/CS/n++') #make sure to use the same sorting as the other output matrices
        diagmat=diags([4*np.pi*n*(n+1)*np.power(rn[n],2) for n,m in nm])

        dictout={"readme":"Diagonal matrix describing the power of the horizontal gradient to Stokes coefficients"}
        dictout.update(get_cfglobal())
        shname=SHindexBase.name
        shname_t=shname+"_"
        daout=xr.DataArray(as_coo(diagmat),dims=[shname,shname_t],coords=dict(nm=([shname],nm,get_cfatts('nm')),nm_=([shname_t],SHindexBase.mi_toggle(nm),get_cfatts('nm'))),attrs=dictout)
        #we're done 
        return daout


    #compute weight matrix
    wmat=scale*np.outer(rn,rn)

    blocks=[]
    nm=[]

    #also include the order m=0 term (all zeros for E-W gradient)


    n1=np.arange(0,nmax+1)
    nm.extend([(n,0) for n in n1])
    if north_south:
        blocks.append(hor_grad_diag)
    else:
        blocks.append(np.zeros([len(n1),len(n1)]))


    for m1 in range(1,nmax+1):
        #build up the elements
        shxlogger.debug(f"Computing gradient power matrix block for order m={m1}") 
        a,b=getab(m1,nmax)
        
        n1=np.arange(m1,nmax+1)
        nm.extend([(n,m1) for n in n1])
        #also add sine component
        nm.extend([(n,-m1) for n in n1])

        #First diagonal term
        pnm1st=pnmbar_overlap_order(m1-1,m1-1,nmax).sel(m1=m1-1,m2=m1-1).sel(n1=n1-1,n2=n1-1)
        #shift degrees to represent output terms
        pnm1st['n1']=pnm1st.n1+1
        pnm1st['n2']=pnm1st.n2+1
        
        T=np.outer(a,a)*pnm1st

        #second diagonal term
        if m1+1 <= nmax:
            pnm2nd=pnmbar_overlap_order(m1+1,m1+1,nmax).sel(m1=m1+1,m2=m1+1).sel(n1=n1[2:]-1,n2=n1[2:]-1)
            pnm2nd['n1']=pnm2nd.n1+1
            pnm2nd['n2']=pnm2nd.n2+1

            #update relevant part of T
            T.loc[pnm2nd.n1,pnm2nd.n2]+=np.outer(b[2:],b[2:])*pnm2nd

            #cross terms (needed twice in transposed form)
            pnm3rd=pnmbar_overlap_order(m1-1,m1+1,nmax).sel(m1=m1-1,m2=m1+1).sel(n1=n1-1,n2=n1[2:]-1)
            pnm3rd['n1']=pnm3rd.n1+1
            pnm3rd['n2']=pnm3rd.n2+1
            T2=np.outer(a,b[2:])*pnm3rd
        
            #update
            T.loc[T2.n1,T2.n2]+=T2.data
            #also update with the transpose
            T.loc[T2.n2,T2.n1]+=T2.data.T
        
        #apply scale and possible weighting
        T=T*wmat[m1:nmax+1,m1:nmax+1]
        if north_south:
            #take the complement relative to the diagonal horizontal gradient
            T=hor_grad_diag[m1:,m1:]-T

        blocks.append(T.data)
        #also add sine block
        blocks.append(T.data)


    mshi=SHindexBase.mi_fromtuples(nm)
    dictout={"readme":"Block diagonal matrix describing the power of the east-west gradient to Stokes coefficients"}
    dictout.update(get_cfglobal())
    shname=SHindexBase.name
    shname_t=shname+"_"
    sparseblk=as_coo(block_diag(blocks,format='coo')) 
    # dsout=xr.Dataset(dict(mat=([shname,shname_t],sparseblk)),coords=dict(nm=([shname],mshi),nm_=([shname_t],SHindexBase.mi_toggle(mshi))),attrs=dictout)
    daout=xr.DataArray(sparseblk,dims=[shname,shname_t],coords=dict(nm=([shname],mshi,get_cfatts('nm')),nm_=([shname_t],SHindexBase.mi_toggle(mshi),get_cfatts('nm'))),attrs=dictout)
    return daout


def gradient_power_degvar(mat,cnm,average=False):
    """
        Compute the degree wise power in when propagating cnm' mat cnm as a  function of the degree only
    """
    nmax=mat.sh.nmax
    nmin=cnm.sh.nmin

    #expand the matrix to dense matrix for speed
    if hasattr(mat.data,'todense'):
        mat.data=mat.data.todense()
    nv=np.arange(nmin,nmax+1)
    # dvar=np.zeros([len(nv)])
    
    #create a base DataArray to fill
    shname=SHindexBase.name
    coout={dim:cnm[dim] for dim in cnm.dims if dim != shname}
    coout={'n':(['n'],nv)}
    dims=[dim for dim in cnm.dims if dim != shname]
    dims.append('n')
    shape=[cnm.sizes[dim] for dim in cnm.dims if dim != shname]
    shape.append(len(nv))

    dvar=xr.DataArray(np.zeros(shape),dims=dims,coords=coout,attrs={"long_name":"Degree wise gradient power","units":"-/m2"})
    for n in nv:
        submat=mat.sel(n=n,n_=n)
        subcnm=cnm.sel(n=n)
        if average:
            scale=1/(2*n+1)
        else:
            scale=1.0
        dvar.loc[dict(n=n)]=scale*(submat.dot(subcnm,dim='m').rename(m_='m')).dot(subcnm,dim='m')

    #turn into xarray
    # dvar=xr.DataArray(dvar,dims=['n'],coords=dict(n=(['n'],nv)),attrs={"long_name":"Degree wise gradient power","units":"-/m2"})
    return dvar

