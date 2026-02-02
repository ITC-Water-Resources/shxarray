# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2025
#
import xarray as xr
from sparse import COO


def coo_serialize(ds):
    """
    Serializes an xarray with COO sparse matrix variables to one with separate variables. which can be writen to a file.
    Parameters
    ----------
    ds : xarray.Dataset or xarray.DataArray
        Input xarray object with sparse COO matrix variables.
    Returns
    -------
    xarray.Dataset
        Modified xarray object with separate variables for COO data and coordinates.
    """

    if type(ds) is not xr.Dataset:
        dsmod=ds.to_dataset(name='data',promote_attrs=True)
    else:
        dsmod=ds.copy()

    for vname,var in dsmod.variables.items():
        if hasattr(var.data,'format'):
            sparse_type=var.data.format
            if sparse_type == 'coo':
                #extract data and index arrays and register them as a set of new variables
                vname_co=f"{vname}_coo_co"
                vname_data=f"{vname}_coo_data"
                vname_dim=f"{vname}_coo_nnz"
                orig_dims=list(var.dims)
                attrs=var.attrs.copy()
                attrs.update(dict(coo_name=f"COO data array for variable {vname}",orig_name=vname,orig_dims=orig_dims))
                dsmod[vname_data]=([vname_dim],var.data.data,attrs)
                dsmod[vname_co]=(['coo_dim',vname_dim],var.data.coords,dict(long_name=f"COO coordinate array for variable {vname}",orig_name=vname,orig_dims=orig_dims))
                #remove the original variable
                dsmod=dsmod.drop_vars(vname)

            else:
                raise ValueError(f"Cannot currently handle sparse type {sparse_type} for variable {vname}")
    return dsmod



def coo_deserialize(ds):
    """
    Deserializes an xarray with separate COO sparse matrix variables to one with COO sparse matrix variables.
    Parameters
    ----------
    ds : xarray.Dataset
        Input xarray object with separate variables for COO data and coordinates.
    Returns
    -------
    xarray.Dataset
        Modified xarray object with sparse COO matrix variables.
    """
    dsmod=ds.copy()

    for vname,var in dsmod.variables.items():
        if vname.endswith('_coo_data'):
            #get the original variable name
            orig_name=var.attrs['orig_name']
            orig_dims=var.attrs['orig_dims']
            vname_co=f"{orig_name}_coo_co"
            vname_dim=f"{orig_name}_coo_nnz"
            shape=[dsmod.sizes[d] for d in orig_dims]
            #reconstruct sparse array
            dmat=COO.from_iter(zip(dsmod[vname_co].data.T,var.data),shape=shape)
            attrs={ky:val for ky,val in var.attrs.items() if ky not in ['orig_name','orig_dims']} 

            dsmod[orig_name]=(orig_dims,dmat,attrs)
            #cleanup
            dsmod=dsmod.drop_vars([vname,vname_co])
    return dsmod




