# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2025
#

def leakage_vishwakarma2016(datws, basinavOp, filterOp):
    """
    Leakage correction as described in Vishwakarma et al. (2016) for the given Total water storage dataset and basins.
    
    Parameters
    ----------
    datws : xarray.Dataset
        Dataset containing the terrestrial water storage data.
    basinavOp :
        Basin averages operator must accept an array of the form datws in its __call__ operator
    filterOp : xarray.DataArray
        Filter operator which can be applied to the data or chained with the basin averages operator itself
    
    Returns
    -------
    xarray.Dataset
        Dataset containing basin averaged leakage corrected data.
    """

     
    pass




