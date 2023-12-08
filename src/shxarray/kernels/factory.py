# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

from shxarray.kernels.gravfunctionals import Stokes2Eqh


class KernelFactory:
    """Provides an interface to generate multiple isotropic kernels"""
    @staticmethod
    def stokes2eqh(knLove):
        return Stokes2Eqh(knLove)

    @staticmethod
    def stokes2geoid():
        pass

    @staticmethod
    def gauss(halfwidthkm):
        pass

    @staticmethod
    def unitload(lon,lat):
        pass
