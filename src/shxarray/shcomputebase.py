class SHComputeBackendBase:
    """base class providing the calling interface for more compute intensive Spherical harmonic operations"""

    def analysis(self,dain,lon,lat):
        self._notImplemented("Analysis")

    def _notImplemented(self,methodname):
        raise RuntimeError(f"Method '{methodname}' is not implemented in backend {self.__class__.__name__}")
