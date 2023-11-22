class SHComputeBackendBase:
    """base class providing the calling interface for more compute intensive Spherical harmonic operations"""

    def analysis(self,dsobj,lon,lat,forcegrid=False):
        self._notImplemented("analysis")

    def _notImplemented(self,methodname):
        raise RuntimeError(f"Method '{methodname}' is not implemented in backend {self.__class__.__name__}")
