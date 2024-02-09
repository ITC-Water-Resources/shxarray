class SHComputeBackendBase:
    """base class providing the calling interface for more compute intensive Spherical harmonic operations"""

    def synthesis(self,*argv):
        self._notImplemented("synthesis")
    
    def analysis(self,*argv):
        self._notImplemented("Analysis")

    def _notImplemented(self,methodname):
        raise RuntimeError(f"Method '{methodname}' is not implemented in backend {self.__class__.__name__}")
