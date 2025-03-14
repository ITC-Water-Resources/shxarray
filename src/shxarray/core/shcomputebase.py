class SHComputeBackendBase:
    """
        Base class providing the calling interface for more compute intensive Spherical harmonic operations.
        The backend can be implemented in another module and registered as an entry_point in the pyproject.toml filei according to e.g.:

        [project.entry-points."shxarray.computebackends"]
        yourbackend = "yourpackage.module:SHComputeBackendclass"
    """
    _credits="Overwrite this statement in your derived class" 
    def synthesis(self,*argv):
        self._notImplemented("synthesis")
    
    def analysis(self,*argv):
        self._notImplemented("analysis")

    def lonlat_grid(self,nmax):
        self._notImplemented("lonlat_grid")

    def _notImplemented(self,methodname):
        raise RuntimeError(f"Method '{methodname}' is not implemented in backend {self.__class__.__name__}")
