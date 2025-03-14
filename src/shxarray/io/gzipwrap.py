# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2024
#

from gzip import GzipFile
import os
from io import TextIOWrapper

try:
    from rapidgzip import RapidgzipFile
except:
    #failed to load (ok)
    RapidgzipFile=None

def gzip_open_r(filename,textmode=False):
    """
        GZip file reading wrapper leveraging parallel decompression speed when rapidgzip is installed on the system
    Parameters
    ----------
    filename : str
        Filename of the gzip archive to open
        
    textmode : 
        Whether to open in textmode (allows iterating over lines)
        
    """
    if RapidgzipFile: 
        gzfid=RapidgzipFile(filename,parallelization=os.cpu_count())

    else:
        gzfid=GzipFile(filename,'rb')


    if textmode:
        return TextIOWrapper(gzfid)
    else:
        return gzfid



