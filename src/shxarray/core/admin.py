# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

import os

def defaultcache(subdir=None):
    path=os.path.expanduser("~/.cache/shxarray_storage")
    if subdir is not None:
        path=os.path.join(path,subdir)
    os.makedirs(path,exist_ok=True)
    return path



