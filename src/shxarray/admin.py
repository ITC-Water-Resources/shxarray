# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

import os

def defaultcache():
    path=os.path.expanduser("~/.shxarray_storage")
    os.makedirs(path,exist_ok=True)
    return path



