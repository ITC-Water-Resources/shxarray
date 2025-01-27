# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

import os
import yaml

def defaultcache(subdir=None):
    path=os.path.expanduser("~/.cache/shxarray_storage")
    if subdir is not None:
        path=os.path.join(path,subdir)
    os.makedirs(path,exist_ok=True)
    return path

def infofile():
    return os.path.join(defaultcache(),"shxarray_userinfo.yaml")


infodictglobal=None

def get_userinfo(infodict=infodictglobal):
    if not infodict:
        userinfo=infofile()
        if os.path.exists(userinfo):
            with open(userinfo) as f:
                infodict=yaml.safe_load(f)
        else:
            infodict={}
    return infodict

def set_userinfo(namecontact:str=None,institution:str=None,write=False):
    """
    Set  user information which will eb stored in xarray outputs
    Parameters
    ----------
    namecontact : str
        e.g. J. Doe (j.doe@email.com)
        
    institution : str
        Affiliation of the author as a string 
        
    write : bool
        Write the info to disk (make it persistent for next sessions)        

    """
    if namecontact is not None:
        infodict["contact"]=namecontact
    
    if institution is not None:
        infodict["institution"]=institution


    if write:
        userinfo=infofile()
        with open(userinfo,"w") as f:
            yaml.dump(infodict,f)



        
