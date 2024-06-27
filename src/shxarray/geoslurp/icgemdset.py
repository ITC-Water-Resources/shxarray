# This file is part of shxarray.
# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2024
# provides a dataset and table for static gravity fields from the icgem website



from geoslurp.dataset import DataSet
from shxarray.geoslurp.gravity import GravitySHTBase
from shxarray.geoslurp.icgem import  Crawler as IcgemCrawler
from shxarray.geoslurp.icgem import icgemMetaExtractor
from geoslurp.datapull.uri import findFiles
import re
from geoslurp.datapull import UriFile
import os

schema="shxarray"

class ICGEMstatic(DataSet):
    """Manages the static gravity fields which are hosted at http://icgem.gfz-potsdam.de/tom_longtime"""
    table=type("ICGEMstaticTable",(GravitySHTBase,), {})
    schema=schema
    stripuri=True
    def __init__(self, dbconn):
        super().__init__(dbconn)
        #initialize postgreslq table
        GravitySHTBase.metadata.create_all(self.db.dbeng, checkfirst=True)
        self.updated=[]
    
    def pull(self,pattern=None,list=False):
        """Pulls static gravity fields from the icgem website
        :param pattern: only download files whose name obeys this regular expression
        :param list (bool): only list available models"""
        self.updated=[]
        crwl=IcgemCrawler()
        if pattern:
            regex=re.compile(pattern)
        outdir=self.dataDir()
        if list:
            print("%12s %5s %4s"%("name","nmax", "year"))
        for uri in crwl.uris():
            if pattern:
                if not regex.search(uri.name):
                    continue
            if list:
                #only list available models
                print("%-12s %5d %4d"%(uri.name,uri.nmax,uri.lastmod.year))
            else:
                tmp,upd=uri.download(outdir,check=True, gzip=True)
                if upd:
                    self.updated.append(tmp)

    def register(self,pattern=None):
        """Register static gravity fields donwloaded in the data director
        :param pattern: only register files whose filename obeys this regular expression
        """
        if not pattern:
            pattern='.*\.gz'
        #create a list of files which need to be (re)registered
        if self.updated:
            files=self.updated
        else:
            files=[UriFile(file) for file in findFiles(self.dataDir(),pattern)]

        #loop over files
        for uri in files:
            urilike=os.path.basename(uri.url)

            if not self.uriNeedsUpdate(urilike,uri.lastmod):
                continue

            meta=icgemMetaExtractor(uri)
            self.addEntry(meta)

        self.updateInvent()


