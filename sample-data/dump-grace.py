#!/usr/bin/env python
# R. Rietbroek
# Extract a subset of GRACE data from a geoslurp database in a sqlite table and tar archive
# This can be used in the testing and tutorial examples

#Note to execute this script a populated database and a geoslurp client (https://github.com/strawpants/geoslurp) is required
from geoslurp.db.exporter import exportQuery
from geoslurp.db import geoslurpConnect
import os

qry2020="select * from gravity.gracecomb_l2_jpl_n60 WHERE date_part('year',time) = 2020"

outputfile="GRACEDataSample_2020.sql"

conn=geoslurpConnect(dbalias="marge")
localdataroot=conn.localdataroot
exportQuery(conn.dbeng.execute(qry2020),outputfile,layer="gracel2",packUriCols=["gsm","gac","gaa","gad"],striproot=os.path.join(localdataroot,"gravity"),localdataroot=localdataroot)

#Also add a static gravity field to subtract
staticqry="select * from gravity.icgem_static where uri LIKE '%%GGM05C%%'"
exportQuery(conn.dbeng.execute(staticqry),outputfile,layer="static",packUriCols=["uri"],striproot=os.path.join(localdataroot,"gravity"),localdataroot=localdataroot)
