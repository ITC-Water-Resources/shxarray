# This file is part of geoslurp.
# geoslurp is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.

# geoslurp is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with Frommle; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

# Author Roelof Rietbroek (r.rietbroek@utwente.nl), 2024

from geoslurp.view.viewBase import TView

schema='shxarray'

def subqry1(table,nmax,tol,apptypes):
    tname=f"{schema}.{table}"
    sbqry="SELECT gsm.tstart+(gsm.tend-gsm.tstart)/2 AS time, gsm.uri AS gsm"
    for apptype in apptypes:
        applo=apptype.lower()
        sbqry+=f", {applo}.uri AS {applo}"
    sbqry+=f" FROM (SELECT tstart,tend,uri from {tname} WHERE type ='GSM' AND nmax = {nmax})AS gsm"

    for apptype in apptypes:
        applo=apptype.lower()
        sbqry+=f" INNER JOIN {tname} {applo} ON ABS(EXTRACT(day FROM gsm.tstart-{applo}.tstart)) < {tol} AND  ABS(EXTRACT(day FROM gsm.tend-{applo}.tend)) < {tol} AND {applo}.type = '{apptype}'"
    return sbqry

def buildGSML2qry(gtable,gfotable,nmax,tol=8,apptypes=["GAA","GAB","GAC","GAD"]):
    qry="%s UNION %s ORDER BY time"%(subqry1(gtable,nmax,tol,apptypes),subqry1(gfotable,nmax,tol,apptypes))
    return qry

class GRACECOMB_L2_JPL_n96(TView):
    schema=schema
    qry=buildGSML2qry("grace_jpl_rl06","grace_fo_jpl_rl063",96) 

class GRACECOMB_L2_JPL_n60(TView):
    schema=schema
    qry=buildGSML2qry("grace_jpl_rl06","grace_fo_jpl_rl063",60) 


class GRACECOMB_L2_GFZ_n96(TView):
    schema=schema
    qry=buildGSML2qry("grace_gfz_rl06","grace_fo_gfz_rl063",96) 

class GRACECOMB_L2_GFZ_n60(TView):
    schema=schema
    qry=buildGSML2qry("grace_gfz_rl06","grace_fo_gfz_rl063",60) 


class GRACECOMB_L2_CSR_n96(TView):
    schema=schema
    qry=buildGSML2qry("grace_csr_rl06","grace_fo_csr_rl063",96,8,["GAC","GAD"]) 


class GRACECOMB_L2_CSR_n60(TView):
    schema=schema
    qry=buildGSML2qry("grace_csr_rl06","grace_fo_csr_rl063",60,8,["GAC","GAD"]) 

def getGRACEviews():
    return [ GRACECOMB_L2_JPL_n96, GRACECOMB_L2_JPL_n60, GRACECOMB_L2_GFZ_n96, GRACECOMB_L2_GFZ_n60, GRACECOMB_L2_CSR_n96, GRACECOMB_L2_CSR_n60]



