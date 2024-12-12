import numpy as np
import xarray as xr
import shxarray
import gzip
from shxarray.core.sh_indexing import SHindexBase
from datetime import datetime,timedelta



def sinex_to_datetime(sn_date):
    year, doy, sec = [int(x) for x in sn_date.split(":")]

    if year < 50:
        year = year + 2000
    else:
        year = year + 1900
        
    return datetime(year, 1, 1) + timedelta(days = doy-1, seconds =sec)


shp = [SHindexBase.name]

def sinex_to_datetime(sn_date):
    year, doy, sec = [int(x) for x in sn_date.split(":")]

    if year < 50:
        year = year + 2000
    else:
        year = year + 1900
        
    return datetime(year, 1, 1) + timedelta(days = doy-1, seconds =sec)


def readblock_solest(fileobj):
    solution_estimate_data = []
    in_block = False
    for line in fileobj:
        line = line.decode('utf-8')
        if "+SOLUTION/ESTIMATE" in line:
            in_block = not in_block
            continue
        if in_block:
            if line.startswith("*"):
                continue
            if line.startswith("-"):
                break

            fields = line.split()
            #index = int(fields[0])
            type_ = str(fields[1])
            n = int(fields[2])
            if type_ == "SN":
                m = -int(fields[4])
            else:
                m = int(fields[4])
            time = sinex_to_datetime(fields[5])
            solution_value = float(fields[8])
            std_dev = float(fields[9])
            solution_estimate_data.append((n, m,solution_value, std_dev))
            
    sl_est = np.array(solution_estimate_data)
    nm = [(sl_est[i, 0], sl_est[i, 1]) for i in range(len(sl_est))]
   # time = sl_est[:, 2]
    solution_estimate = sl_est[:, 2].astype(float)
    std_dev = sl_est[:, 3].astype(float)
    return nm, time, solution_estimate, std_dev


def readblock_solapriori(fileobj):
    solution_apriori_data = []
    in_block = False
    for line in fileobj:
        line = line.decode('utf-8')
        if "+SOLUTION/APRIORI" in line:
            in_block = not in_block
            continue
        if in_block:
            if line.startswith("*"):
                continue
            if line.startswith("-"):
                break

            fields = line.split()
            #index = int(fields[0])
            #type_ = str(fields[1])
            #n = int(fields[2])
            #m = int(fields[4])
            apriori_value = float(fields[8])
            std_dev = float(fields[9])
            solution_apriori_data.append((apriori_value, std_dev))

    sl_apr = np.array(solution_apriori_data)
    #nm = [(sl_apr[i, 2], sl_apr[i, 3]) for i in range(len(sl_apr))]
    solution_apriori = sl_apr[:, 0].astype(float)
    std_dev = sl_apr[:, 1].astype(float)

    return solution_apriori, std_dev


def readblock_normal_equation_vector(fileobj):
    normal_equation_vector_data = []
    in_block = False
    for line in fileobj:
        line = line.decode('utf-8')
        if "+SOLUTION/NORMAL_EQUATION_VECTOR" in line:
            in_block = not in_block
            continue
        if in_block:
            if line.startswith("*"):
                continue
            if line.startswith("-"):
                break

            fields = line.split()
            #index = int(fields[0])
            #type_ = str(fields[1])
            #n = int(fields[2])
            #m = int(fields[4])
            solution_value = float(fields[8])
            normal_equation_vector_data.append(( solution_value))

    normal_equation_vector = np.array(normal_equation_vector_data).astype(float)
    #nm = [(normal_eq_vec[i, 2], normal_eq_vec[i, 3]) for i in range(len(normal_eq_vec))]
    #normal_equation_vector = normal_eq_vec[:, 0].astype(float)

    return  normal_equation_vector

def block_normal_eq_matrix(fileobj):
    normal_equation_matrix_u = []  
    in_block = False
    for line in fileobj:
        line  = line.decode('utf-8')
        if "+SOLUTION/NORMAL_EQUATION_MATRIX U" in line:
            in_block = not in_block
            continue

        if in_block:
            if line.startswith("*"):
                continue
            if line.startswith("-"):
                break

            fields = line.split()

            para1 = int(fields[0]) - 1
            para2 = int(fields[1]) - 1
            values = list(map(float, fields[2:4])) # create a list of floats using map object (iterator)
    
                    
            while len(normal_equation_matrix_u) <= para1:  
                normal_equation_matrix_u.append([])  ## Append the matrix size to include new row of data

            row = para1
            first_column = para2


                 # while len(normal_equation_matrix_u[row]) < first_column:
                    #     normal_equation_matrix_u[row].append(0.0)

            for i in range(len(values)):
                col_idx = first_column + i
                value = values[i]

                while len(normal_equation_matrix_u[row]) <= col_idx:
                    normal_equation_matrix_u[row].append(0.0)  ## add new columns to reach the desired col_idx


                        #if col_idx >= row:
                normal_equation_matrix_u[row][col_idx] = value  
                        # else:
                        #     normal_equation_matrix_u[row][col_idx] = 0.0


    max_row_length = max(len(row) for row in normal_equation_matrix_u) # Return the maximun row length of the matrix
    for row in normal_equation_matrix_u:
        row.extend([0.0] * (max_row_length - len(row)))
            
                
    normal_eq_matrix_u = np.array( normal_equation_matrix_u)
    num_rows = np.arange(normal_eq_matrix_u.shape[0])
    num_columns = np.arange(normal_eq_matrix_u.shape[1])
    
    symmetric_normal_eq_matrix = normal_eq_matrix_u + normal_eq_matrix_u.T - np.diag(np.diag(normal_eq_matrix_u)) # copy upper triangle to lower triangle in a python matrix
    return symmetric_normal_eq_matrix



def block_statistics(fileobj):

    inheader = False
    statistics = {}
    for ln in fileobj:
        ln = ln.decode('utf-8')
        if "+SOLUTION/STATISTICS" in ln:
            inheader = True
            continue
        if "-SOLUTION/STATISTICS" in ln:
            break

        spl = ln.split()
        if len(spl) >= 2:
            if "NUMBER OF DEGREES OF FREEDOM " in ln:
                statistics["degree_of_freedom"] = int(spl[-1])
            elif "NUMBER OF OBSERVATIONS" in ln:
                statistics["obs_num"] = int(spl[-1])
            elif "NUMBER OF UNKNOWNS" in ln:
                statistics["unknowns"] = int(spl[-1])
            elif "WEIGHTED SQUARE SUM OF O-C" in ln:
                statistics["wssoc"] = float(spl[-1])

    return statistics
            
            
    try:
        attr1["degree of freedom"] = int(spl[0])
        attr1["obs_num"] = int(spl[1])
        attr1["unknowns"] = int(spl[2])
        attr["wssoc"] = float(spl[5]) 
    except KeyError:
    #some values may not be present but that is ok
        pass 
        
    return attr1

def read_shsinex(fileobj, nmaxstop=96, logger=None):
    ds_sinex1=[]
    ds_sinex = [] 
    timee=[]
    needs_closing = False
    if isinstance(fileobj, str):
        needs_closing = True
        if fileobj.endswith('.gz'):
            fileobj = gzip.open(fileobj, 'rb')
        else:
            fileobj = open(fileobj, 'rb')
            
            
    #first read the icgem header 
    inheader=False
    hdr={}
    for ln in fileobj:
        ln = ln.decode('utf-8')
        if "+FILE/COMMENT" in ln:
            inheader=True
            continue
        if "-FILE/COMMENT" in ln:
            break
        
        spl=ln.split()
        if len(spl) == 2:
            #insert name value pairs in the hdr dict
            hdr[spl[0]]=spl[1]
    #extract relevant parameters from the header
    attr={}
    try:
        nmaxsupp=int(hdr["max_degree"])
        attr["nmaxfile"]=nmaxsupp
        if nmaxsupp < nmaxstop:
            attr["nmax"]=nmaxsupp
        else:
            attr["nmax"]=nmaxstop
        nmax=attr["nmax"]

        if nmax > nmaxsupp:
            logger.warning("Nmax ({nmax}) requested larger than supported, higher degree coefficients will be set to zero")


        if 'format' in hdr:
            attr["format"]=hdr['format']
        else:
            attr["format"]="sinex"
        
        if "norm" in hdr:
            attr["norm"]=hdr["norm"]
        
        attr["gm"]=float(hdr["earth_gravity_constant"])
        attr["re"]=float(hdr["radius"])
        attr["modelname"]=hdr["modelname"]
        attr["tidesystem"]=hdr["tide_system"]
        attr["errors"]=hdr["errors"]
        
        
    except KeyError:
    #some values may not be present but that is ok
        pass 
    
        


    attr1 = block_statistics(fileobj)
    attrr= {**attr, **attr1}
    nm,time, solution_estimate, solution_estimate_std_dev = readblock_solest(fileobj)
    solution_apriori, solution_apriori_std_dev = readblock_solapriori(fileobj)
    normal_equation_vector = readblock_normal_equation_vector(fileobj)
    normal_equation_matrix = block_normal_eq_matrix(fileobj)

    
    if needs_closing:
        fileobj.close()
    
#     if time:
#         shp=["time",SHindexBase.name]
#         coords={SHindexBase.name:SHindexBase.mi_fromtuples(nm),"time":time}
#         #also expand variables
#         solution_estimate=np.expand_dims(solution_estimate, axis=0)
#         solution_apriori_std_dev=np.expand_dims(solution_apriori_std_dev,axis=0)
#         solution_apriori=np.expand_dims(solution_apriori, axis=0)
#         normal_equation_vector=np.expand_dims(normal_equation_vector,axis=0)
        
#     else:
#         shp=[SHindexBase.name]
#         coords={SHindexBase.name:SHindexBase.mi_fromtuples(nm)}
#         solution_estimate=solution_estimate
#         solution_apriori_std_dev=solution_apriori_std_dev
#         solution_apriori=solution_apriori
#         normal_equation_vector=normal_equation_vector

    
    coords = {SHindexBase.name: SHindexBase.mi_fromtuples(nm)}

    timee.append(time)
    ds_sinex1 = xr.Dataset(data_vars = dict(sol_est = (shp, solution_estimate),sol_est_std_dev =(shp, solution_estimate_std_dev),sol_apriori = (shp, solution_apriori)
    ,normal_equation_vector = (shp, normal_equation_vector)), coords = coords,attrs=attrr)
    #.assign_coords(dict(time=("time", timee)))
    
    
    ds_sinex.append(ds_sinex1)
    ds = xr.concat(ds_sinex, dim="time").assign_coords(dict(time=("time", timee)))

    return ds
    
    
    