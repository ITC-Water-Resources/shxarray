import numpy as np
import xarray as xr
import shxarray
import os
import matplotlib.pyplot as plt
import gzip
from shxarray.core.sh_indexing import SHindexBase
from datetime import datetime,timedelta

class Sinexread:
    def __init__(self):
        #Initialize an Sinexread object with four attributes
        self.solution_estimate = []
        self.solution_apriori = []
        self.normal_equation_vector = []
        self.normal_equation_matrix = []
        

    def read_file(self, fileobj):
        needsClosing=False
        if type(fileobj) == str:
            needsClosing=True
            if fileobj.endswith('.gz'):
                fileobj=gzip.open(fileobj,'rb')
            else:
                fileobj=open(fileobj,'rb')


        shp = [SHindexBase.name]  ## check if all variables have the same dimension/ define the dimensions of the data variables
 
        solution_estimate_data = []
        in_block = False
        for line in fileobj:
            line  = line.decode('utf-8') # Decode bytes to string
            if "+SOLUTION/ESTIMATE" in line:
                in_block = not in_block
                continue
            if in_block:
                if line.startswith("*"):
                    continue
                if line.startswith("-"):
                    break

                fields = line.split()
                index = int(fields[0])
                type_ = str(fields[1])
                n = int(fields[2])
                solution_value = float(fields[8])
                std_dev = float(fields[9])
                if type_ == "SN":
                    m = -int(fields[4])
                else:
                    m = int(fields[4])
                    #nm = (n,m)
                    
                    
                solution_estimate_data.append((index, type_, n, m, solution_value, std_dev))


        sl_est = np.array(solution_estimate_data)
            #nm = sl_est[:,0] 
        index = sl_est[:,0].astype(int)
        type_ = sl_est[:,1]
        n = sl_est[:,2].astype(int)
        m = sl_est[:,3].astype(int)
        nm = [(sl_est[i, 2],sl_est[i, 3]) for i in range(len(sl_est))]
        solution_estimate = sl_est[:,4].astype(float)
        std_dev = sl_est[:,5].astype(float)
    
        coords = {SHindexBase.name: SHindexBase.mi_fromtuples(nm)}
        ds_solution_estimate = xr.Dataset(data_vars = dict(solution_estimate = (shp, solution_estimate), std_dev = (shp, std_dev)), coords = coords)          
        self.solution_estimate =  ds_solution_estimate           

        solution_apriori_data = []
        nm = []; type_ = []; n = []; m = []; index =[]; coords =[]; std_dev =[]
        in_block = False
        for line in fileobj:
            line  = line.decode('utf-8')
            if "+SOLUTION/APRIORI" in line:
                in_block = not in_block
                continue
            if in_block:
                if line.startswith("*"):
                    continue
                if line.startswith("-"):
                    break

                fields = line.split()

                index = int(fields[0])
                type_ = str(fields[1])
                n = int(fields[2])
                apriori_value = float(fields[8])
                std_dev = float(fields[9])
                if type_ == "SN":
                    m = -int(fields[4])
                else:
                    m = int(fields[4])

                    #nm = (n,m)
                solution_apriori_data.append((index, type_, n, m, apriori_value, std_dev))
                    #self.solution_apriori.append((index, type_, n, m, apriori_value, std_dev))
                    
            
            
            
        sl_apr = np.array(solution_apriori_data)
        index = sl_apr[:,0].astype(int)
        type_ = sl_apr[:,1]
        n = sl_apr[:,2].astype(int)
        m = sl_apr[:,3].astype(int)
        nm = [(sl_apr[i, 2],sl_apr[i, 3]) for i in range(len(sl_est))]
        solution_apriori = sl_apr[:,4].astype(float)
        std_dev = sl_apr[:,5].astype(float)
    
        coords = {SHindexBase.name: SHindexBase.mi_fromtuples(nm)}
        ds_solution_apriori = xr.Dataset(data_vars = dict(solution_apriori = (shp, solution_apriori), std_dev = (shp,std_dev)), coords = coords)
        self.solution_apriori = ds_solution_apriori
            
        normal_equation_vector_data = []
        nm = []; type_ = []; n = []; m = []; index =[]; coords =[]; std_dev =[]
        in_block = False
        for line in fileobj:
            line  = line.decode('utf-8')
            if "+SOLUTION/NORMAL_EQUATION_VECTOR" in line:
                in_block = not in_block  
                continue

            if in_block:
                if line.startswith("*"):
                    continue
                if line.startswith("-"): 
                    break

                fields = line.split()

                index = int(fields[0])
                type_ = str(fields[1])
                n = int(fields[2])
                solution_value = float(fields[8])
                if type_ == "SN":
                    m = -int(fields[4])
                else:
                    m = int(fields[4])

                    #nm = (n,m)
                normal_equation_vector_data.append((index, type_, n, m, solution_value))

                    
        normal_eq_vec = np.array(normal_equation_vector_data)
        index = normal_eq_vec[:,0].astype(int)
        type_ = normal_eq_vec[:,1]
        n = normal_eq_vec[:,2].astype(int)
        m = normal_eq_vec[:,3].astype(int)
        nm =[(normal_eq_vec[i,2], normal_eq_vec[i,3]) for i in range(len(normal_eq_vec))]
        normal_equation_vector = normal_eq_vec[:,4].astype(float)
     
        coords = {SHindexBase.name: SHindexBase.mi_fromtuples(nm)}
        ds_normal_equation_vector = xr.Dataset(data_vars = dict(normal_equation_vector = (shp, normal_equation_vector)), coords = coords)
    
        self.normal_equation_vector = ds_normal_equation_vector
    
                    
        in_block = False
        normal_equation_matrix_u = []  

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


            #self.normal_equation_matrix = normal_equation_matrix_u   
            
            
            
        normal_eq_matrix_u = np.array( normal_equation_matrix_u)
        num_rows = np.arange(normal_eq_matrix_u.shape[0])
        num_columns = np.arange(normal_eq_matrix_u.shape[1])
    
        symmetric_normal_eq_matrix = normal_eq_matrix_u + normal_eq_matrix_u.T - np.diag(np.diag(normal_eq_matrix_u)) # copy upper triangle to lower triangle in a python matrix
    

            #nm = [(i, j) for i in n for j in m]
            #coords = {SHindexBase.name: SHindexBase.mi_fromtuples(nm)}
        ds_symmetric_normal_equation_matrix = xr.Dataset(data_vars = dict(symmetric_normal_eq_matrix = (["rows","columns"], symmetric_normal_eq_matrix)), coords = dict(rows = ("rows",num_rows), columns = ("columns", num_columns)))
        self.normal_equation_matrix = ds_symmetric_normal_equation_matrix
        
        
#         if needsClosing:
#             fileobj.close()

#         if time:
#             shp=["time",SHindexBase.name]
#             coords={SHindexBase.name:SHindexBase.mi_fromtuples(nm),"time":time}
#             cnm=np.expand_dims(cnm[0:ncount], axis=0)
#             sigcnm=np.expand_dims(sigcnm[0:ncount],axis=0)

    

        
        return self
    
    