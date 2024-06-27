# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2024
#


from geoslurp.dataset import DataSet
from geoslurp.datapull.ftp import Uri as ftp
from geoslurp.datapull.http import Uri as http
from datetime import datetime,timedelta
from shxarray.core.time import decyear2dt
import os
import re
from shxarray.geoslurp.gravity import GravitySHTBase,GravitySHinDBTBase,Trig,JSONSHArchive
import xarray as xr

schema="shxarray"

# class geocenter_Rietbroeketal2016upd(DataSet):
    # fout="Geocenter_dec2017.tgz"
    # sqrt3timesRE=11047256.4063275
    # schema=schema
    # table=type("geocenter_Rietbroeketal2016updTable", (GravitySHTBase,), {})
    # def __init__(self,dbconn):
        # super().__init__(dbconn)
        # # super().__init__(direc=direc,uri='https://wobbly.earth/data/Geocenter_dec2017.tgz',order=['c10','c11','s11'],lastupdate=datetime(2018,10,16))
    
    # def pull(self):
        # """Pulls the geocenter ascii files in the cache"""
        
        # uri=http("https://wobbly.earth/data/"+self.fout,lastmod=datetime(2018,10,16)).download(self.cacheDir(),check=True)


    # def register(self):
        # self.truncateTable()
        # #set general settings
        # self._dbinvent.data={"citation":"Rietbroek, R., Brunnabend, S.-E., Kusche, J., Schröter, J., Dahle, C., 2016. " \
                                      # "Revisiting the Contemporary Sea Level Budget on Global and Regional Scales. " \
                                      # "Proceedings of the National Academy of Sciences 201519132. " \
                                      # "https://doi.org/10.1073/pnas.1519132113"}
        
        # with tarfile.open(os.path.join(self.cacheDir(),self.fout),'r:gz') as tf:

            # metacomb=[]

            # files=['Geocenter/GeocentCM-CF_Green.txt',
                    # 'Geocenter/GeocentCM-CF_Antarctica.txt',
                   # 'Geocenter/GeocentCM-CF_Hydrology.txt',
                   # 'Geocenter/GeocentCM-CF_LandGlaciers.txt',
                   # 'Geocenter/GeocentCM-CF_GIA.txt',
                   # 'Geocenter/GeocentCM-CF_TotSurfload.txt']

            # order=[(1,1,Trig.c),(1,1,Trig.s),(1,0,Trig.c)]
            # lastupdate=datetime.now()
            # for file in files:
                # #get files
                # with tf.extractfile(file) as fid:
                    # for ln in fid:
                        # shar=JSONSHArchive(1)
                        # lnspl=ln.decode('utf-8').split()
                        # tcent=decyear2dt(float(lnspl[0]))
                        # tstart,tend=dt2monthlyinterval(tcent)
                        
                        # meta={"type":file.split('_')[-1][:-4],"time":tcent,"tstart":tstart,"tend":tend,"lastupdate":lastupdate,"nmax":1,"omax":1,"origin":"CF","format":"JSONB","uri":"self:data","gm":0.3986004415e+15,"re":0.6378136460e+07
# }
                        
                        # for el,val in zip(order,lnspl[1:4]):
                            # # import pdb;pdb.set_trace()
                            # shar["cnm"][shar.idx(el)]=float(val)/self.sqrt3timesRE

                        # #also add sigmas 
                        # for el,val in zip(order,lnspl[4:7]):
                            # shar["sigcnm"][shar.idx(el)]=float(val)/self.sqrt3timesRE
                        # meta["data"]=shar.dict
                        # self.addEntry(meta)
            # self.updateInvent()


def parseGSMDate(dtstr):
    """Parse datestr as found in GSM files (yyyymmdd.00000)"""
    return datetime(int(dtstr[0:4]),int(dtstr[4:6]),int(dtstr[6:8]))  

class geocenter_GRCRL06_TN13(DataSet):
    schema=schema
    rooturl="http://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/grace/open/docs/"
    def __init__(self,dbconn):
        self.table=type(self.__class__.__name__.lower().replace('-',"_")+"Table", (GravitySHinDBTBase,), {})
        super().__init__(dbconn)
    
    def pull(self):
        """Pulls the geocenter ascii files in the cache"""
        uri=http(self.rooturl+self.fout,lastmod=datetime(2019,12,1)).download(self.cacheDir(),check=True)


    def register(self):
        self.truncateTable()
        #set general settings
        self._dbinvent.data={"citation":"GRACE technical note 13"}
        lastupdate=datetime.now()
        with open(os.path.join(self.cacheDir(),self.fout),'r') as fid:
            #skip header
            for ln in fid:
                if re.search("end of header",ln):
                    break
            
            nv=[]
            mv=[]
            tv=[]
            cnmv=[]
            sigcnmv=[]
            #loop over entry lines 
            for ln in fid:
                
                tp,n,m,cnm,snm,sigcnm,sigsnm,ts,te=ln.split()
                
                #Append cosine coefficients
                n=int(n)
                m=int(m)
                nv.append(n)
                mv.append(m)
                tv.append(0)
                cnmv.append(float(cnm))
                sigcnmv.append(float(sigcnm))

                #append sine coefficients
                if m > 0:
                    nv.append(n)
                    mv.append(m)
                    tv.append(1)
                    cnmv.append(float(snm))
                    sigcnmv.append(float(sigsnm))

                if m == 1: 
                    #register the accumulated entry
                    tstart=parseGSMDate(ts)
                    tend=parseGSMDate(te)
                    #get the central time
                    tcent=tstart+(tend-tstart)/2
                    #snap the central epoch to the 15th of the month of the central time
                    # tcent=datetime(tstart.year,tstart.month,15)
                
                    meta={"type":"GSM","time":tcent,"tstart":tstart,"tend":tend,"lastupdate":lastupdate,"nmax":1,"omax":1,"origin":"CF","format":"JSONB","gm":0.3986004415e+15,"re":0.6378136460e+07}
                    meta["data"]=xr.Dataset(data_vars=dict(cnm=(["shg"],cnmv),sigcnm=(["shg"],sigcnmv)),coords=dict(n=(["shg"],nv),m=(["shg"],mv),t=(["shg"],tv)))
                    
                    self.addEntry(meta)
                    nv=[]
                    mv=[]
                    tv=[]
                    cnmv=[]
                    sigcnmv=[]

            self.updateInvent()
        





class GeocenterRIESCFCM(DataSet):
    fout30="GCN_L1_L2_30d_CF-CM.txt"
    fout60="GCN_L1_L2_60d_CF-CM.txt"
    #note also embed mm to meter conversion in here (e3)
    sqrt3timesRE=11047256.23312e3
    schema=schema
    table=type("geocenter_ries_cfcmTable", (GravitySHTBase,), {})
    def __init__(self,dbconn):
        super().__init__(dbconn)
        # super().__init__(direc=direc,uri='https://wobbly.earth/data/Geocenter_dec2017.tgz',order=['c10','c11','s11'],lastupdate=datetime(2018,10,16))
    
    def pull(self):
        """Pulls the geocenter ascii files in the cache"""
        
        uri=http("http://download.csr.utexas.edu/pub/slr/geocenter/"+self.fout30).download(self.cacheDir())
        uri=http("http://download.csr.utexas.edu/pub/slr/geocenter/"+self.fout60).download(self.cacheDir())


    def register(self):
        self.truncateTable()
        #set general settings
        self._dbinvent.data={"citation":"Ries, J.C., 2016. Reconciling estimates of annual geocenter motion from space geodesy, in: Proceedings of the 20th International Workshop on Laser Ranging, Potsdam, Germany. pp. 10–14."}
        
            
        self.extractSLR(os.path.join(self.cacheDir(),self.fout30))
        self.extractSLR(os.path.join(self.cacheDir(),self.fout60))
        self.updateInvent()

    def extractSLR(self,filen):
            lastupdate=datetime.now()
            order=[(1,1,Trig.c),(1,1,Trig.s),(1,0,Trig.c)]
            if re.search('30d',filen):
                dt=timedelta(days=15)
            else:
                dt=timedelta(days=30)

            with open(filen,'r') as fid:
                #skip header
                fid.readline()
                for ln in fid:
                    shar=JSONSHArchive(1)
                    lnspl=ln.split()
                    #note little hack as the 60d file version has an empty line at the back
                    if len(lnspl) == 0:
                        break
                    tcent=decyear2dt(float(lnspl[0]))
                    tstart=tcent-dt
                    tend=tcent+dt
                    
                    meta={"type":"GSM"+os.path.basename(filen)[10:13],"time":tcent,"tstart":tstart,"tend":tend,"lastupdate":lastupdate,"nmax":1,"omax":1,"origin":"CF","format":"JSONB","uri":"self:data","gm":0.3986004415e+15,"re":0.6378136460e+07
}
                    
                    for el,val in zip(order,lnspl[1:4]):
                        shar["cnm"][shar.idx(el)]=float(val)/self.sqrt3timesRE

                    #also add sigmas 
                    for el,val in zip(order,lnspl[4:7]):
                        shar["sigcnm"][shar.idx(el)]=float(val)/self.sqrt3timesRE
                    meta["data"]=shar.dict
                    self.addEntry(meta)


class  TN14SLRGSFC(DataSet):
    schema=schema
    rooturl="ftp://isdcftp.gfz-potsdam.de/grace-fo/DOCUMENTS/TECHNICAL_NOTES/"
    fout="TN-14_C30_C20_SLR_GSFC.txt"
    def __init__(self,dbconn):
        self.table=type(self.__class__.__name__.lower().replace('-',"_")+"Table", (GravitySHinDBTBase,), {})
        super().__init__(dbconn)
    
    def pull(self):
        """Pulls the C20 Technical note 14"""
        uri=ftp(self.rooturl+self.fout,lastmod=datetime(2019,12,1)).download(self.cacheDir(),check=True)


    def register(self):
        self.truncateTable()
        #set general settings
        self._dbinvent.data={"citation":"GRACE technical note 14"}
        lastupdate=datetime.now()
        mjd00=datetime(1858,11,17)
        nmax=2
        omax=0
        with open(os.path.join(self.cacheDir(),self.fout),'r') as fid:
            #skip header
            for ln in fid:
                if re.search("^Product:",ln):
                    break
            
            #loop over entry lines 
            for ln in fid:
                
                mjd0,decy0,c20,dc20,sigc20,c30,dc30,sigc30,mjd1,decyr1=ln.split()
                
                nv=[]
                mv=[]
                tv=[]
                cnmv=[]
                dcnmv=[]
                sigcnmv=[]
                
                #Append c20 coefficients
                nv.append(2)
                mv.append(0)
                tv.append(0)
                cnmv.append(float(c20))
                dcnmv.append(float(dc20)*1e-10)
                sigcnmv.append(float(sigc20))
                if c30 != "NaN":
                    nmax=3
                    nv.append(3)
                    mv.append(0)
                    tv.append(0)
                    cnmv.append(float(c30))
                    dcnmv.append(float(dc30)*1e-10)
                    sigcnmv.append(float(sigc30))

                #register the accumulated entry
                tstart=mjd00+timedelta(days=float(mjd0))
                tend=mjd00+timedelta(days=float(mjd1))
                tcent=tstart+(tend-tstart)/2
            
                meta={"type":"GSM","time":tcent,"tstart":tstart,"tend":tend,"lastupdate":lastupdate,"nmax":nmax,"omax":omax,"format":"JSONB","gm":0.3986004415e+15,"re":0.6378136460e+07}
                meta["data"]=xr.Dataset(data_vars=dict(cnm=(["shg"],cnmv),dcnm=(["shg"],dcnmv),sigcnm=(["shg"],sigcnmv)),coords=dict(n=(["shg"],nv),m=(["shg"],mv),t=(["shg"],tv)))
                
                self.addEntry(meta)
            self.updateInvent()

def getDeg1n2corrDsets(conf):
    out=[]
    for center in ["CSR", "GFZ", "JPL"]:
        out.append(type("geocenter_"+center+"RL06_TN13",(geocenter_GRCRL06_TN13,),{"fout":"TN-13_GEOC_"+center+"_RL06.txt"}))
    out.append(GeocenterRIESCFCM)
    out.append(TN14SLRGSFC)
    return out

# class Sun2017Comb(LowdegreeSource):
    # def __init__(self,direc,uri=None):
        # if not uri:
            # uri='https://d1rkab7tlqy5f1.cloudfront.net/CiTG/Over%20faculteit/Afdelingen/' \
                # 'Geoscience%20%26%20Remote%20sensing/Research/Gravity/models%20and%20data/3_Deg1_C20_CMB.txt'
        # lastupdate=datetime(2018,1,1)
        # super().__init__(direc=direc,uri=uri,lastupdate=lastupdate)
    # def extract(self):
        # singlemeta=self.meta

        # citation="Reference: Yu Sun, Pavel Ditmar, Riccardo Riva (2017), Statistically optimal" \
                                       # " estimation of degree-1 and C20 coefficients based on GRACE data and an " \
                                       # "ocean bottom pressure model Geophysical Journal International, 210(3), " \
                                       # "1305-1322, 2017. doi:10.1093/gji/ggx241."


        # valorder=['c10','c11','s11']
        # covorder=[('c10','c10'),('c10','c11'),('c10','s11'),('c10','c20'),('c11','c11'),('c11','s11'),('c11','c20'),('s11','s11'),('s11','c20'),('c20','c20')]

        # singlemeta['data']["citation"]=citation
        # with open(os.path.join(self.direc,self.fout),'rt') as fid:
            # dataregex=re.compile('^ +[0-9]')
            # time=[]
            # for ln in fid:
                # if dataregex.match(ln):
                    # lnspl=ln.split()
                    # time.append(decyear2dt(float(lnspl[0])))
                    # d12=[0]*4
                    # for el,val in zip(valorder,lnspl[1:4]):
                        # d12[self.dindex(el)]=val
                    # singlemeta['data']['d12'].append(d12)


                    # covUpper=[0]*10
                    # for (el1,el2),val in zip(covorder,lnspl[5:]):
                        # covUpper[self.covindex(el1,el2)]=val
                    # singlemeta['data']['covUpper'].append(covUpper)

        # singlemeta["data"]["time"]=[dt.isoformat() for dt in time]
        # singlemeta['tstart']=min(time)
        # singlemeta['tend']=max(time)
        # singlemeta['lastupdate']=datetime.now()

        # return [singlemeta]

# class Sun2017Comb_GIArestored(Sun2017Comb):
        # def __init__(self,direc):
            # uri="https://d1rkab7tlqy5f1.cloudfront.net/CiTG/Over%20faculteit/Afdelingen/" \
                     # "Geoscience%20%26%20Remote%20sensing/Research/Gravity/models%20and%20data/" \
                     # "4_Deg1_C20_CMB_GIA_restored.txt"
            # super().__init__(direc,uri)

# class GeocTable(GeocTBase):
    # """Defines the Geocenter motion table"""
    # __tablename__='deg1n2'
    # id=Column(Integer,primary_key=True)
    # name=Column(String,unique=True)
    # lastupdate=Column(TIMESTAMP)
    # tstart=Column(TIMESTAMP)
    # tend=Column(TIMESTAMP)
    # origin=Column(String)
    # data=Column(JSONB)

    # "ftp://ftp.csr.utexas.edu/pub/slr/geocenter/"
    #     # {'name':'Rietbroeketal2016updated','uri':'https://wobbly.earth/data/Geocenter_dec2017.tgz','lastupdate':datetime(2018,10,16)},
    #     # {'name':'SwensonWahr2008','uri':'ftp://podaac.jpl.nasa.gov/allData/tellus/L2/degree_1/deg1_coef.txt','lastupdate':datetime(2018,10,16)},
    #     # {'name':'Sun2017Comb','uri':'https://d1rkab7tlqy5f1.cloudfront.net/CiTG/Over%20faculteit/Afdelingen/Geoscience%20%26%20Remote%20sensing/Research/Gravity/models%20and%20data/3_Deg1_C20_CMB.txt'
    #     #                             ,'lastupdate':datetime(2018,10,16)},
    # ]


