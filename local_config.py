import os

dfm_lnx64="/opt/delft3dfm_2021.03/lnx64"                                                                                                              
os.environ['LD_LIBRARY_PATH']=os.path.join(dfm_lnx64,'lib')                                                                                           
  
class LocalConfig(object):
    num_procs=1
    dfm_bin_dir=os.path.join(dfm_lnx64,"bin")
    waq_proc_def=os.path.join(dfm_lnx64,"share/delft3d/proc_def.def")
    
