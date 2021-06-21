from __future__ import print_function
from root_pandas import read_root

path = '/afs/cern.ch/user/a/adsalvad/work/Trigger/AnalysisTrigger/'
datafile='signalB.root'
treename = 'tree_vec'
columns = ["EventNumber", "Coretrack_signalv0","Coretrack_signalv1", "Coretrack_pt", "Coretrack_z0", "Coretrack_d0", "Coretrack_nPiHits", "Coretrack_nSiHoles", "Coretrack_ratioptCalo", "Coretrack_dR", "Coretrack_dRleadtrk", "Coretrack_CaloEMpt", "Coretrack_CaloHadpt"]
selection = '( tau_offl_isMedium && Coretrack_ratioptCalo > 0 )'

df_mc = read_root(path+datafile, treename,columns=columns,where=selection,flatten=columns[1:])
df_mc.reset_index(inplace=True)
df_mc.to_feather("pandas_trigger_v0.feather")
print (df_mc.shape )
