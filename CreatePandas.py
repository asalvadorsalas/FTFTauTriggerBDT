from __future__ import print_function
from root_pandas import read_root

path = '/nfs/at3/scratch/salvador/TriggerR21/FTFTauTriggerBDT/datafiles/'
datafile='mc16_13TeV.425200.Pythia8EvtGen_A14NNPDF23LO_Gammatautau_MassWeight.root'
treename = 'tree_vec'
columns = ["EventNumber", "Coretrack_signalv0","Coretrack_signalv1", "Coretrack_pt", "Coretrack_z0", "Coretrack_d0", "Coretrack_nPiHits", "Coretrack_nSiHoles", "Coretrack_ratioptCalo", "Coretrack_dR", "Coretrack_dRleadtrk", "Coretrack_CaloEMpt", "Coretrack_CaloHadpt"]
selection = '( tau_offl_isMediumRNN && Coretrack_ratioptCalo > 0 )'

df_mc = read_root(path+datafile, treename,columns=columns,where=selection,flatten=columns[1:])
df_mc.reset_index(inplace=True)
df_mc.to_feather(path+"pandas_trigger_v0.feather")
print (df_mc.shape )
