#!/usr/bin/env python
#from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.reweighting.vhbb_analysis import *  
#from PhysicsTools.NanoAODTools.postprocessing.analysis.higgs.vhbb.KinfitProducer import * 
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
from  PhysicsTools.NanoAODTools.postprocessing.modules.jme.mht import *
from  PhysicsTools.NanoAODTools.postprocessing.modules.common.PrefireCorr import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import *



from importlib import import_module
import os
import sys
import ROOT
import argparse
ROOT.PyConfig.IgnoreCommandLineOptions = True
# soon to be deprecated
# new way of using jme uncertainty

#parser = argparse.ArgumentParser("")
#parser.add_argument('-isMC', '--isMC', type=int, default=1, help="")
#parser.add_argument('-isVjets', '--isVjets', type=int, default=0, help="")
#parser.add_argument('-jobNum', '--jobNum', type=int, default=1, help="")
#parser.add_argument('-era', '--era', type=str, default="2017", help="")
#parser.add_argument('-dataRun', '--dataRun', type=str, default="X", help="")
#args = parser.parse_args()
#print "args = ",args
isMC = True #args.isMC
era = "2018"#args.era
isVjets = False#args.isVjets
dataRun = "X"#args.dataRun

print "isMC = ",isMC,"era = ",era, "dataRun = ",dataRun
desy = '/pnfs/desy.de/cms/tier2/store//user/anigamov//'
cern = '/eos/cms/store/mc/RunIIAutumn18NanoAODv6/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/'

fnames = [
cern+'NANOAODSIM/Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/230000/17227660-5AAA-7844-A6B9-B38362ABE23F.root',
cern+'NANOAODSIM/Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/230000/38D746F3-783A-064E-9AC3-AB2D85032250.root',
cern+'NANOAODSIM/Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/230000/5CE71A7A-B111-DE40-907A-7ABCCDC6315C.root',
cern+'NANOAODSIM/Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/230000/688060C7-151C-4944-82F2-5E691AED6AC9.root',
cern+'NANOAODSIM/Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/230000/99F692D0-AB49-A840-9D79-F8A062D1441C.root',
cern+'NANOAODSIM/Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/230000/DFD2630A-096A-DB4F-932C-C97497CA6C2F.root',
cern+'NANOAODSIM/Nano25Oct2019_102X_upgrade2018_realistic_v20_ext1-v1/230000/3F17120C-D93E-DB47-81E9-37F3080590CE.root',
cern+'NANOAODSIM/Nano25Oct2019_102X_upgrade2018_realistic_v20_ext1-v1/230000/C810B8DD-6B80-C444-96B9-43237A8219B2.root',
cern+'NANOAODSIM/Nano25Oct2019_102X_upgrade2018_realistic_v20_ext1-v1/240000/4D79ED9F-BA4A-7F4D-AEA9-3DC55901686A.root',
cern+'NANOAODSIM/Nano25Oct2019_102X_upgrade2018_realistic_v20_ext1-v1/240000/89DF89EF-660B-1045-8C44-80274AF0B4EA.root',
cern+'NANOAODSIM/Nano25Oct2019_102X_upgrade2018_realistic_v20_ext1-v1/240000/D40F0261-714C-5540-BE83-D9AEC9A02AC3.root',
cern+'NANOAODSIM/Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/20000/0694B929-669C-314B-A92B-1F133051723B.root',
]
jmeCorrections2018MC = createJMECorrector(True, "2018", "A", "Merged", "AK4PFchs", False,splitJER=False)
jmeCorrections2018MCAll = createJMECorrector(True, "2018", "A", "All", "AK4PFchs", False,splitJER=True)
jmeCorrections2018AK8MC = createJMECorrector(True, "2018", "A", "Merged", "AK8PFPuppi", False,splitJER=False)
jmeCorrections2018AK8MCAll = createJMECorrector(True, "2018", "A", "All", "AK8PFPuppi", False,splitJER=True)

mhtVHbb = lambda : mhtProducer( lambda j : j.pt > 30,
                            lambda mu : mu.pt > 5 and mu.pfRelIso04_all < 0.4,
                            lambda el : el.pt > 5 and el.pfRelIso03_all < 0.4 )

selection = ""

p=PostProcessor(".",fnames,cut=selection.replace('\n',' '),branchsel="keep_and_drop_input.txt",modules=[jmeCorrections2018MC(),jmeCorrections2018MCAll(),jmeCorrections2018AK8MC(),jmeCorrections2018AK8MCAll(),mhtVHbb(),vhbb2018()],provenance=True,outputbranchsel="keep_and_drop_output.txt")

#p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",[puWeight_2018(),jmeCorrections2018MC(),jmeCorrections2018MCAll(),jmeCorrections2018AK8MC(),jmeCorrections2018AK8MCAll(),muonScaleRes2018(),mhtVHbb(),btagSFProducer("2018","deepcsv"),vhbb2018()],provenance=True)

#p = PostProcessor(".", fnames, "Jet_pt>150", "", [
#                  jmeCorrections(), exampleModuleConstr()], provenance=True)
p.run()
print "DONE"
#os.system("ls -lR")

