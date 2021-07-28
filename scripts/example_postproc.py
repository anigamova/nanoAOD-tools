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

# Function parameters
# (isMC=True, dataYear=2016, runPeriod="B", jesUncert="Total", redojec=False, jetType = "AK4PFchs", noGroom=False)
# All other parameters will be set in the helper module


#fnames = ["/eos/cms/store/mc/RunIISummer16NanoAODv5/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7_ext2-v1/120000/FF69DF6E-2494-F543-95BF-F919B911CD23.root"]
fnames = ["/eos/cms/store/mc/RunIIAutumn18NanoAODv6/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/NANOAODSIM/Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/20000/0694B929-669C-314B-A92B-1F133051723B.root"]
#srm://cmsdcatape.fnal.gov:8443/srm/managerv2?SFN=/11/store/mc/RunIIAutumn18NanoAODv4/ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8/NANOAODSIM/Nano14Dec2018_102X_upgrade2018_realistic_v16-v1/110000/8F9B988E-CEA9-644E-9A39-1D6ECE83E513.root"]
#/eos/cms/store/mc/RunIISummer16NanoAODv5/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7_ext2-v1/120000/FF69DF6E-2494-F543-95BF-F919B911CD23.root"]

# p=PostProcessor(".",fnames,"Jet_pt>150","",[jetmetUncertainties2016(),exampleModuleConstr()],provenance=True)
jmeCorrections2018MC = createJMECorrector(True, "2018", "A", "Merged", "AK4PFchs", False,splitJER=False)
jmeCorrections2018MCAll = createJMECorrector(True, "2018", "A", "All", "AK4PFchs", False,splitJER=True)
jmeCorrections2018AK8MC = createJMECorrector(True, "2018", "A", "Merged", "AK8PFPuppi", False,splitJER=False)
jmeCorrections2018AK8MCAll = createJMECorrector(True, "2018", "A", "All", "AK8PFPuppi", False,splitJER=True)

mhtVHbb = lambda : mhtProducer( lambda j : j.pt > 30,
                            lambda mu : mu.pt > 5 and mu.pfRelIso04_all < 0.4,
                            lambda el : el.pt > 5 and el.pfRelIso03_all < 0.4 )

selection = ""
p=PostProcessor(".",fnames,selection.replace('\n',' '),"keep_and_drop_output.txt",[jmeCorrections2018MC(),jmeCorrections2018MCAll(),jmeCorrections2018AK8MC(),jmeCorrections2018AK8MCAll(),mhtVHbb(),vhbb2018()],provenance=True)

#p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",[puWeight_2018(),jmeCorrections2018MC(),jmeCorrections2018MCAll(),jmeCorrections2018AK8MC(),jmeCorrections2018AK8MCAll(),muonScaleRes2018(),mhtVHbb(),btagSFProducer("2018","deepcsv"),vhbb2018()],provenance=True)

#p = PostProcessor(".", fnames, "Jet_pt>150", "", [
#                  jmeCorrections(), exampleModuleConstr()], provenance=True)
p.run()
print "DONE"
#os.system("ls -lR")

