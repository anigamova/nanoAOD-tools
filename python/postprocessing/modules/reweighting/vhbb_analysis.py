import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import sys
import copy
import math

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import * #deltaR, matching etc..
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
from  PhysicsTools.NanoAODTools.postprocessing.modules.jme.mht import *
from PhysicsTools.NanoAODTools.postprocessing.modules.reweighting.applysmearing import SmearApplicator

class VHbbProducer(Module):
    def __init__(self, isMC, era, useCMVA=False,isVjets=False,genOnly=False):
        self.era = era
        self.isMC = isMC
        self.useCMVA = useCMVA
        self.isVjets = isVjets
        self.genOnly = genOnly
        self.smearapplicator = SmearApplicator(year=era, isdata=(not isMC), useV13=False, doscaling=False)
        self.cutDict = {
            "common":{
                "mjj_low":90,
                "mjj_high":150,
                "tag_wpL":0.1241,
                "tag_wpM":0.4184,
                "tag_wpT":0.7527,
                "fatJetPtCut":250,
            },
            0: {
                "ptjj":120,
                "j1ptCut":60,
                "j2ptCut":35,
                "HLTs":["HLT_PFMET120_PFMHT120_IDTight"],
                "METFilters":["Flag_goodVertices","Flag_globalSuperTightHalo2016Filter","Flag_HBHENoiseFilter","Flag_HBHENoiseIsoFilter","Flag_EcalDeadCellTriggerPrimitiveFilter"],
        
            },
            1: {
                "ptjj":100,
                "j1ptCut":25,
                "j2ptCut":25,
                "HLTs":["HLT_Ele32_WPTight_Gsf","HLT_IsoMu24"],
        
            },
            2: {
                "ptjj":0.,
                "j1ptCut":20.,
                "j2ptCut":20.,
                "vmass_low":75,
                "vmass_high":105,
                "hv_dphi_cut":2.5,
                "HLTs":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"],
            }    
        }
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def EquationSolve(self, a, b, c, d):


      result = []



      if (a != 0):

        q = (3*a*c-b*b)/(9*a*a)
        r = (9*a*b*c - 27*a*a*d - 2*b*b*b)/(54*a*a*a)
        Delta = q*q*q + r*r

        rho=0.
        theta=0.

        if( Delta<=0):
          rho = math.sqrt(-(q*q*q))

          theta = acos(r/rho)

          s = complex(math.sqrt(-q)*cos(theta/3.0),math.sqrt(-q)*sin(theta/3.0))
          t = complex(math.sqrt(-q)*cos(-theta/3.0),math.sqrt(-q)*sin(-theta/3.0))

        if(Delta>0):
          #print r, math.sqrt(Delta)
          if (r+math.sqrt(Delta) > 0):
              s = complex(pow((r+math.sqrt(Delta)),(1./3)),0)
          else:
              s = complex(-pow(abs((r+math.sqrt(Delta))),(1./3)),0)
          if (r-math.sqrt(Delta) > 0):
              t = complex(pow((r-math.sqrt(Delta)),(1./3)),0)
          else:
              t = complex(-pow(abs((r-math.sqrt(Delta))),(1./3)),0)

        i = complex(0.,1.0)

        x1 = s+t+complex(-b/(3.0*a),0)
        x2 = (s+t)*complex(-0.5,0)-complex(b/(3.0*a),0)+(s-t)*i*complex(math.sqrt(3)/2.0,0)
        x3 = (s+t)*complex(-0.5,0)-complex(b/(3.0*a),0)-(s-t)*i*complex(math.sqrt(3)/2.0,0)

        if(abs(x1.imag)<0.0001): result.append(x1.real)
        if(abs(x2.imag)<0.0001): result.append(x2.real)
        if(abs(x3.imag)<0.0001): result.append(x3.real)

        #print x1,x2,x3
        return result
      else:
          return result

      return result
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Reco_Cat",   "I");
        self.out.branch("Vtype",   "I");
        self.out.branch("V_pt",    "F");
        self.out.branch("V_eta",    "F");
        self.out.branch("V_phi",   "F");
        self.out.branch("V_mass",  "F");
        self.out.branch("V_mt",    "F");
        self.out.branch("H_pt",  "F");
        self.out.branch("H_eta",  "F");
        self.out.branch("H_phi",  "F");
        self.out.branch("H_mass",  "F");
        self.out.branch("FatJetCandidate_pt",  "F");
        self.out.branch("FatJetCandidate_eta",  "F");
        self.out.branch("FatJetCandidate_phi",  "F");
        self.out.branch("FatJetCandidate_mass",  "F");
        self.out.branch("Jet_lepFilter", "F", 1, "nJet");
        self.out.branch("HJ1_HJ2_dPhi","F")
        self.out.branch("HJ1_HJ2_dEta","F")
        self.out.branch("HJ1_HJ2_dR","F")
        self.out.branch("hJets_leadingPt","F")
        self.out.branch("hJets_subleadingPt","F")
        self.out.branch("hJets_btagWP_0","F")
        self.out.branch("hJets_btagWP_1","F")
        self.out.branch("nAddJets302p4_puid","F")
        self.out.branch("otherJetsHighestPt","F")
        self.out.branch("otherJetsHighestBtag","F")
        self.out.branch("minDPhiFromOtherJets","F")
        self.out.branch("SA5","F")

        self.out.branch("HVdR","F")
        self.out.branch("HVdEta","F")
        self.out.branch("HVdPhi","F")
        self.out.branch("VPtjjRatio","F")
        self.out.branch("jjVPtRatio","F")
        self.out.branch("MET_Pt","F")
        self.out.branch("lepMetDPhi","F")
        self.out.branch("Top1_mass_fromLepton_regPT_w4MET","F")


#        self.out.branch("Jet_CvsL", "F", 1, "nJet")
#        self.out.branch("Jet_CvsB", "F", 1, "nJet")
#        self.out.branch("Jet_DeepFlavCvsL", "F", 1, "nJet")
#        self.out.branch("Jet_DeepFlavCvsB", "F", 1, "nJet")
#        self.out.branch("MET_Pt","F");
#        self.out.branch("MET_Phi","F");
#        self.out.branch("GenJetWithNeutrinos_pt","F",1, "nGenJet")
#        self.out.branch("GenJetWithNeutrinos_eta","F",1, "nGenJet")
#        self.out.branch("GenJetWithNeutrinos_phi","F",1, "nGenJet")
#        self.out.branch("GenJetWithNeutrinos_mass","F",1, "nGenJet")
#
        ## for the boosted analysis
#        self.out.branch("Pt_fjidx",  "I");        
#        self.out.branch("Msd_fjidx",  "I");
#        self.out.branch("Hbb_fjidx",  "I");
#        
#        self.out.branch("SAptfj_HT",  "F");
#        self.out.branch("SAptfj5",  "F");
#        self.out.branch("SAmfj_HT",  "F");
#        self.out.branch("SAmfj5",  "F");
#        self.out.branch("SAhbbfj_HT",  "F");
#        self.out.branch("SAhbbfj5",  "F");
#        
        self.out.branch("FatJet_lepFilter",  "O", 1, "nFatJet");
        self.out.branch("FatJetCandidate_deepTagMD_bbvsLight",  "F");
#        self.out.branch("FatJetCandidate_Pt", "F", 1, "nFatJet");
#        self.out.branch("FatJetCandidate_Msoftdrop", "F", 1, "nFatJet");
#        
#        self.out.branch("FatJetCandidate_FlavourComposition", "I", 1, "nFatJet"); #55 bb, #54 bc, #5 b, #4 c, #44 cc, #1 other
#        
#        self.out.branch("FatJetCandidate_HiggsProducts", "O", 1, "nFatJet");
#        self.out.branch("FatJetCandidate_WProducts", "O", 1, "nFatJet");
#        self.out.branch("FatJetCandidate_ZProducts", "O", 1, "nFatJet");
#        
#        ## Gen information
#        self.out.branch("nGenStatus2bHad", "I");
        self.out.branch("GenBJ1_pt", "F");
        self.out.branch("GenBJ1_genjet_pt", "F");
        self.out.branch("GenBJ1_eta", "F");
        self.out.branch("GenBJ1_phi", "F");
        self.out.branch("GenBJ1_mass", "F");
        self.out.branch("GenBJ1_genjet_mass", "F");
        self.out.branch("GenBJ1_index", "I");
        self.out.branch("GenBJ1_genjet_index", "I");
        self.out.branch("GenBJ2_pt", "F");
        self.out.branch("GenBJ2_genjet_pt", "F");
        self.out.branch("GenBJ2_eta", "F");
        self.out.branch("GenBJ2_phi", "F");
        self.out.branch("GenBJ2_mass", "F");
        self.out.branch("GenBJ2_genjet_mass", "F");
        self.out.branch("GenBJ2_index", "I");
        self.out.branch("GenBJ2_genjet_index", "I");
        self.out.branch("GenBJJ_pt", "F");
        self.out.branch("GenBJJ_genjet_pt", "F");
        self.out.branch("GenBJJ_eta", "F");
        self.out.branch("GenBJJ_phi", "F");
        self.out.branch("GenBJJ_mass", "F");
        self.out.branch("GenBJJ_genjet_mass", "F");
        self.out.branch("GenBJJ_dPhi", "F");
        self.out.branch("GenBJJ_dR", "F");
        self.out.branch("GenBJJ_dEta", "F");
#        self.out.branch("GenLepIndex1","I");
#        self.out.branch("GenLepIndex2","I");
#        self.out.branch("GenLep_GenBJ1_dR","F");
#        self.out.branch("GenLep_GenBJ1_dEta","F");
#        self.out.branch("GenLep_GenBJ1_dPhi","F");
#        self.out.branch("GenTop1_mass","F");
#        self.out.branch("GenLep_GenBJ2_dR","F");
#        self.out.branch("GenLep_GenBJ2_dEta","F");
#        self.out.branch("GenLep_GenBJ2_dPhi","F");
#        self.out.branch("GenTop2_mass","F");
#        self.out.branch("nGenTop","I")
#        self.out.branch("nW","I")
#        self.out.branch("nWlep","I")
#        self.out.branch("nGenVBosons","I")
#        self.out.branch("LeadGenVBoson_pt","F")
#        self.out.branch("LeadGenVBoson_eta","F")
#        self.out.branch("LeadGenVBoson_phi","F")
#        self.out.branch("LeadGenVBoson_pdgId","F")
#        self.out.branch("LeadNeutrinoFromGenVBoson_pt", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def matchSoftActivity(self,jets,saJets,dR=0.4) :
      matched=set()
      for saj in saJets:
          for j in jets :
              if deltaR(saj,j) < dR :
                        matched.add(saj)
      return matched
    
    def matchSoftActivityFSR(self,jet1,jet2,saJets,dR=0.4) :
        matched=set()
        drjj = deltaR(jet1,jet2)
        sumDeltaRMin = drjj + 2*dR
        for saj in saJets:
            dr1 = deltaR(saj,jet1)
            dr2 = deltaR(saj,jet2)
            if ((dr1+dr2) < sumDeltaRMin):
                matched.add(saj)
        return matched
            
    def pt(self, jet, isMC, noReg=False, sysVar=0):
        ## the MC has JER smearing applied which has output branch Jet_pt_nom which should be compared 
        ## with data branch Jet_pt. This essentially aliases the two branches to one common jet pt variable.
        if noReg:
            if isMC:
                return jet.pt_nom
            else: #jet_pt_nom contains re-corrected jet pT if jecRecalibrator has been run
                if hasattr(jet, 'pt_nom'):
                    return jet.pt_nom
                else:
                    return jet.pt    
        else:
#            genJet=None
#            if isMC and jet.genJetIdx >=0 and  jet.genJetIdx < len(self.genJetsWithNeutrinos) :
#                genJet=self.genJetsWithNeutrinos[jet.genJetIdx]

            jet_pt = jet.pt_nom if hasattr(jet, 'pt_nom') else jet.pt
            return jet_pt
            # Not using rho at the moment (luckily?)
            # If initialized for data, does not do any smearing
#            pt_variations = self.smearapplicator.get_smear(jet_pt, jet.bRegCorr, genJet.Pt() if genJet else 0)
#
#            if sysVar==0: # nominal
#                return pt_variations.nominal
#
#            elif sysVar==1: # up
#                return pt_variations.up
#
#            elif sysVar==-1: # down
#                return pt_variations.down
    
    def met(self, met, isMC):
        ## the MC has JER smearing applied which has output branch met_[pt/phi]_nom which should be compared 
        ## with data branch MET_[pt/phi]. This essentially aliases the two branches to one common variable.
        if isMC:
            if hasattr(met, 'pt_nom'):
                return (met.pt_nom,met.phi_nom)
            else:
                return (met.pt,met.phi)
        else:
            if hasattr(met, 'pt_nom'):#pt_nom contains re-corrected pT if jecRecalibrator has been run
                return (met.pt_nom,met.phi_nom)
            else:
                return (met.pt,met.phi)
   
    def msoftdrop(self, jet, isMC):
        ## the MC has JER smearing applied which has output branch Jet_pt_smeared which should be compared 
        ## with data branch Jet_pt. This essentially aliases the two branches to one common jet pt variable.
        if isMC:
            return jet.msoftdrop_nom*0.9
        else:
            return jet.msoftdrop
 
    def btag(self, jet):
        if (self.useCMVA):
            return jet.btagCMVA
        else:
            return jet.btagDeepB

    def cvsltag(self, jet):
        btagDeepL = 1.-(jet.btagDeepC+jet.btagDeepB)
        if jet.btagDeepB >= 0. and jet.btagDeepB < 1. and jet.btagDeepC >= 0. and btagDeepL >= 0.:
            return jet.btagDeepC/(1.-jet.btagDeepB)
        else:
            return -1

    def cvsbtag(self, jet):
        btagDeepL = 1.-(jet.btagDeepC+jet.btagDeepB)
        if jet.btagDeepB > 0. and jet.btagDeepC > 0. and btagDeepL >= 0.:
            return jet.btagDeepC/(jet.btagDeepC+jet.btagDeepB)
        else:
            return -1

    def deepflavcvsltag(self, jet):
        if not hasattr(jet, 'btagDeepFlavC'):           # only Nano V5 onwards
            return -99. 
        btagDeepFlavL = 1.-(jet.btagDeepFlavC+jet.btagDeepFlavB)
        if jet.btagDeepFlavB >= 0. and jet.btagDeepFlavB < 1. and jet.btagDeepFlavC >= 0. and btagDeepFlavL >= 0.:
            return jet.btagDeepFlavC/(1.-jet.btagDeepFlavB)
        else:
            return -1

    def deepflavcvsbtag(self, jet):
        if not hasattr(jet, 'btagDeepFlavC'):           # only Nano V5 onwards
            return -99.
        btagDeepFlavL = 1.-(jet.btagDeepFlavC+jet.btagDeepFlavB)
        if jet.btagDeepFlavB > 0. and jet.btagDeepFlavC > 0. and btagDeepFlavL >= 0.:
            return jet.btagDeepFlavC/(jet.btagDeepFlavC+jet.btagDeepFlavB)
        else:
            return -1

    def elid(self, el, wp):
        if (wp == "80"):
            return el.mvaFall17V2Iso_WP80
        elif (wp == "90"):
            return el.mvaFall17V2Iso_WP90

    def computeTopMass_new(self, lep, met, jet, METmethod):
        new_top = ROOT.TLorentzVector()
        jet_transverse = ROOT.TLorentzVector()
        lep_transverse = ROOT.TLorentzVector()
        if METmethod == 1:
           jet_transverse.SetPxPyPzE(jet.Px(),jet.Py(),0,math.sqrt((jet.M())**2 + (jet.Pt())**2))   
           lep_transverse.SetPxPyPzE(lep.Px(),lep.Py(),0,math.sqrt((lep.M())**2 + (lep.Pt())**2))   

           new_top = jet_transverse + lep_transverse + met
        elif METmethod == 2:
           neutrino = ROOT.TLorentzVector()
           neutrino = self.getNu4Momentum(lep, met)
           new_top = lep + jet + neutrino
        elif METmethod == 3:
           new_top = lep + jet + met

        return new_top.M()
    def get_topMass(self, jets, lepton,met):
      minDR = 999
      closestJet_4v = ROOT.TLorentzVector()
      closestJet = -1;
      for jet in jets:
          if (self.pt(jet,self.isMC,noReg=False)>30 and self.btag(jet)>self.cutDict["common"]['tag_wpL'] and jet.lepFilter and ((jet.puId>0 and self.pt(jet,self.isMC,noReg=True)>50) or jet.puId>6) and jet.jetId>0): 
              dR = jet.DeltaR(lepton)
              if (dR < minDR):
                  minDR = dR
                  closestJet = jet
      if closestJet!=-1: 
          closestJet_4v.SetPtEtaPhiM(self.pt(closestJet,self.isMC,noReg=False), closestJet.eta, closestJet.phi, closestJet.mass)
          top_mass = self.computeTopMass_new(lepton, met, closestJet_4v, 2)
      else: top_mass=-99.
      return top_mass
    def getNu4Momentum(self, TLepton, TMET, debug=False):

      branch = -1

      Lepton = ROOT.TLorentzVector()
      Lepton.SetPxPyPzE(TLepton.Px(), TLepton.Py(), TLepton.Pz(), TLepton.E());
      MET = ROOT.TLorentzVector()
      MET.SetPxPyPzE(TMET.Px(), TMET.Py(), 0., TMET.E());

      mW = 80.38;

      result = []

      MisET2 = (MET.Px()*MET.Px() + MET.Py()*MET.Py());
      mu = (mW*mW)/2 + MET.Px()*Lepton.Px() + MET.Py()*Lepton.Py();
      a  = (mu*Lepton.Pz())/(Lepton.Energy()*Lepton.Energy() - Lepton.Pz()*Lepton.Pz());
      a2 = pow(a,2);
      b  = (pow(Lepton.Energy(),2.)*(MisET2) - pow(mu,2.))/(pow(Lepton.Energy(),2) - pow(Lepton.Pz(),2));
      pz1 = 0.
      pz2 = 0.
      pznu = 0.
      nNuSol = 0

      p4nu_rec = ROOT.TLorentzVector()
      p4W_rec = ROOT.TLorentzVector()
      p4b_rec = ROOT.TLorentzVector()
      p4Top_rec = ROOT.TLorentzVector()
      p4lep_rec = ROOT.TLorentzVector()

      p4lep_rec.SetPxPyPzE(Lepton.Px(),Lepton.Py(),Lepton.Pz(),Lepton.Energy());

      #print a2,b
      if(a2-b > 0 ):
        root = math.sqrt(a2-b);
        pz1 = a + root;
        pz2 = a - root;
        nNuSol = 2;

        pznu = pz1;
        branch = 1
        if abs(pz2)<abs(pz1):
            pznu = pz2
            branch = 2

        Enu = math.sqrt(MisET2 + pznu*pznu);

        p4nu_rec.SetPxPyPzE(MET.Px(), MET.Py(), pznu, Enu);

        result.append(p4nu_rec);

      else:

        ptlep = Lepton.Pt()
        pxlep=Lepton.Px()
        pylep=Lepton.Py()
        metpx=MET.Px()
        metpy=MET.Py()

        EquationA = 1.
        EquationB = -3.*pylep*mW/(ptlep)
        EquationC = mW*mW*(2*pylep*pylep)/(ptlep*ptlep)+mW*mW-4*pxlep*pxlep*pxlep*metpx/(ptlep*ptlep)-4*pxlep*pxlep*pylep*metpy/(ptlep*ptlep)
        EquationD = 4.*pxlep*pxlep*mW*metpy/(ptlep)-pylep*mW*mW*mW/ptlep

        solutions = self.EquationSolve(EquationA,EquationB,EquationC,EquationD)

        solutions2 = self.EquationSolve(EquationA,-EquationB,EquationC,-EquationD);

        deltaMin = 14000*14000
        zeroValue = -mW*mW/(4*pxlep)
        minPx=0
        minPy=0

        for i in range(len(solutions)):
            if(solutions[i]<0 ): continue
            p_x = (solutions[i]*solutions[i]-mW*mW)/(4*pxlep)
            p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x -mW*ptlep*solutions[i])/(2*pxlep*pxlep)
            Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy)

            if(Delta2< deltaMin and Delta2 > 0):
                deltaMin = Delta2
                minPx=p_x
                minPy=p_y
                branch = 3

        for i in range(len(solutions2)):
            if(solutions2[i]<0 ): continue
            p_x = (solutions2[i]*solutions2[i]-mW*mW)/(4*pxlep)
            p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x +mW*ptlep*solutions2[i])/(2*pxlep*pxlep)
            Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy)
            if(Delta2< deltaMin and Delta2 > 0):
                deltaMin = Delta2
                minPx=p_x
                minPy=p_y
                branch = 4

        pyZeroValue= ( mW*mW*pxlep + 2*pxlep*pylep*zeroValue)
        delta2ZeroValue= (zeroValue-metpx)*(zeroValue-metpx) + (pyZeroValue-metpy)*(pyZeroValue-metpy)

        if(deltaMin==14000*14000):
            if debug:
                branch = 6
                return TLorentzVector(0,0,0,0), branch
            else:
                return TLorentzVector(0,0,0,0)

        if(delta2ZeroValue < deltaMin):
            deltaMin = delta2ZeroValue
            minPx=zeroValue
            minPy=pyZeroValue
            branch = 5

        mu_Minimum = (mW*mW)/2 + minPx*pxlep + minPy*pylep
        a_Minimum  = (mu_Minimum*Lepton.Pz())/(Lepton.Energy()*Lepton.Energy() - Lepton.Pz()*Lepton.Pz())
        pznu = a_Minimum

        Enu = math.sqrt(minPx*minPx+minPy*minPy + pznu*pznu)
        p4nu_rec.SetPxPyPzE(minPx, minPy, pznu , Enu)
        result.append(p4nu_rec)

      if debug:
          return result[0], branch
      else:
          return result[0]

    def statusFlags_dict(self, bit):
        Dict={0 : "isPrompt", 1 : "isDecayedLeptonHadron", 2 : "isTauDecayProduct", 3 : "isPromptTauDecayProduct", 4 : "isDirectTauDecayProduct", 5 : "isDirectPromptTauDecayProduct", 6 : "isDirectHadronDecayProduct", 7 : "isHardProcess", 8 : "fromHardProcess", 9 : "isHardProcessTauDecayProduct", 10 : "isDirectHardProcessTauDecayProduct", 11 : "fromHardProcessBeforeFSR", 12 : "isFirstCopy", 13 : "isLastCopy", 14 : "isLastCopyBeforeFSR" }
        return Dict[bit] 
    def select_leptons(self, event):
        electrons = list(Collection(event, "Electron"))
        muons = list(Collection(event, "Muon"))
        wElectrons = [x for x in electrons if self.elid(x,"80") and x.pt > 25 and x.pfRelIso03_all < 0.12]
        wMuons = [x for x in muons if x.pt > 25 and x.tightId >= 1 and x.pfRelIso04_all < 0.15 and abs(x.dxy) < 0.05 and abs(x.dz) < 0.2]
        zElectrons = [x for x in electrons if x.pt > 20 and self.elid(x,"90") and x.pfRelIso03_all < 0.15]
        zMuons = [x for x in muons if x.pt > 20 and x.pfRelIso04_all < 0.25 and abs(x.dxy) < 0.05 and abs(x.dz) < 0.2] # muons already preselected with looseId requirement
        zMuons.sort(key=lambda x:x.pt,reverse=True)
        zElectrons.sort(key=lambda x:x.pt,reverse=True)
        return zMuons,wMuons,zElectrons,wElectrons 
    def select_vboson(self, event, jets):
        electrons = list(Collection(event, "Electron"))
        muons = list(Collection(event, "Muon"))  
        wElectrons = [x for x in electrons if self.elid(x,"80") and x.pt > 25 and x.pfRelIso03_all < 0.12]
        wMuons = [x for x in muons if x.pt > 25 and x.tightId >= 1 and x.pfRelIso04_all < 0.15 and abs(x.dxy) < 0.05 and abs(x.dz) < 0.2]
        zElectrons = [x for x in electrons if x.pt > 20 and self.elid(x,"90") and x.pfRelIso03_all < 0.15]
        zMuons = [x for x in muons if x.pt > 20 and x.pfRelIso04_all < 0.25 and abs(x.dxy) < 0.05 and abs(x.dz) < 0.2] # muons already preselected with looseId requirement
        zMuons.sort(key=lambda x:x.pt,reverse=True)
        zElectrons.sort(key=lambda x:x.pt,reverse=True)
        wMuons.sort(key=lambda x:x.pt,reverse=True)
        wElectrons.sort(key=lambda x:x.pt,reverse=True)
        mht = sum([self.pt(x,self.isMC) for x in jets if x.lepFilter and x.puId>0 and x.jetId>0 and self.pt(x,self.isMC)>30 and abs(x.eta)<2.5])
        met = Object(event, "MET")
        metPt,metPhi = self.met(met,self.isMC)
        #self.out.fillBranch("MET_Pt",metPt)
        #self.out.fillBranch("MET_Phi",metPhi)
        Vtype = -1
        V_mt = -99.
        Top1_mass_fromLepton_regPT_w4MET = -99.
        vLeptons = [] # decay products of V
        vLidx = [-1,-1] # indices in lepton collection of selected leptons
        if len(zMuons) >= 2:
            if zMuons[0].pt > 20:
                for i in xrange(1,len(zMuons)):
                    if zMuons[0].charge * zMuons[i].charge < 0:
                        Vtype = 0
                        vLeptons = [zMuons[0],zMuons[i]]
                        vLidx[0] = muons.index(zMuons[0])
                        vLidx[1] = muons.index(zMuons[i])
                        break
        elif len(zElectrons) >= 2:
            if zElectrons[0].pt > 20:
                for i in xrange(1,len(zElectrons)):
                    if zElectrons[0].charge * zElectrons[i].charge < 0:
                        Vtype = 1
                        vLeptons = [zElectrons[0],zElectrons[i]]
                        vLidx[0] = electrons.index(zElectrons[0])
                        vLidx[1] = electrons.index(zElectrons[i])
                        break
        elif len(wElectrons) + len(wMuons) == 1:
            lep_met_dPhi = -1;
            if len(wMuons) == 1:
                Vtype = 2
                vLeptons = [wMuons[0]]
                vLidx[0] = muons.index(wMuons[0])
                lep_met_dPhi = deltaPhi(wMuons[0].phi,metPhi)
            if len(wElectrons) == 1:
                Vtype=3
                vLeptons = [wElectrons[0]]
                vLidx[0] = electrons.index(wElectrons[0])
                lep_met_dPhi = deltaPhi(wElectrons[0].phi,metPhi)
            n_addLeptons = (len([x for k,x in enumerate(electrons) if x.pt > 15 and x.eta<2.5 and self.elid(x,"90") and x.pfRelIso03_all < 0.1 and k!=vLidx[0]]) + len([x for i, x in enumerate(muons) if x.pt >15 and x.eta<2.5 and x.pfRelIso04_all < 0.1  and i!=vLidx[0]])) 
            if not (n_addLeptons==0 and lep_met_dPhi>2.):Vtype-1
            self.out.fillBranch("lepMetDPhi",lep_met_dPhi)            
        elif len(zElectrons) + len(zMuons) > 0:
            Vtype = 5
        else:
            if ((metPt> 170) and min(metPt,mht)>100 and  self.pass_filters(event)):
                nJetsCloseToMET = len([x for x in jets if x.lepFilter and x.puId>0 and x.jetId>0 and self.pt(x,self.isMC)>30 and abs(x.eta)<2.5 and deltaPhi(x.phi,metPhi)<0.5])
                n_addLeptons = (len([x for x in electrons if x.pt > 15 and x.eta<2.5 and self.elid(x,"90") and x.pfRelIso03_all < 0.1]) + len([x for x in muons if x.pt >15 and x.eta<2.5 and x.pfRelIso04_all < 0.1])) 
                if not (nJetsCloseToMET==0 and n_addLeptons==0): Vtype = -1
                else: Vtype = 4
        if (Vtype>=0 and Vtype<5):
             if not self.pass_trigger(event, Vtype):return [0,0,-1,vLeptons]
        else: return [0,0,-1,vLeptons]
        V = ROOT.TLorentzVector()
        met_4vec = ROOT.TLorentzVector()
        met_4vec.SetPtEtaPhiM(metPt,0.,metPhi,0.) # only use met vector to derive transverse quantities
        for vLepton in vLeptons:
            vLepton_4vec = ROOT.TLorentzVector()
            vLepton_4vec.SetPtEtaPhiM(vLepton.pt,vLepton.eta,vLepton.phi,vLepton.mass)
            if Vtype==3 or Vtype==2: Top1_mass_fromLepton_regPT_w4MET = self.get_topMass(jets,vLepton_4vec,met_4vec)
            else: Top1_mass_fromLepton_regPT_w4MET = -99.
            cosPhi12 = (vLepton_4vec.Px()*met_4vec.Px() + vLepton_4vec.Py()*met_4vec.Py()) / (vLepton_4vec.Pt() * met_4vec.Pt())
            V_mt= math.sqrt(2*vLepton_4vec.Pt()*met_4vec.Pt() *max(0.0, (1 - cosPhi12)));
            V = V + vLepton_4vec
        if Vtype >=2 and Vtype<=4:
            V = V + met_4vec
        self.out.fillBranch("V_mt",V_mt)
        self.out.fillBranch("V_pt",V.Pt())
        self.out.fillBranch("V_eta",V.Eta())
        self.out.fillBranch("V_phi",V.Phi())
        self.out.fillBranch("V_mass",V.M())
        self.out.fillBranch("Vtype",Vtype)
        self.out.fillBranch("Top1_mass_fromLepton_regPT_w4MET",Top1_mass_fromLepton_regPT_w4MET)
        return V,met_4vec,Vtype,vLeptons

    def select_higgs_jets(self, event,jets,channel,vLeptons):
        sa = Collection(event, "SoftActivityJet")
        resolved= -1
        failed = [-1,-1,-1,-1]
        jetsForHiggs = [x for x in jets if x.lepFilter and x.puId>0 and x.jetId>0 and self.pt(x,self.isMC)>20 and abs(x.eta)<2.5]
        hjets_btag = []
        if (len(jetsForHiggs) >= 2):
            hJets = sorted(jetsForHiggs, key = lambda jet : self.btag(jet), reverse=True)[0:2]
            hJidx = [jets.index(x) for x in hJets]
            hjets_btag = [self.btag(x) for x in hJets]
            hjets_pt = [self.pt(x,self.isMC,noReg=False) for x in hJets]
            if not (hjets_btag[0]>self.cutDict["common"]['tag_wpM'] and hjets_btag[1]>self.cutDict["common"]['tag_wpL'] and hjets_pt[0]>self.cutDict[channel]['j1ptCut'] and hjets_pt[1]>self.cutDict[channel]['j2ptCut']): 
                resolved = -1
                return failed
            if not (hjets_pt[0]>self.cutDict[channel]['j1ptCut'] and hjets_pt[1]>self.cutDict[channel]['j2ptCut']): return failed
            ## Save a few basic reco. H kinematics
            hj1 = ROOT.TLorentzVector()
            hj2 = ROOT.TLorentzVector()
            hj1.SetPtEtaPhiM(self.pt(jets[hJidx[0]],self.isMC),jets[hJidx[0]].eta,jets[hJidx[0]].phi,jets[hJidx[0]].mass)
            hj2.SetPtEtaPhiM(self.pt(jets[hJidx[1]],self.isMC),jets[hJidx[1]].eta,jets[hJidx[1]].phi,jets[hJidx[1]].mass)
            hbb = hj1 + hj2
            ## try to recover FSR
            jetsFromFSR = []
            for ijet in xrange(len(jets)):
                if ijet == hJidx[0] or ijet == hJidx[1]: continue
                jet = jets[ijet]
                if self.pt(jet,self.isMC,noReg=True)>20 and abs(jet.eta)<3.0 and ((jet.puId>0 and self.pt(jet,self.isMC,noReg=True)>50) or jet.puId>6) and jet.jetId>0 and jet.lepFilter:
                   if min(deltaR(jet,jets[hJidx[0]]),deltaR(jet,jets[hJidx[1]])) < 0.8:
                       jetsFromFSR.append(jet)
            HFSR = hbb
            for jet in jetsFromFSR:
                fsrJetToAdd = ROOT.TLorentzVector()
                fsrJetToAdd.SetPtEtaPhiM(self.pt(jet,self.isMC,noReg=True),jet.eta,jet.phi,jet.mass)
                HFSR = HFSR + fsrJetToAdd
            if not (HFSR.M()>self.cutDict["common"]['mjj_low'] and HFSR.M()<self.cutDict["common"]['mjj_high']):
                resolved = -1
                return failed
            HJ1_HJ2_dPhi = abs(hj1.DeltaPhi(hj2));
            HJ1_HJ2_dEta = abs(hj1.Eta() - hj2.Eta())
            HJ1_HJ2_dR   = hj1.DeltaR(hj2)

            resolved = 1

            toVeto = hJets + vLeptons
            matchedSAJets=self.matchSoftActivity(toVeto,sa)
            matchedSAJetsPt5=[x for x in matchedSAJets if x.pt>5]
            softActivityJetNjets5=event.SoftActivityJetNjets5-len(matchedSAJetsPt5)

            addJetspt30eta2p5 = [x for i,x in enumerate(jets) if x.lepFilter and x.puId>0 and x.jetId>0 and self.pt(x,self.isMC)>30 and abs(x.eta)<2.5 and i!=hJidx[0] and i!=hJidx[1]]
            addJetspt30eta2p4 = [x for i,x in enumerate(jets) if x.lepFilter and x.puId>0 and x.jetId>0 and self.pt(x,self.isMC)>30 and abs(x.eta)<2.4 and i!=hJidx[0] and i!=hJidx[1]]
            otherJets = [x for i,x in enumerate(jets) if x.lepFilter and ((x.puId>0 and self.pt(x,self.isMC,noReg=True)>50) or x.puId>6) and x.jetId>0 and self.pt(x,self.isMC)>30 and i!=hJidx[0] and i!=hJidx[1]]
            if (len(otherJets)>0):
                otherJetsHighestBtag = self.btag(sorted(otherJets, key = lambda x : self.btag(x), reverse=True)[0])
                otherJetsHighestPt = self.pt(sorted(otherJets, key = lambda x : self.pt(x,self.isMC,noReg=True), reverse=True)[0],self.isMC,noReg=True)
                met = Object(event, "MET")
                metPhi = self.met(met,self.isMC)[1]
                minDphi = abs(deltaPhi((sorted(otherJets, key = lambda x : abs(deltaPhi(x,metPhi)), reverse=False)[0]),metPhi))
                self.out.fillBranch("otherJetsHighestPt",otherJetsHighestPt)
                self.out.fillBranch("otherJetsHighestBtag",otherJetsHighestBtag)
                self.out.fillBranch("minDPhiFromOtherJets",minDphi)
#                print(otherJetsHighestBtag, otherJetsHighestPt, minDphi)
#                if (len(otherJets)>1):
#                    print otherJetsHighestBtag, self.btag(sorted(otherJets, key = lambda x : self.btag(x), reverse=True)[1])#,self.btag(sorted(otherJets, key = lambda x : self.btag(x), reverse=True)[2])
#                    print otherJetsHighestPt,self.pt(sorted(otherJets, key = lambda x : self.pt(x,self.isMC,noReg=True), reverse=True)[1],self.isMC,noReg=True)#,self.pt(sorted(otherJets, key = lambda x : self.pt(x,self.isMC,noReg=True)[2]),self.isMC,noReg=True)
##                    print minDphi, sorted(otherJets, key = lambda x : x.deltaPhi(metPhi), reverse=False)[1].deltaPhi(metPhi)#,sorted(otherJets, key = lambda x : x.phi.deltaPhi(metPhi), reverse=False)[2].phi.delta
            else:
                self.out.fillBranch("otherJetsHighestPt",-99.)
                self.out.fillBranch("otherJetsHighestBtag",-99.)
                self.out.fillBranch("minDPhiFromOtherJets",-99.)
    
            self.out.fillBranch("H_pt",HFSR.Pt())
            self.out.fillBranch("H_phi",HFSR.Phi())
            self.out.fillBranch("H_eta",HFSR.Eta())
            self.out.fillBranch("H_mass",HFSR.M())
            self.out.fillBranch("HJ1_HJ2_dPhi",HJ1_HJ2_dPhi)
            self.out.fillBranch("HJ1_HJ2_dEta",HJ1_HJ2_dEta)
            self.out.fillBranch("HJ1_HJ2_dR",HJ1_HJ2_dR)
            self.out.fillBranch("hJets_leadingPt",hj1.Pt())
            self.out.fillBranch("hJets_btagWP_0",hjets_btag[0])
            self.out.fillBranch("hJets_btagWP_1",hjets_btag[1])
            self.out.fillBranch("hJets_subleadingPt",hj2.Pt())
            self.out.fillBranch("nAddJets302p4_puid",len(addJetspt30eta2p4))
            self.out.fillBranch("SA5",softActivityJetNjets5)
            return resolved,HFSR,len(addJetspt30eta2p5),len(addJetspt30eta2p4)

        else:
            return failed
            self.out.fillBranch("H_pt",-99.)
            self.out.fillBranch("H_phi",-99.)
            self.out.fillBranch("H_eta",-99.)
            self.out.fillBranch("H_mass",-99.)
            self.out.fillBranch("HJ1_HJ2_dPhi",-99.)
            self.out.fillBranch("HJ1_HJ2_dEta",-99.)
            self.out.fillBranch("HJ1_HJ2_dR",-99.)
            self.out.fillBranch("hJets_leadingPt",-99.)
            self.out.fillBranch("hJets_subleadingPt",-99.)
            self.out.fillBranch("hJets_btagWP_0",-99.)
            self.out.fillBranch("hJets_btagWP_1",-99.)
            self.out.fillBranch("nAddJets302p4_puid",-99.)
            self.out.fillBranch("SA5",-99.)

    def select_fat_jets(self,fatjets,jets): 
        boosted = -1;BH=0;
        fatjetPts = [-99.]*len(fatjets)
        fatjetMSD = [-99.]*len(fatjets)
        fatjet_bbvsLight = [-99.]*len(fatjets)
        for i in xrange(len(fatjets)):
            fatjetPts[i] = self.pt(fatjets[i],self.isMC,True)
            fatjetMSD[i] = self.msoftdrop(fatjets[i],self.isMC)
            fatjet_bbvsLight[i] = fatjets[i].deepTagMD_bbvsLight

        fatjetsForHiggs = [x for x in fatjets if x.lepFilter and x.jetId>0 and self.pt(x,self.isMC,True)>250 and self.msoftdrop(x,self.isMC)>50 and abs(x.eta)<2.5 and  x.deepTagMD_bbvsLight>0.8]
  
        if (len(fatjetsForHiggs) >= 1):

            jh = sorted(fatjetsForHiggs, key = lambda jet : jet.deepTagMD_bbvsLight, reverse=True)
            bbvsLight_idx = fatjets.index(jh[0])

            BH = ROOT.TLorentzVector()
            BH.SetPtEtaPhiM(self.pt(jh[0],self.isMC,True),jh[0].eta,jh[0].phi,self.msoftdrop(jh[0],self.isMC))
            self.out.fillBranch("FatJetCandidate_pt",BH.Pt())
            self.out.fillBranch("FatJetCandidate_phi",BH.Phi())
            self.out.fillBranch("FatJetCandidate_eta",BH.Eta())
            self.out.fillBranch("FatJetCandidate_mass",BH.M())
            self.out.fillBranch("FatJetCandidate_deepTagMD_bbvsLight",jh[0].deepTagMD_bbvsLight)
            boosted = 1;
            bJetsOutsideFatJet = []
            for ijet in xrange(len(jets)):
                jet = jets[ijet]
                if self.pt(jet,self.isMC,noReg=False)>25 and abs(jet.eta)<2.5 and jet.lepFilter and self.btag(jet)>self.cutDict["common"]['tag_wpM']:
                   if deltaR(jet,jh[0]) > 0.8:
                       bJetsOutsideFatJet.append(jet)
            nBJetsOutsideFatJet = len(bJetsOutsideFatJet)
            return boosted,BH,nBJetsOutsideFatJet
        else:
            return -1,-1,-1
            self.out.fillBranch("FatJetCandidate_pt",-99.)
            self.out.fillBranch("FatJetCandidate_phi",-99.)
            self.out.fillBranch("FatJetCandidate_eta",-99.)
            self.out.fillBranch("FatJetCandidate_mass",-99.)
            self.out.fillBranch("FatJetCandidate_deepTagMD_bbvsLight",-99.)
    def categorise_event(self, lepton_channel,higgs,V,MET,n_add_jets,resolved,boosted):
        channel=-1;channel_STXS=-1
        V_pt = V.Pt()
        V_mass = V.M()
        MET_Pt = MET.Pt()
        H_pt = higgs.Pt()
        H_mass = higgs.M()
        jjVPtRatio = H_pt/V_pt
        VPtjjRatio = V_pt/H_pt
        HVdPhi = abs(higgs.DeltaPhi(V));
        HVdEta = abs(higgs.Eta() - V.Eta())
        HVdR   = higgs.DeltaR(V)
        if resolved >0:
          if not (V_pt>75):
              channel=-1
              return channel,channel_STXS
          channel=lepton_channel
          if channel==2 and V_pt>75 and HVdPhi>self.cutDict[channel]["hv_dphi_cut"] and V_mass>self.cutDict[channel]["vmass_low"] and V_mass<self.cutDict[channel]["vmass_high"] : 
              if (V_pt>75 and V_pt<150 ): channel_STXS=channel*10 + 1
              elif ( V_pt>150 and V_pt<250 and n_add_jets<1): channel_STXS = channel*10 + 2
              elif ( V_pt>150 and V_pt<250 and n_add_jets>0): channel_STXS = channel*10 + 3
              elif (V_pt>250):channel_STXS = channel*10 + 4
              else: 
                  channel=-1
                  channel_STXS = -1
          elif channel==1: 
              if (n_add_jets<2):
                  if (V_pt<250): channel_STXS=channel*10 + 1
                  elif (V_pt>250):channel_STXS=channel*10 + 2
                  else: channel=-1
              else: channel=-1
          elif (channel==0): 
              if (n_add_jets<2):
                  if (V_pt<250 and n_add_jets<1): channel_STXS=channel*10 + 1
                  elif (V_pt<250 and n_add_jets>0): channel_STXS=channel*10+ 2
                  elif (V_pt>250):channel_STXS=channel*10+3
                  else: 
                      channel=-1
                      channel_STXS=-1
              else: channel=-1
          else: return -1,-1
        elif boosted>0:
            if not (HVdPhi>1.57 and V_pt>250): return channel,channel_STXS
            channel=lepton_channel
            if channel==2:
              if not (V_mass>self.cutDict[channel]["vmass_low"] and V_mass<self.cutDict[channel]["vmass_high"]):
                  return -1,-1
              if V_pt<400: channel_STXS = 100+channel*10+5
              else: channel_STXS = 100 + channel*10 + 6
            elif channel==0 or channel==1 and (n_add_jets==0):
              if V_pt<400: channel_STXS = 100+channel*10+5
              else: channel_STXS = 100 + channel*10 + 6
        if channel>-1 and channel_STXS>-1:
            self.out.fillBranch("HVdR",HVdR)
            self.out.fillBranch("HVdEta",HVdEta)
            self.out.fillBranch("HVdPhi",HVdPhi)
            self.out.fillBranch("VPtjjRatio",VPtjjRatio)
            self.out.fillBranch("jjVPtRatio",jjVPtRatio)
            self.out.fillBranch("MET_Pt",MET_Pt)
        else:
            self.out.fillBranch("HVdR",-99.)
            self.out.fillBranch("HVdEta",-99.)
            self.out.fillBranch("HVdPhi",-99.)
            self.out.fillBranch("VPtjjRatio",-99.)
            self.out.fillBranch("jjVPtRatio",-99.)
            self.out.fillBranch("MET_Pt",-99.)
        return channel,channel_STXS

    def lepton_channel(self,Vtype):
        lepton_ch = -1
        if Vtype>-1 and Vtype<2: lepton_ch = 2
        elif Vtype>1 and Vtype<4: lepton_ch = 1
        elif Vtype==4: lepton_ch = 0
        else: lepton_ch= -1
        return lepton_ch
    def pass_trigger(self, event, Vtype):
        channel = self.lepton_channel(Vtype)
        hlt_list = [getattr(event, path) for path in self.cutDict[channel]["HLTs"]]  
        if any(hlt_list):return True
        else: return False
    def pass_filters(self, event):
        filters_list = [getattr(event, path) for path in self.cutDict[0]["METFilters"]]
        if all(filters_list):return True
        else: return False
    def get_gen_info(self, event):
        genParticles = Collection(event, "GenPart")
        if self.isMC:
            ##Info on b-quarks from top or Higgs decay
            genb_id_1=-999;
            genb_id_2=-999;
            genBParts = [];
            genBPartsShallow = [];
            GenBJ1 = ROOT.TLorentzVector()
            GenBJ2 = ROOT.TLorentzVector()
            GenBJ1_genjet = ROOT.TLorentzVector()
            GenBJ2_genjet = ROOT.TLorentzVector()
            for genP in genParticles:
                if(abs(genP.pdgId)==5 and genP.genPartIdxMother > -1 and (abs(genParticles[genP.genPartIdxMother].pdgId)==25 or abs(genParticles[genP.genPartIdxMother].pdgId)==6)):
                    genBParts.append(copy.copy(genP))
                    genBPartsShallow.append(genP)#Just to extract the IDs in the loop again
            genBParts.sort(key=lambda x:x.pt,reverse=True)
            if len(genBPartsShallow) > 1:
                for i in range(0,len(genParticles)):
                    if(genBPartsShallow[0]==genParticles[i]): genb_id_1=i
                    if(genBPartsShallow[1]==genParticles[i]): genb_id_2=i
            self.out.fillBranch("GenBJ1_index",genb_id_1)
            self.out.fillBranch("GenBJ2_index",genb_id_2)
            for part in genBParts:
                part.mass=4.2
            if len (genBParts) > 1:
                GenBJ1.SetPtEtaPhiM(genBParts[0].pt,genBParts[0].eta,genBParts[0].phi,genBParts[0].mass)
                GenBJ2.SetPtEtaPhiM(genBParts[1].pt,genBParts[1].eta,genBParts[1].phi,genBParts[1].mass)
                self.out.fillBranch("GenBJ1_pt",genBParts[0].pt)
                self.out.fillBranch("GenBJ1_eta",genBParts[0].eta)
                self.out.fillBranch("GenBJ1_phi",genBParts[0].phi)
                self.out.fillBranch("GenBJ1_mass",genBParts[0].mass)
                self.out.fillBranch("GenBJ2_pt",genBParts[1].pt)
                self.out.fillBranch("GenBJ2_eta",genBParts[1].eta)
                self.out.fillBranch("GenBJ2_phi",genBParts[1].phi)
                self.out.fillBranch("GenBJ2_mass",genBParts[1].mass)
                self.out.fillBranch("GenBJJ_pt",(GenBJ1+GenBJ2).Pt())
                self.out.fillBranch("GenBJJ_eta",(GenBJ1+GenBJ2).Eta())
                self.out.fillBranch("GenBJJ_phi",(GenBJ1+GenBJ2).Phi())
                self.out.fillBranch("GenBJJ_mass",(GenBJ1+GenBJ2).M())
                self.out.fillBranch("GenBJJ_dEta",abs(genBParts[0].eta-genBParts[1].eta))
                self.out.fillBranch("GenBJJ_dR", deltaR(genBParts[0],genBParts[1]))
                self.out.fillBranch("GenBJJ_dPhi", deltaPhi(genBParts[0],genBParts[1]))
                return True
            else:
                self.out.fillBranch("GenBJ1_pt",-999)
                self.out.fillBranch("GenBJ1_eta",-999)
                self.out.fillBranch("GenBJ1_phi",-999)
                self.out.fillBranch("GenBJ1_mass",-999)
                self.out.fillBranch("GenBJ2_pt",-999)
                self.out.fillBranch("GenBJ2_eta",-999)
                self.out.fillBranch("GenBJ2_phi",-999)
                self.out.fillBranch("GenBJ2_mass",-999)
                self.out.fillBranch("GenBJJ_pt",-999)
                self.out.fillBranch("GenBJJ_eta",-999)
                self.out.fillBranch("GenBJJ_phi",-999)
                self.out.fillBranch("GenBJJ_mass",-999)
                self.out.fillBranch("GenBJJ_dEta",-999)
                self.out.fillBranch("GenBJJ_dR",-999)
                self.out.fillBranch("GenBJJ_dPhi",-999)
                return False


    def analyze(self, event):
        if self.get_gen_info(event) and self.genOnly:
            return True
        """process event, return True (go to next module) or False (fail, go to next event)"""
        resolved = -1; boosted = -1;
        jets = list(Collection(event, "Jet"))
        fatjets = list(Collection(event, "FatJet"))
        met = Object(event, "MET")
        
        zMuons,wMuons,zElectrons,wElectrons = self.select_leptons(event)

        allLeptons = zElectrons[:]
        allLeptons.extend(zMuons)
        allLeptons.extend(wElectrons)
        allLeptons.extend(wMuons)

        jetFilterFlags = [True]*len(jets)
        for lepton in allLeptons:
            jetInd = lepton.jetIdx
            if jetInd >= 0:
                jetFilterFlags[jetInd] = False
        self.out.fillBranch("Jet_lepFilter",jetFilterFlags)
        fatjetFilterFlags = [True]*len(fatjets)
        for fatjet in fatjets:
            fatjet.jetFilter = True
            for lepton in allLeptons:
               if deltaR(fatjet,lepton) < 0.8:
                  fatjetFilterFlags[fatjets.index(fatjet)] = False
        self.out.fillBranch("FatJet_lepFilter",fatjetFilterFlags)


        V_4vec,MET,Vtype,vLeptons = self.select_vboson(event,jets)
        if Vtype<0: return False
        lepton_ch = self.lepton_channel(Vtype)
        resolved, higgs_cand,n_add_jets,n_add_jets_un = self.select_higgs_jets(event,jets,lepton_ch,vLeptons)
        if resolved<0: 
            boosted, higgs_cand_boost, n_bjets_outside_fat_jet = self.select_fat_jets(fatjets,jets)
            if boosted<0: return False 
        print resolved, higgs_cand,n_add_jets,n_add_jets_un 
        if resolved>0:
          channel,stxs_cat = self.categorise_event(lepton_ch,higgs_cand,V_4vec,MET,n_add_jets,resolved,boosted) 
        elif boosted>0: 
          channel,stxs_cat = self.categorise_event(lepton_ch,higgs_cand_boost,V_4vec,MET,n_bjets_outside_fat_jet,resolved,boosted)
          print "STXS boosted = ", stxs_cat
        else: stxs_cat=-1;
        print "STXS = ", stxs_cat 
        if stxs_cat<0: return False 
        self.out.fillBranch("Reco_Cat",stxs_cat)
        return True
        
                

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

vhbb2016 = lambda : VHbbProducer(True,"2016") 
vhbb2017 = lambda : VHbbProducer(True,"2017") 
vhbb2018 = lambda : VHbbProducer(True,"2018",False) 
vhbb2018_gen = lambda : VHbbProducer(True,"2018",True,False,True) 
jmeCorrections2018MC = createJMECorrector(True, "2018", "A", "Merged", "AK4PFchs", False,splitJER=False)
jmeCorrections2018MCAll = createJMECorrector(True, "2018", "A", "All", "AK4PFchs", False,splitJER=True)
jmeCorrections2018AK8MC = createJMECorrector(True, "2018", "A", "Merged", "AK8PFPuppi", False,splitJER=False)
jmeCorrections2018AK8MCAll = createJMECorrector(True, "2018", "A", "All", "AK8PFPuppi", False,splitJER=True)
mhtVHbb = lambda : mhtProducer( lambda j : j.pt > 30,
                            lambda mu : mu.pt > 5 and mu.pfRelIso04_all < 0.4,
                            lambda el : el.pt > 5 and el.pfRelIso03_all < 0.4 )
