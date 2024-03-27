import uproot
import numpy as np
import json
from collections import OrderedDict
import os
import sys
import bisect
import math
import pickle
import copy
import ROOT
import cmsstyle as CMS
CMS.SetExtraText("Simulation Preliminary")
CMS.SetLumi("")


ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPalette(56)
ROOT.gROOT.SetBatch(ROOT.kTRUE)

from optparse import OptionParser

def get_scaling(i_req, j_req, matrx, npars):
    k = 0
    for i in range(npars):
        for j in range(i+1):
            if (i==i_req and j==j_req):
                print(i_req, j_req, i, j, matrx[k])
                return matrx[k]
            k = k+1


parser = OptionParser(usage="%prog EFT2Obs_config.json inputFile1 inputFile2... [options] ")
parser.add_option("--parametrisation", dest="par_file", default="parametrisation.json",help = "")
parser.add_option("--datacard-path", dest="datacard", default="datacard.root",help = "path to combined datacard")
parser.add_option("--coeff", dest="coeff", default="datacard.root",help = "path to combined datacard")
parser.add_option("--category", dest="category", default="",help = "path to combined datacard")
parser.add_option("--process", dest="process", default="",help = "path to combined datacard")

_MYW = 800; _MYH = 600
_MYT = 0.05*_MYH; _MYB = 0.3*_MYH;
_MYL = 0.1*_MYW; _MYR = 0.05*_MYW


(options, args) = parser.parse_args()

with open(options.par_file, "r") as f:
    par_json = json.load(f)
    print(par_json[0])

datacard_root = ROOT.TFile(options.datacard, "read") 

color = {"chj3":ROOT.kBlue+2,"chj1":ROOT.kMagenta+7, "chw":ROOT.kMagenta+3}


# go through categories
for d in par_json:
    category = d["channel"]
    if options.category != "" and not options.category in category: continue
    process = d["process"]
    if options.process != "" and not options.process in process: continue
    if "FWDH" in process: continue 
    parameters = d["parameters"]
    n_pars = len(parameters)
    scaling = d["scaling"]
    n_bins = len(scaling)
    category_hist = datacard_root.Get(category + "/" +process+"125" ) # nominal hist
    print("extracting hist for category: ", category, "process: ",process)
    try: 
        nBins_datacard = category_hist.GetNbinsX()
    except:
        print("missing hist for category: ", category, "process: ",process)
        continue
    if (n_bins != category_hist.GetNbinsX()): 
        print("ERROR: bins from parametrisation file does not match the bins histograms\n Check your inputs")
    category_hist.SetLineColor(ROOT.kGray)
    category_hist.SetLineWidth(2)
    category_hist.SetTitle("")
    hist_modified = category_hist.Clone()
    hist_modified.SetLineColor(color[options.coeff.split("=")[0]])
    hist_modified.SetLineWidth(2)
    for ibin, matrx in enumerate(scaling):
        if (len(matrx) != (n_pars*(n_pars+1)/2)): 
            print("ERROR: unexpected length of scaling array")
        ipar = parameters.index(options.coeff.split("=")[0] + "[0,-1,1]" ) 
        print(options.coeff.split("=")[1])
        val = float(options.coeff.split("=")[1])
        print(n_pars)
        linear = get_scaling(n_pars-1, ipar, matrx, n_pars)
        quadr = get_scaling(ipar, ipar, matrx, n_pars)   
        print(linear,quadr)
        mu_ = 1 + linear*val + quadr*val*val # scale current bin with 2nd order equation
        print("channel = ", category, "; process = ", process, "; bin = ", ibin, "; linear = ", linear, "; quadr = ",quadr)
        hist_modified.SetBinContent(ibin + 1, mu_*category_hist.GetBinContent(ibin + 1))
    maximum = 1.4*max(x.GetMaximum() for x in [category_hist, hist_modified])
    hist_modified.SetAxisRange(0,maximum,'Y')
    cB=0;
    cB=ROOT.TCanvas("cB","cB",_MYW,_MYH);
    
    rp = ROOT.TRatioPlot(hist_modified,category_hist,"divsym");
    rp.SetH2DrawOpt("hist");rp.SetH1DrawOpt("hist");
    rp.GetLowYaxis().SetNdivisions(505)
    rp.SetLeftMargin( _MYL/_MYW );  rp.SetRightMargin( _MYR/_MYW );
    rp.SetUpTopMargin( _MYT/_MYH );   rp.SetLowBottomMargin( _MYB/_MYH );
    rp.Draw("nogrid");
    rp.GetLowerRefYaxis().SetTitle("Ratio to SM");

    leg1 = ROOT.TLegend(0.12, 0.8, 0.3, 0.95);
    leg1.SetTextFont(42)
    leg1.SetTextSize(0.04)
    leg1.SetBorderSize(0)
    leg1.AddEntry(0,"Category: "+category,"")
    leg1.AddEntry(0,"STXS: "+process,"")
    leg1.SetFillStyle(0);
    leg1.Draw('same');

    leg = ROOT.TLegend(0.12, 0.65, 0.3, 0.8);
    leg.SetFillStyle(0);
    leg.SetTextFont(42)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.05)
    leg.AddEntry(category_hist,"c_{i}=0 (SM)", "l");
    leg.AddEntry(hist_modified,options.coeff, "l");
    leg.Draw('same');

    rp.GetLowerPad().cd()
    line = ROOT.TLine(0,1,1,1)
    #line.SetLineStyle(kDashed)
    line.SetLineColor(ROOT.kGray)
    line.SetLineWidth(2)
    line.Draw("same")

    rp.GetLowerRefGraph().SetMinimum(0.)
    rp.GetLowerRefGraph().SetMaximum(3.)
    cB.Update()
    cB.SaveAs(process+"_"+category+"_"+options.coeff.split("=")[0]+".pdf")#plot_var+'_EFT_allCats_'+fine+'stxs'+str(options.stxs)+"_"+pth_str+'.pdf');

    


                
                












