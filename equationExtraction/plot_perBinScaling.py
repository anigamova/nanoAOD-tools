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

from optparse import OptionParser

def get_scaling(i_req, j_req, matrx, npars):
    k = 0
    for i in range(n_pars):
        for j in range(i+1):
            if (i==i_req and j==j_req):
                return matrx[k]
            k = k+1


parser = OptionParser(usage="%prog EFT2Obs_config.json inputFile1 inputFile2... [options] ")
parser.add_option("--parametrisation", dest="par_file", default="parametrisation.json",help = "")
parser.add_option("--datacard-path", dest="datacard", default="datacard.root",help = "path to combined datacard")
parser.add_option("--coeff", dest="coeff", default="datacard.root",help = "path to combined datacard")



(options, args) = parser.parse_args()

with open(options.par_file, "r") as f:
    par_json = json.load(f)
    print(par_json[0])

datacard_root = ROOT.TFile(options.datacard, "read") 

# go through categories
for d in par_json:
    category = d["channel"]
    process = d["process"]
    if "FWDH" in process: continue 
    parameters = d["parameters"]
    n_pars = len(parameters)
    scaling = d["scaling"]
    n_bins = len(scaling)
    category_hist = datacard_root.Get(category + "/" +process+"125" ) # nominal hist
    if (n_bins != category_hist.GetNbinsX()): 
        print("ERROR: bins from parametrisation file does not match number histograms\n Check your inputs")
    for ibin, matrx in enumerate(scaling):
        if (len(matrx) != (n_pars*(n_pars+1)/2)): 
            print("ERROR: unexpected length of scaling array")
        ipar = parameters.index(options.coeff[0])
        value = options.coeff[1]
        linear = get_scaling(ipar, npar, matrx, npars)
        quadr = get_scaling(ipar, ipar, matrx, npars)
        mu_ = 1 + linear*val + quadr*val**2
        category_hist.SetBinContent(ibin + 1, mu_*category_hist.GetBinContent(ibin + 1))

                
                












