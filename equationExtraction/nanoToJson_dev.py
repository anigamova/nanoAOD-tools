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

def tqdm(iterable):
  for i, item in enumerate(iterable):
    j = i + 1
    sys.stdout.write('\r'+"%d/%d"%(j, len(iterable)))
    sys.stdout.flush()
    if j == len(iterable): print("")
    yield item


def getTerms_from_Trweights(rw, j, x):
  """
  Derive A and B coefficients for scaling equation.
  sm: array of event weights at SM benchmark
  rw1: array of event reweights at C=x/2
  rw2: array of event reweights at C=x
  """
  sm = rw[:,0].sum()
  #print('val,ratio,index',x,rw[:,j].sum()/sm,j)
  if (j % 2) != 0: mu1 = rw[:,j].sum() / (sm)
  #if (j % 2) != 0: mu1 = rw[:,j].sum() / (x*sm)
  else: mu1 = rw[:,j].sum() / (sm)
  #else: mu1 = rw[:,j].sum() / (x*x*sm)
  return mu1

def getCrossTerm_from_Trweights(rw,  ic, x):
  """
  Derive B_ij cross term for scaling equation.
  sm: array of event weights at SM benchmark
  rw: array of event reweights at C1=x C2=x
  A1, B1, A2, B2: linear and quad scaling eqn coefficients for C1 and C2
  """
  sm = rw[:,0].sum()
  mu = rw[:,ic].sum() / (sm)#*x[0]*x[1])
  return mu



def fill_term(val, ind, unc, n, rw, def_val):
  term = {}
  term["val"] = val
  term["ind"] = ind
  term["unc"] = unc
  term["N_Events"] = n
  term["rws"] = rw #.tolist()
  term["def_val"] = def_val #.tolist()
  return term


def processConfig(config_path):
  with open(config_path, "r") as f:
    config_dict = json.load(f)

#  def_vals = config_dict['parameter_defaults']['val']
  params = [param['name'] for param in config_dict['parameters']]
  def_vals = [param['val'] for param in config_dict['parameters']]
  return def_vals, params 

#helper function for reading json and converting keys to integers
def keysToInt(x):
  od = OrderedDict()
  for k, v in x:
    od[int(k)] = v
  return od

def getKeys(options):
  key, key2, key3 = None, None, None

  if options.key != None:
    with open(options.key, "r") as f:
      key = json.load(f, object_pairs_hook=keysToInt)
  if options.key2 != None:
    with open(options.key2, "r") as f:
      key2 = json.load(f, object_pairs_hook=keysToInt)

  if options.key3 != None:
    with open(options.key3, "r") as f:
      key3 = json.load(f, object_pairs_hook=keysToInt)

  return key, key2, key3

def vhbb_cond(key):

  return "_9_" in key and ("Wmn" in key or "Wen" in key)

def define_cat_dnn_bin(dnn_vals, category, key, options):
  with open(options.obs_binning, "r") as f:
    binning_json = json.load(f)
    #print(category)
    #vhbb_debug_cond = "_9_" in key[x] and ("Wmn" in key[x] or "Wen" in key[x])
    bins = np.array([binning_json[key[x]] if not vhbb_cond(key[x]) else binning_json[key[x].replace("_9_","_5_")] for x in category ])
    binind = np.array([ (bisect.bisect_left(bins[i],x)-1) for i,x in enumerate(dnn_vals)])
    return binind

def flatten_arrays(first, second, key=None, key2=None):
  max_tag_2 = max(second)
  min_tag_2 = min(second)
  spacing = max_tag_2 - min_tag_2 + 1
  new_tag = (first*spacing + (second-min_tag_2)).astype(int)
  if key is not None and key2 is not None:
    new_key = OrderedDict()
    for i in np.unique(first):
      for j in np.unique(second):
        tag = (i*spacing) + (j-min_tag_2)
        new_key[tag] = "%s__%s"%(key[i], key2[j])
    return new_tag, new_key
  else: 
    return new_tag



def cat_ind(options,f):
  year = '2018'#TODO make configurable
  dictionary = {
    'htt':{
      'cats' : [200,201,202,203,100,101,102,103,104,105,106,107,108,109,110],
      'channels' : ['em','et','mt','tt'] 
    },
    'vhbb':{
      'cats': [1,5,9,21,22,23,24],
      'channels' : ['Znn','Wen','Wmn','Zee','Zmm'] 
    },
    'ttH':{
      'ljets_ge6j_ge4t_classification_ttH_ttmb_vs_slike_times_ttHbb_STXS_0':6,
      'ljets_ge6j_ge4t_classification_ttH_ttmb_vs_slike_times_ttHbb_STXS_1':7,
      'ljets_ge6j_ge4t_classification_ttH_ttmb_vs_slike_times_ttHbb_STXS_2':8,
      'ljets_ge6j_ge4t_classification_ttH_ttmb_vs_slike_times_ttHbb_STXS_3':9,
      'ljets_ge6j_ge4t_classification_ttH_ttmb_vs_slike_times_ttHbb_STXS_4':10,
      'ljets_5j_ge4t_classification_ttH_ttmb_vs_slike_times_ttHbb_STXS_0':1,
      'ljets_5j_ge4t_classification_ttH_ttmb_vs_slike_times_ttHbb_STXS_1':2,
      'ljets_5j_ge4t_classification_ttH_ttmb_vs_slike_times_ttHbb_STXS_2':3,
      'ljets_5j_ge4t_classification_ttH_ttmb_vs_slike_times_ttHbb_STXS_3':4,
      'ljets_5j_ge4t_classification_ttH_ttmb_vs_slike_times_ttHbb_STXS_4':5
    }

  }
  if options.analysis == 'vhbb': 
    channel = 1*(1-(f['isZee'].array() + f['isWmunu'].array() + f['isWenu'].array() + f['isZmm'].array())) + 2*f['isWmunu'].array() + 3*f['isWenu'].array() + 4*f['isZee'].array() + 5*f['isZmm'].array()
  elif options.analysis == 'htt': 
    channel = f['channel'].array()

  if options.analysis != 'ttH':
    if options.analysis != 'vhbb':cats_key = {i:str(k)+'_'+year for i,k in enumerate(dictionary[options.analysis]['cats'])}
    else: cats_key = {i:str(k)+'_13TeV'+year for i,k in enumerate(dictionary[options.analysis]['cats'])}
    chan_key = {(i+1):(options.analysis+'_'+k) for i,k in enumerate(dictionary[options.analysis]['channels'])}
    category =  f['category'].array()
    category_id = np.array([dictionary[options.analysis]['cats'].index(x) for x in f['category'].array()])
    chan_name = np.array([dictionary[options.analysis]['channels'][int(x)-1] for x in channel])
    global_cat_index,cat_key = flatten_arrays(channel,category_id,chan_key,cats_key)
    cat_key = {i:k.replace('__','_') for i,k in cat_key.items()}
  else: 
    global_cat_index = f['category'].array()
    cat_key = {v: 'ttH_'+year+'_'+k for k, v in dictionary['ttH'].items()}
  print(global_cat_index)
  return global_cat_index,cat_key

def readNanoFile(nano_file, options, key, key2, key3):
  f = uproot.open(nano_file)["Events"]
  reweights_original = f["Reweights"].array()
  ratio_to_sm  = (reweights_original/reweights_original[:,0][:,np.newaxis])<50
  cut = np.all(ratio_to_sm,axis=1)
  
  reweights = f["Reweights_transformed"].array()
  reweights = reweights[cut]
  if options.weight_branch != None:
    original_weights = f[options.weight_branch].array()
    sf = reweights / reweights[:,0,np.newaxis]
    reweights = sf * original_weights[...,np.newaxis]

  if options.inclusive:
    return reweights, np.ones(len(reweights)), {1.0: "inclusive"}

  tag_branch = f[options.tag_branch].array()
  tag_branch = tag_branch[cut]
  if options.use_fine_stxs and options.analysis=='vhbb':  #TODO:move to a separate method
    stxs_fine = f[options.tag_branch_fine].array()
    stxs_fine = stxs_fine[cut]
    tag_branch = np.where(tag_branch%100==4, 1, tag_branch) 
    tag_branch = np.where(tag_branch%100<5, tag_branch, stxs_fine)

    tag_branch = np.where(tag_branch%100==4, (tag_branch/100)*100 + 27, tag_branch)
    tag_branch = np.where(tag_branch%100==9, (tag_branch/100)*100 + 27, tag_branch)
    tag_branch = np.where(tag_branch%100==14, (tag_branch/100)*100 + 27, tag_branch)

    tag_branch = np.where(tag_branch%100==5, (tag_branch/100)*100 + 30, tag_branch)
    tag_branch = np.where(tag_branch%100==10, (tag_branch/100)*100 + 30, tag_branch)
    tag_branch = np.where(tag_branch%100==15, (tag_branch/100)*100 + 30, tag_branch)

    tag_branch = np.where(tag_branch==1, 304, tag_branch)

    
  if key == None:
    key = {}
    for i in np.unique(tag_branch):
      key[i] = str(i)

  if options.stage0 != None:
    #rebin to stage0
    tags = [k for k in key.keys() if ((options.stage0 in key[k]) and ("FWDH" not in key[k]))] #all keys corresponding to stage0 process which are not FWDH
    selection = np.isin(tag_branch, tags)
    return reweights[selection], np.ones(len(reweights[selection])), {1.0: options.stage0}

  if options.tag_branch2 == None and options.tag_branch3 == None:
    return reweights, tag_branch, key
  else:
    tag_branch2,key2 = cat_ind(options, f)
    tag_branch2 = tag_branch2[cut]
    if key2 is None:
      key2 = {}
      for i in np.unique(tag_branch2):
        i = int(i)
        key2[i] = str(i)
    new_key12 = OrderedDict()
    new_tag_branch_12,new_key12 = flatten_arrays(tag_branch,tag_branch2, key, key2)
    
    if options.tag_branch3 != None:

      tag_branch3 = define_cat_dnn_bin(f[options.tag_branch3].array()[cut], tag_branch2, key2, options)
      if key3 == None:
        key3 = {}
        for i in np.unique(tag_branch3):
          key3[i] = 'bin'+str(i)
      new_tag_branch,new_key = flatten_arrays(new_tag_branch_12,tag_branch3,new_key12,key3)

    else: 
      new_tag_branch = new_tag_branch_12
      new_key = new_key12
    dev = True
    #print("level one(stxs)",tag_branch)
    #print("level two(stxs+cat)",new_tag_branch_12)
    #print("level three(stxs+cat+bin)",new_tag_branch)

    if not dev: return reweights, new_tag_branch, new_key  
    else: return reweights, tag_branch, key, new_tag_branch_12, new_key12, new_tag_branch, new_key

def deriveEqns_alllevels(params, def_val, reweights,tag1,tag2,tag3,key1,key2,key3,options):
  new_eqn = OrderedDict()
  for j in range(len(params)):
    for tag_level in tag1,tag2,tag3:
      for k in np.unique(tag_level):
        if k in key3.keys():
          key = key3[k]
        elif k in key2.keys(): 
          key = key2[k]
        elif k in key1.keys():
          key = key1[k]
        else: print("error")
        rw = reweights[tag_level==k]
        n_events = len(rw)
        #print("key:", key, "; nevents=", n_events)
        #if key=="QQ2HLNU_FWDH__vhbb_Wmn_24_13TeV2018": print("QQ2HLNU_FWDH__vhbb_Wmn_24_13TeV2018 nweights = ", n_events)#debug
        resamples = np.random.randint(0, n_events, size=(options.n_bootstrap,n_events))

        A = getTerms_from_Trweights(rw, 2*j+1, def_val[j])
        B = getTerms_from_Trweights(rw, 2*j+2, def_val[j])
        A_values = [getTerms_from_Trweights(rw[s], 2*j+1, def_val[j]) for s in resamples]
        B_values = [getTerms_from_Trweights(rw[s], 2*j+2, def_val[j]) for s in resamples]
        
        term_lin = {};term_quadr = {}
        
        term_lin = fill_term(float(A),2*j+1,float(np.std(A_values)),n_events,rw,[def_val[j]])
        term_quadr = fill_term(float(B),2*j+2,float(np.std(B_values)),n_events,rw,[def_val[j]])
        
        if not "A_%s"%params[j] in new_eqn.keys(): new_eqn["A_%s"%params[j]] = {key:term_lin}
        else: new_eqn["A_%s"%params[j]].update({key:term_lin})
        if not "B_%s_2"%params[j] in new_eqn.keys(): new_eqn["B_%s_2"%params[j]] = {key:term_quadr}
        else: new_eqn["B_%s_2"%params[j]].update({key:term_quadr})
        
  ic = 0
  for j in range(len(params)):
    for k in range(j+1, len(params)):
      for tag_level in tag1,tag2,tag3:
        for tag in np.unique(tag_level):
          rw = reweights[tag_level==tag]
          n_events = len(rw)
          resamples = np.random.randint(0, n_events, size=(options.n_bootstrap,n_events))

          if tag in key3.keys():
            key = key3[tag]
          elif tag in key2.keys():
            key = key2[tag]
          elif tag in key1.keys():
            key = key1[tag]
          else: print("error")

          #print("key:", key, "; nevents=", n_events)
          #if key=="QQ2HLNU_FWDH__vhbb_Wmn_24_13TeV2018": print("QQ2HLNU_FWDH__vhbb_Wmn_24_13TeV2018 nweights = ", n_events)#debug

          B = getCrossTerm_from_Trweights(rw, 2*len(params)+ 1 + ic,[def_val[j],def_val[k]])
          B_values = [getCrossTerm_from_Trweights(rw[s], 2*len(params)+ 1 + ic, [def_val[j],def_val[k]]) for s in resamples]

          term = fill_term(float(B),2*len(params)+ 1 + ic,float(np.std(B_values)),n_events,rw,[def_val[j],def_val[k]])
          if "B_%s_%s"%(params[j], params[k]) in new_eqn.keys():
             new_eqn["B_%s_%s"%(params[j], params[k])].update({key:term})
          else: new_eqn["B_%s_%s"%(params[j], params[k])] = {key:term}

      ic += 1
  return new_eqn;

def make_stxs_cat_bin_dict(key1,key2,key3,eqns):
#  print
  term = eqns[eqns.keys()[0]]
  stxs = [i for i in key1.values() if 'QQ2HLNU' in i ]
  #stxs = [i for i in key1.values() if "QQ2HLL" in i or 'QQ2HLNU' in i ]
  stxs_cat = key2.values()
  stxs_cat_bin = key3.values()
  out = {}
  for i in stxs:
    if i not in term.keys(): continue;
    out[i] = {}
    for j in stxs_cat:
      if j not in term.keys(): continue;
      bins = []
      for k in stxs_cat_bin:
        if j in k and k in term.keys(): bins.append(k)
      if i in j: 
        out[i].update({j:bins})
  return out


def cleanUp(eqns, options):
  for tag in eqns.keys():
    for bin in eqns[tag].keys():
      params = eqns[tag][bin].keys()
      params = filter(lambda x: (x[0] != 'u' and x[0]!='b'), params)
      for param in params:
        if abs(eqns[tag][bin][param]) < options.absolute_threshold:
          del eqns[tag][bin][param]
          del eqns[tag][bin]["u_"+param]
  return eqns

def check_and_merge_eqns(eqns, tags):
  bins_to_merge = {}
  merge_flag = False
  eqns_copy = eqns.copy()
  for param in eqns_copy.keys():
    print("param = ", param)
    term = eqns_copy[param]
    for i in tags.keys():
      for j in tags[i].keys():
        stxs_eq_cat = False
        bins = []
        stxs_val=term[i]['val']
        cat_val=term[j]['val']
        cat_unc=term[j]['unc']
        if abs(stxs_val-cat_val)<cat_unc : stxs_eq_cat = True
        bins = list(tags[i][j])
        merged_prev = False
        found_different_term = False
        nbins = len(bins)

        for k in range(nbins-1):
          if k+1>len(bins)-1:break;
          this = term[bins[k]]['val']
          this_unc = term[bins[k]]['unc']
          param_ind = term[bins[k]]['ind']  
          param_val = term[bins[k]]['def_val']  
          next = term[bins[k+1]]['val']
          next_unc = term[bins[k+1]]['unc']
          total_unc = math.sqrt(this_unc**2 + next_unc**2)
          if abs(this-next)>total_unc: 
            merged_prev = False 
            if abs(this - cat_val)>this_unc:#check if compatible with cat term from cat within unc  
              found_different_term = True;
            else: del eqns_copy[param][bins[k]]  
          else: 
            bin1 = bins[k].split('__')[-1];bin2 = bins[k+1].split('__')[-1];
            bin_new = j+'__'+bin1+bin2
            rw1 = term[bins[k]]['rws']
            rw2 = term[bins[k+1]]['rws']
            rws_new = np.vstack((rw1, rw2))
            n_events = len(rws_new)
            resamples = np.random.randint(0, n_events, size = (options.n_bootstrap,n_events))
            if len(param_val)==1:
              A = getTerms_from_Trweights(rws_new, param_ind, param_val[0])
              A_values = [getTerms_from_Trweights(rws_new[s], param_ind, param_val[0]) for s in resamples]
            else: 
              A = getCrossTerm_from_Trweights(rws_new, param_ind, param_val)
              A_values = [getCrossTerm_from_Trweights(rws_new[s], param_ind, param_val) for s in resamples]

            del eqns_copy[param][bins[k]]
            del eqns_copy[param][bins[k+1]]

            if abs(A - cat_val)>this_unc:#check if compatible with cat term from cat within unc  
              found_different_term = True;
              eqns_copy[param][bin_new] = fill_term(float(A),param_ind,float(np.std(A_values)),n_events,rws_new,param_val)# TODO (not sure)   
              bins.insert(k+1,bin_new)
              bins.pop(k+2)

            else: 
              bins.pop(k+1)#ignore this bin and the next, don't add the merged to dict and remove from bins list
            merged_prev=True
  #after bin merging is completed no need to store rwgt array
  for par in eqns_copy.keys():
    for bin in eqns_copy[par].keys():
      del eqns_copy[par][bin]['rws']

  return eqns_copy    

def reshape_eqns_dict_forIntModel(eqns, tag_dict, options, pars, bin_dict):
  reshaped_eqns = []
  pars.append("sm")
  for key1 in tag_dict.keys():
    tag = key1 
    process = key1
    for key2 in tag_dict[key1].keys():
      d = OrderedDict();
      category = key2.split('__')[-1]
      #take into account the missing bins.
      bins_from_datacards = bin_dict[category]
      if len(tag_dict[key1][key2])<len(bins_from_datacards):
        print("bin numbers don't match, datacard length = ", len(bins_from_datacards), "par: = ", len(tag_dict[key1][key2]))
      bin_scaling = []
      for b in range(len(bins_from_datacards)-1):
        tag = ""
        scaling = []
        for i in range(len(pars)):
          for j in range(i+1):
            par_i = pars[i]
            par_j = pars[j]
            if i == j: term = "B_%s_2"%(par_i)
            elif j < i and par_i != 'sm': 
              term = "B_%s_%s"%(par_i,par_j) 
              if not term in eqns.keys(): 
                term = "B_%s_%s"%(par_j,par_i)
                if not term in eqns.keys(): print('term is missing')
            elif i==len(pars)-1: 
              term = "A_%s"%(par_j) 
            else: print("invalid term")
            print(term)
            if term == 'B_sm_2': s = 1.
            else: 
              tag = define_tag(eqns[term] ,key1, key2, b, bins_from_datacards) #TODO
              s = eqns[term][tag]["val"]
            if s < options.absolute_threshold: s = 0
            scaling.append(s)
        bin_scaling.append(scaling)             
      if (len(bin_scaling)) != (len(bins_from_datacards)-1):
        print("final bin numbers don't match, datacard length = ", len(bins_from_datacards), "par: = ", len(bin_scaling))
      else: print("final bin numbers are fine")
      d["channel"] = category
      d["process"] = process+'_'+options.final_state
      d["parameters"] = ["%s[0,-1,1]"%pars[k] if pars[k] != 'sm' else "%s[1]"%pars[k] for k in range(len(pars))]
      d["scaling"] = bin_scaling
      reshaped_eqns.append(d) 
  return reshaped_eqns

def define_tag(term, stxs, stxs_cat, b, bins_from_datacards):
  found_merged_bin_par = False
  found_bin_par = False
  full_tag = stxs_cat + "__bin" + str(b)
  #print("full tag from define_tag = ", full_tag)
  if full_tag in term.keys(): #found not merged bin
    tag = full_tag
    found_bin_par = True
  else:
    tmp = [i for i in term.keys() if stxs_cat in i and (i.endswith(full_tag.split('__')[-1]) or full_tag.split('__')[-1] + 'bin' in i)]
    if len(tmp)==1:
      tag = tmp[0]
      found_bin_par = True
      found_merged_bin_par = True
    elif len(tmp)==0: found_bin_par = False
    else: 
      print('ERROR', b,tmp)
  if (found_bin_par == False and found_merged_bin_par == False):
    if stxs_cat in term.keys():
      tag = stxs_cat
    else: tag = stxs
  return tag
  

def getOptions():
  from optparse import OptionParser
  parser = OptionParser(usage="%prog EFT2Obs_config.json inputFile1 inputFile2... [options] ")
  parser.add_option("--tag-branch", dest="tag_branch", default="HTXS_stage1_2_cat_pTjet30GeV",
                    help="Specify the branch used to tag the events, e.g. HTXS_stage1_2_cat_pTjet30GeV")
  parser.add_option("--tag-branch-fine", dest="tag_branch_fine", default="HTXS_stage1_2_fine_cat_pTjet30GeV",
                    help="Specify the branch used to tag the events, e.g. HTXS_stage1_2_cat_pTjet30GeV")

  parser.add_option("--tag-branch2", dest="tag_branch2", default=None,
                    help="Specify a second tag branch so the events will be binned by branch1 x branch2")
  parser.add_option("--tag-branch3", dest="tag_branch3", default=None,
                    help="Specify a second tag branch so the events will be binned by branch1 x branch2 x branch3")

  parser.add_option("--inclusive", dest="inclusive", default=False, action="store_true",
                    help="Dump all events into a single bin. Do no splitting by e.g. STXS bin.")
  parser.add_option("--fine-stxs", dest="use_fine_stxs", default=False, action="store_true",
                    help="use HTXS_stage1_2_fine_cat_pTjet30GeV to modify the defaule HTXS_stage1_2_cat_pTjet30GeV bins")

  parser.add_option("--weight-branch", dest="weight_branch", default=None,
                    help="By default, standalone reweighting using the 'genWeight' branch as the original weights to scale. Use this option to specify a different branch to act as the original weight.")
  
  parser.add_option("--key", dest="key", default="equationExtraction/tagKeys/HTXS_stage1_2_cat_pTjet30GeV.json",
                    help="tagKey.json file that provides names for each bin in tag_branch.")  
  parser.add_option("--key2", dest="key2", default=None,
                    help="tagKey.json file that provides names for each bin in tag_branch2.") 
  parser.add_option("--key3", dest="key3", default=None,
                    help="tagKey.json file that provides names for each bin in tag_branch3.")

  parser.add_option("--obs-binning", dest="obs_binning", default=None,
                    help="tagKey.json file that provides binning for each analysis category")

  parser.add_option("--outputDir", dest="outputDir", default="",
                    help="Specify the directory to save the output json files.")

  parser.add_option("--absoluteThreshold", dest="absolute_threshold", type=float, default=0,
                    help="Remove terms that are smaller than this number")
  parser.add_option("--minNevents", dest="nevents_threshold", type=float, default=0,
                    help="Remove terms that are smaller than this number")

  parser.add_option("--stage0", dest="stage0", default=None)
  parser.add_option("--filter-params", dest="filter_params", default="")
  parser.add_option("--analysis", dest="analysis", default="htt")
  parser.add_option("--final-state", dest="final_state", default="htt")

  parser.add_option("--nBootstrap", dest="n_bootstrap", type="int", default=5, 
                    help="The number of resamples used in the boostrapping")

  (options, args) = parser.parse_args()

  if len(args) < 2:
    parser.print_help()
    sys.exit(1)
  EFT2Obs_config = args[0]
  nano_files = args[1:]

  return EFT2Obs_config, nano_files, options

def main(EFT2Obs_config, nano_files, options):
  def_val, params = processConfig(EFT2Obs_config)
  key1, key2, key3 = getKeys(options)

  for nano_file in nano_files:
    nano_name = nano_file.split("/")[-1].split(".")[-2]
    json_path = os.path.join(options.outputDir, nano_name+".json")
    with open(options.obs_binning, "r") as f:
      binning_json = json.load(f)

    print(">> Reading %s and outputting %s"%(nano_file, json_path))
    reweights, tag_branch_lvl1, key_lvl1, tag_branch_lvl2, key_lvl2, tag_branch_lvl3, key_lvl3 = readNanoFile(nano_file, options, key1, key2, key3)
    #print("key2's",key_lvl2)
    print(">> Finished reading")
    eqns = deriveEqns_alllevels(params, def_val, reweights, tag_branch_lvl1, tag_branch_lvl2, tag_branch_lvl3, key_lvl1, key_lvl2, key_lvl3, options)
    #print("raw eqns: ",eqns)
    stxs_cat_bin_dict = make_stxs_cat_bin_dict(key_lvl1, key_lvl2, key_lvl3, eqns)
    #print(stxs_cat_bin_dict)
    merged = check_and_merge_eqns(eqns, stxs_cat_bin_dict)
    
    with open(json_path.replace('.json', "_merged.json"), "w") as f:
      json.dump(merged, f, indent = 4)
    reshaped = reshape_eqns_dict_forIntModel(merged, stxs_cat_bin_dict, options, params, binning_json)  
    with open(json_path.replace('.json', "_reshaped.json"), "w") as f:
      json.dump(reshaped, f, indent = 4)


if __name__=="__main__":
  EFT2Obs_config, nano_files, options = getOptions()
  main(EFT2Obs_config, nano_files, options)

