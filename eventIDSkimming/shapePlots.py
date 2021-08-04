import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json
from collections import OrderedDict


def replaceNan(arr, num):
  for i in range(len(arr)):
    if np.isnan(arr[i]) or np.isinf(arr[i]):
      arr[i] = num
  return arr

def hist(shape, reweights, i, options, density=False):
  n_sm, edges = np.histogram(shape, options.bin_number, options.bin_range, weights = reweights[:,0], density=density)
  n_eft, edges = np.histogram(shape, options.bin_number, options.bin_range, weights = reweights[:,(2+i*2)], density=density)
  return n_sm, n_eft, edges

def getRatioAndUncert(shape, reweights, i, options):
  #get ratio using all events
  n_sm, n_eft, edges = hist(shape, reweights, i, options)
  ratio = n_eft/n_sm -1

  overall_scaling = sum(n_eft) / sum(n_sm) - 1

  ratio = replaceNan(ratio, overall_scaling)
  #print(ratio)

  #separate events into 3 and calculate ratios
  ratios = []
  nsplits = 10
  indicies = np.arange(0, len(shape), 1)
  np.random.shuffle(indicies)
  splits = np.array_split(indicies, nsplits)

  for j in range(nsplits):
    n_sm, n_eft, edges = hist(shape[splits[j]], reweights[splits[j]], i, options)
    ratios.append(replaceNan(n_eft/n_sm -1, overall_scaling))
  ratios = np.array(ratios)
  #print(ratios)

  uncert = np.std(ratios, axis=0) / np.sqrt(nsplits)
  #print(uncert[25])
  return ratio, uncert, overall_scaling

def plot(shape, reweights, i, options, params, defval):
  n_sm_normed, n_eft_normed, edges = hist(shape, reweights, i, options, density=True)
  ratio, uncert, overall_scaling = getRatioAndUncert(shape, reweights, i, options)

  f, axs = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

  #print(len(edges))
  #print(edges)
  #print(n_sm_normed)
  #print(n_eft_normed)

  axs[0].hist(edges[:-1], edges, weights=n_sm_normed, histtype="step", label="SM")
  axs[0].hist(edges[:-1], edges, weights=n_eft_normed, histtype="step", label="%s = %.2f"%(params[i], defval))
  axs[0].set_ylabel("No. Events (normalised)")

  std_down = ratio-uncert
  std_up = ratio+uncert

  axs[1].hist(edges[:-1], edges, weights=ratio, histtype="step", color="tab:orange")
  axs[1].fill_between(edges, list(std_down)+[std_down[-1]], list(std_up)+[std_up[-1]], alpha=0.3, step="post", color="grey", lw=1)
  line = axs[1].plot([edges[0], edges[-1]], [overall_scaling, overall_scaling], 'g--', label="overall scaling")
  axs[1].set_xlabel(options.xlabel)
  axs[1].set_ylabel("Ratio to SM -1")
  axs[1].set_ylim(min(ratio-uncert)-max(uncert)*0.2, max(ratio+uncert)+max(uncert)*0.2)

  handles, labels = axs[0].get_legend_handles_labels()
  handles.append(line[0])
  axs[0].legend(handles=handles)

  f.savefig("%s/%s.png"%(options.savepath, params[i]))
  f.savefig("%s/%s.pdf"%(options.savepath, params[i]))


from optparse import OptionParser
parser = OptionParser(usage="%prog reweightedNanoAOD branch")

parser.add_option("--weight-branch", dest="weight_branch", default=None,
                  help="By default, standalone reweighting using the 'genWeight' branch as the original weights to scale. Use this option to specify a different branch to act as the original weight.")
parser.add_option("--config", "-c", dest="config", default=None,
                  help="The json config used to generate the reweighting module. If included, parameter names will be shown in plots.")
parser.add_option("--range", "-r", dest="bin_range", default="0,200",
                  help="Bin range for the histogram.")
parser.add_option("--binNumber", "-n", dest="bin_number", default=50, type="int",
                  help="Number of bins.")
parser.add_option("--xlabel", "-x", dest="xlabel", default=None,
                  help="Label for the x axis. Default will be the name of the branch")                 
parser.add_option("--save-path", "-s", dest="savepath", default="./",
                  help="Directory to save the plots to.")

(options, args) = parser.parse_args()

if len(args) < 2:
  parser.print_help()
  sys.exit(1)
filename = args[0]
branch = args[1]

f = uproot.open(filename)["Events"]

shape = f[branch].array(library='np')
reweights = f['Reweights'].array(library='np')

select = abs(reweights)[:,4].argsort()[0:int(0.999*len(shape))]
shape = shape[select]
reweights = reweights[select]

print(shape)

#if using different original weight, find scaling and apply to new original weights
if options.weight_branch is not None:
  original_weights = f[options.weight_branch].array(library='np')
  sf = reweights / reweights[:,0,np.newaxis]
  reweights = sf * original_weights[...,np.newaxis]

"""
select = reweights[:,4] < 10**5
shape = shape[select]
reweights = reweights[select]
"""

n_rw = len(reweights[0])
n_pars = int((-3 + np.sqrt(9+8*(n_rw-1)))/2)

if options.config is not None:
  with open(options.config, "r") as f:
    options.config = json.loads(f.read(), object_pairs_hook=OrderedDict)
  params = [entry['name'] for entry in options.config['parameters']]
  defval = options.config['parameter_defaults']['val']
else:
  params = [i for i in range(n_pars)]
  defval = 0.01

options.bin_range = [float(x) for x in options.bin_range.split(",")]

if options.xlabel is None:
  options.xlabel = branch

for i in range(n_pars):
  plot(shape, reweights, i, options, params, defval)
