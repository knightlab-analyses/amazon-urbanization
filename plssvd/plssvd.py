# usage: python me.py metadata_fp microbes_fp metabolites_fp
# [microbes_label_fp] [metabolites_label_fp]

import sys
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import qiime2
from skbio import TreeNode
from gneiss.util import match
from sklearn.cross_decomposition import PLSSVD
from skbio.stats.composition import clr, centralize, multiplicative_replacement
from biplot import make_biplot

plt.rcParams['svg.fonttype'] = 'none'

args = sys.argv[1:]

mapping = pd.read_table(args[0], index_col=0, sep='\t')
microbes = qiime2.Artifact.load(args[1]).view(pd.DataFrame)
metabolites = qiime2.Artifact.load(args[2]).view(pd.DataFrame)

# do this match thing twice to make sure they are all matched
mapping, microbes = match(mapping, microbes)
mapping, metabolites = match(mapping, metabolites)
microbes, metabolites = match(microbes, metabolites)

mapping, microbes = match(mapping, microbes)
mapping, metabolites = match(mapping, metabolites)
microbes, metabolites = match(microbes, metabolites)

catdict = {i + 1: val for i,
           val in enumerate(sorted(mapping['category'].unique().tolist()))}

n = mapping.shape[0]
print('Number of samples: %d' % n)
print('Number of microbes: %d' % microbes.shape[1])
print('Number of metabolites: %d' % metabolites.shape[1])

# drop features represented in <10% samples
min_n = int(n * 0.1)

microbes = microbes.loc[:, microbes.mean(axis=0) > 20]
microbes = microbes.loc[:, (microbes > 0).sum(axis=0) > min_n]
microbes = microbes[(microbes.T != 0).any()]
print('Number of retained microbes: %d' % microbes.shape[1])

metabolites = metabolites.loc[:, metabolites.mean(axis=0) > 2e-3]
metabolites = metabolites.loc[:, (metabolites > 0).sum(axis=0) > min_n]
metabolites = metabolites[(metabolites.T != 0).any()]
print('Number of retained metabolites: %d' % metabolites.shape[1])

mapping, microbes = match(mapping, microbes)
mapping, metabolites = match(mapping, metabolites)
microbes, metabolites = match(microbes, metabolites)

mapping, microbes = match(mapping, microbes)
mapping, metabolites = match(mapping, metabolites)
microbes, metabolites = match(microbes, metabolites)
print('Number of retained samples: %d' % mapping.shape[0])

microbes.to_csv('microbes.csv')
metabolites.to_csv('metabolites.csv')
mapping.to_csv('mapping.csv')

# call external R script indval.R to perform statistics
p = subprocess.Popen('Rscript indval.R', shell=True, stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE)
p.communicate()

microbe_iv = pd.read_csv('microbe.indicator.values.csv', index_col=0)
metabolite_iv = pd.read_csv('metabolite.indicator.values.csv', index_col=0)

microbe_iv['group'] = microbe_iv['group'].map(catdict)
metabolite_iv['group'] = metabolite_iv['group'].map(catdict)

# highlight features with p-value <= 0.001
max_pval = 0.001

microbe_iv.loc[microbe_iv.pval > max_pval, 'group'] = 'None'
print('Number of significant microbes: %d'
      % microbe_iv[microbe_iv['group'] != 'None'].shape[0])

metabolite_iv.loc[metabolite_iv.pval > max_pval, 'group'] = 'None'
print('Number of significant metabolites: %d'
      % metabolite_iv[metabolite_iv['group'] != 'None'].shape[0])

plssvd = PLSSVD(n_components=3)
plssvd.fit(X=clr(centralize(multiplicative_replacement(microbes))),
           Y=clr(centralize(multiplicative_replacement(metabolites))))


def standardize(A):
    A = (A - np.mean(A, axis=0)) / np.std(A, axis=0)
    return A


pls_microbes = pd.DataFrame(standardize(plssvd.x_weights_),
                            columns=['PCA1', 'PCA2', 'PCA3'],
                            index=microbes.columns)
pls_metabolites = pd.DataFrame(standardize(plssvd.y_weights_),
                               columns=['PCA1', 'PCA2', 'PCA3'],
                               index=metabolites.columns)

color_map = {'Checherta': '#FF1600', 'Puerto Almendras': '#028001',
             'Iquitos': '#0122FF', 'Manaus low': '#F2BC04',
             'Manaus middle': '#F27305', 'None': '#C0C0C0'}
cat_order = ['Checherta', 'Puerto Almendras', 'Iquitos', 'Manaus low',
             'Manaus middle', 'None']
cat_order = {val: i + 1 for i, val in enumerate(cat_order)}

kwds = {}
if len(args) > 3:
    with open(args[3], 'r') as f:
        kwds['sample_labels'] = dict(x.split('\t')
                                     for x in f.read().splitlines())
if len(args) > 4:
    with open(args[4], 'r') as f:
        kwds['feature_labels'] = dict(x.split('\t')
                                      for x in f.read().splitlines())

fig = make_biplot(pls_microbes, pls_metabolites,
                  sample_metadata=microbe_iv, feature_metadata=metabolite_iv,
                  sample_color_category='group',
                  feature_color_category='group',
                  sample_color_dict=color_map, feature_color_dict=color_map,
                  sample_zorder=cat_order, feature_zorder=cat_order,
                  samp_alpha=0.75, feat_alpha=0.75,
                  feature_order=1,
                  **kwds)

fig[0].savefig('plot.svg')
