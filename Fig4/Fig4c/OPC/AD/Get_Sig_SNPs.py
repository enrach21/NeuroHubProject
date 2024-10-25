# Author: Ian Jones
# Goal: Given GKM_Predictions return significant snps

# Load libs

import sys
import math
import random
import statistics
from decimal import Decimal

# Data / ML / Stats Libraries

import numpy as np
import pandas as pd
import scipy
from scipy.stats import *

from plotnine import *
from matplotlib import pyplot as plt
# from viz_preprocess import *
# from viz_sequence import *
import warnings
warnings.filterwarnings('ignore')
plt.style.use('default')


# Read inputs:

SNPS = '/shen/shenlabstore3/ijones1/GKM_explain_test/MG_AD_Snps/Accessible.AD.8.1.24.csv'

# GKM input
GKM_ref = 'MG.AD.Filt.allATAC.SNPs.200bp.V2.txt'
GKM_alt = 'MG.AD.Filt.allATAC.ALT.SNPs.200bp.V2.txt'
GKM_ref_shuf = 'Ref_Shuffled.MG.AD.Filt.allATAC.SNPs.200bp.V2.txt'
GKM_alt_shuf = 'Alt_Shuffled.MG.AD.Filt.allATAC.SNPs.200bp.V2.txt'

# Read in delta SVM inputs
delta = "Delta_AD.Filt.MG.ATAC.200bp.txt"
delta_shuf = "Delta_Shuf_AD.Filt.MG.ATAC.200bp.txt"

# Read in delta SVM inputs
ISM_ref = 'ISM.MG.AD.Filt.allATAC.SNPs.200bp.V2.txt'
ISM_alt = 'ISM.MG.AD.Filt.allATAC.ALT.SNPs.200bp.V2.txt'
ISM_ref_shuf = 'ISM.Ref_Shuffled.MG.AD.Filt.allATAC.SNPs.200bp.V2.txt'
ISM_alt_shuf = 'ISM.Alt_Shuffled.MG.AD.Filt.allATAC.SNPs.200bp.V2.txt'

# Output csv
output ='OPC.AD.8.14.24.csv'

#### Read in SNPs
df = pd.read_csv(SNPS)
print(df.shape)
print(df.head)
df


#### Read in Explain scores

#read in the importance scores

impscores_ref = [
    np.array( [[float(z) for z in y.split(",")] for y in x.rstrip().split("\t")[2].split(";")])
    for x in open(GKM_ref)
]

impscores_alt = [
    np.array( [[float(z) for z in y.split(",")] for y in x.rstrip().split("\t")[2].split(";")])
    for x in open(GKM_alt)
]

shuf_impscores_ref = [
    np.array( [[float(z) for z in y.split(",")] for y in x.rstrip().split("\t")[2].split(";")])
    for x in open(GKM_ref_shuf)
]

shuf_impscores_alt = [
    np.array( [[float(z) for z in y.split(",")] for y in x.rstrip().split("\t")[2].split(";")])
    for x in open(GKM_alt_shuf)
]

#### Get Observed Magnitude and Prominence Scores for postive reigons

observed_confidence_dict = {'observed_seqlet_start': [], 'observed_seqlet_end': [],
                            'observed_seqlet_effect': [], 'observed_seqlet_noneffect': [],
                            'observed_bg_effect': [], 'observed_bg_noneffect': [],
                            'observed_magnitude_score': [], 'observed_prominence_score': [],
                            'observed_active_allele': [], 'observed_inactive_allele': [],
                            'observed_seqlet_len': []}

score = []
for x in np.arange(len(shuf_impscores_ref)):
    # get tested element
    y = shuf_impscores_ref[x]
    # get score of the bps
    for z in np.arange(len(y)):
        score.append(y[z][np.nonzero(y[z])][0])

null_upper_thresh = np.quantile(score, 0.975)
print('Per-base Importance Score Threshold: ', '\t', null_upper_thresh)

for seq in range(len(impscores_ref)):
    observed_per_seq_dict = {'observed_scores_active': [], 'observed_scores_inactive': [],
                             'observed_scores_effect': [], 'observed_scores_noneffect': [],
                             'observed_active_allele': ''}
    observed_seqlet_start = 0
    observed_seqlet_end = 0
    observed_seqlet_effect = 0
    observed_seqlet_noneffect = 0
    observed_bg_effect = 0
    observed_bg_noneffect = 0
    
    observed_effect_total = np.sum(np.array([max(0, np.sum(i)) for i in impscores_ref[seq][0:200]]))
    observed_noneffect_total = np.sum(np.array([max(0, np.sum(i)) for i in impscores_alt[seq][0:200]]))
    observed_per_seq_dict['observed_scores_effect'] = [np.sum(j) for j in impscores_ref[seq]]
    observed_per_seq_dict['observed_scores_noneffect'] = [np.sum(j) for j in impscores_alt[seq]]
    if observed_effect_total > observed_noneffect_total:
        observed_per_seq_dict['observed_active_allele'] = 'effect'
        observed_per_seq_dict['observed_inactive_allele'] = 'noneffect'
        observed_confidence_dict['observed_active_allele'].append('effect')
        observed_confidence_dict['observed_inactive_allele'].append('noneffect')
    else:
        observed_per_seq_dict['observed_active_allele'] = 'noneffect'
        observed_per_seq_dict['observed_inactive_allele'] = 'effect'
        observed_confidence_dict['observed_active_allele'].append('noneffect')
        observed_confidence_dict['observed_inactive_allele'].append('effect')
    observed_per_seq_dict['observed_scores_active'] = observed_per_seq_dict['observed_scores_'+observed_per_seq_dict['observed_active_allele']][0:min(len(observed_per_seq_dict['observed_scores_effect']),len(observed_per_seq_dict['observed_scores_noneffect']))]
    observed_per_seq_dict['observed_scores_inactive'] = observed_per_seq_dict['observed_scores_'+observed_per_seq_dict['observed_inactive_allele']][0:min(len(observed_per_seq_dict['observed_scores_effect']),len(observed_per_seq_dict['observed_scores_noneffect']))]
    
    start = 99
    end = 100
    while True:
        if observed_per_seq_dict['observed_scores_active'][start - 1] <= null_upper_thresh:
            if observed_per_seq_dict['observed_scores_active'][start - 2] <= null_upper_thresh:
                break
            else:
                start -= 1
        else:
            start -= 1
    while True:
        if observed_per_seq_dict['observed_scores_active'][end] <= null_upper_thresh:
            if observed_per_seq_dict['observed_scores_active'][end + 1] <= null_upper_thresh:
                break
            else:
                end += 1
        else:
            end += 1
    if start != (end - 1):
        observed_seqlet_len = end - start
        if (observed_seqlet_len) < 7:
            if (end - 100) > (99 - start):
                observed_seqlet_start = start - math.ceil((7 - observed_seqlet_len) / 2)
                observed_seqlet_end = end + math.floor((7 - observed_seqlet_len) / 2)
            else:           
                observed_seqlet_start = start - math.floor((7 - observed_seqlet_len) / 2)
                observed_seqlet_end = end + math.ceil((7 - observed_seqlet_len) / 2)
        else:
            observed_seqlet_start = start
            observed_seqlet_end = end
    else:
        observed_seqlet_start = start - 3
        observed_seqlet_end = end + 3
    
    observed_confidence_dict['observed_seqlet_start'].append(observed_seqlet_start)
    observed_confidence_dict['observed_seqlet_end'].append(observed_seqlet_end)
    observed_confidence_dict['observed_seqlet_len'].append(observed_seqlet_end - observed_seqlet_start)
    
    for i,j in enumerate(observed_per_seq_dict['observed_scores_active'][observed_seqlet_start:observed_seqlet_end]):
        if observed_per_seq_dict['observed_scores_effect'][i+observed_seqlet_start] >= 0:
            observed_seqlet_effect += observed_per_seq_dict['observed_scores_effect'][i+observed_seqlet_start]
        if observed_per_seq_dict['observed_scores_noneffect'][i+observed_seqlet_start] >= 0:
            observed_seqlet_noneffect += observed_per_seq_dict['observed_scores_noneffect'][i+observed_seqlet_start]
        assert j == observed_per_seq_dict['observed_scores_active'][i+observed_seqlet_start]
    for i,j in enumerate(observed_per_seq_dict['observed_scores_active']):
        if observed_per_seq_dict['observed_scores_effect'][i] >= 0:
            observed_bg_effect += observed_per_seq_dict['observed_scores_effect'][i]
        if observed_per_seq_dict['observed_scores_noneffect'][i] >= 0:
            observed_bg_noneffect += observed_per_seq_dict['observed_scores_noneffect'][i]
        assert j == observed_per_seq_dict['observed_scores_active'][i]
    
    observed_magnitude_score = observed_seqlet_effect - observed_seqlet_noneffect
    observed_prominence_score = (observed_seqlet_effect / observed_bg_effect) - (observed_seqlet_noneffect / observed_bg_noneffect)

    observed_confidence_dict['observed_seqlet_effect'].append(observed_seqlet_effect)
    observed_confidence_dict['observed_seqlet_noneffect'].append(observed_seqlet_noneffect)
    observed_confidence_dict['observed_bg_effect'].append(observed_bg_effect)
    observed_confidence_dict['observed_bg_noneffect'].append(observed_bg_noneffect)
    observed_confidence_dict['observed_magnitude_score'].append(observed_magnitude_score)
    observed_confidence_dict['observed_prominence_score'].append(observed_prominence_score)

print('Mean observed seqlet length:', '\t', '\t', statistics.mean(observed_confidence_dict['observed_seqlet_len']))
print('Median observed seqlet length:', '\t', '\t', statistics.median(observed_confidence_dict['observed_seqlet_len']))
print('St. Dev observed seqlet length:', '\t', statistics.stdev(observed_confidence_dict['observed_seqlet_len']))
print('Mode observed seqlet length:', '\t', '\t', statistics.mode(observed_confidence_dict['observed_seqlet_len']))
print('Max observed seqlet length:', '\t', '\t', max(observed_confidence_dict['observed_seqlet_len']))
print('Min observed seqlet lenght:', '\t', '\t', min(observed_confidence_dict['observed_seqlet_len']))

# Add observed scores to df
df2 = pd.DataFrame(observed_confidence_dict)
df2['name'] = df['name']
df2['alts'] = df['alts']
df = pd.merge(df, df2)
df.shape
df2.shape


# Code the same  Null back ground

null_confidence_dict = {'null_seqlet_start': [], 'null_seqlet_end': [],
                        'null_seqlet_effect': [], 'null_seqlet_noneffect': [],
                        'null_bg_effect': [], 'null_bg_noneffect': [],
                        'null_magnitude_score': [], 'null_prominence_score': [],
                        'null_active_allele': [], 'null_inactive_allele': [],
                        'null_seqlet_len': []}

score = []
for x in np.arange(len(shuf_impscores_ref)):
    # get tested element
    y = shuf_impscores_ref[x]
    # get score of the bps
    for z in np.arange(len(y)):
        score.append(y[z][np.nonzero(y[z])][0])


null_upper_thresh = np.quantile(score, 0.975)
print('Per-base Importance Score Threshold: ', '\t', null_upper_thresh)

merged_null_effect_scores = []
merged_null_noneffect_scores = []
for i in range(len(shuf_impscores_ref)):
    merged_null_effect_scores.append(shuf_impscores_ref[i].max(axis=1))
    merged_null_noneffect_scores.append(shuf_impscores_alt[i].max(axis=1))
    
for seq in range(len(shuf_impscores_ref)):
    null_per_seq_dict = {'null_scores_active': [], 'null_scores_inactive': [],
                         'null_scores_effect': [], 'null_scores_noneffect': [],
                         'null_active_allele': ''}
    null_seqlet_start = 0
    null_seqlet_end = 0
    null_seqlet_effect = 0
    null_seqlet_noneffect = 0
    null_bg_effect = 0
    null_bg_noneffect = 0
    
    null_effect_total = sum([max(0,i) for i in merged_null_effect_scores[seq][0:200]])
    null_noneffect_total = sum([max(0,i) for i in merged_null_noneffect_scores[seq][0:200]])
    null_per_seq_dict['null_scores_effect'] = merged_null_effect_scores[seq]
    null_per_seq_dict['null_scores_noneffect'] = merged_null_noneffect_scores[seq]
    if null_effect_total > null_noneffect_total:
        null_per_seq_dict['null_active_allele'] = 'effect'
        null_per_seq_dict['null_inactive_allele'] = 'noneffect'
        null_confidence_dict['null_active_allele'].append('effect')
        null_confidence_dict['null_inactive_allele'].append('noneffect')
    else:
        null_per_seq_dict['null_active_allele'] = 'noneffect'
        null_per_seq_dict['null_inactive_allele'] = 'effect'
        null_confidence_dict['null_active_allele'].append('noneffect')
        null_confidence_dict['null_inactive_allele'].append('effect')
    null_per_seq_dict['null_scores_active'] = null_per_seq_dict['null_scores_'+null_per_seq_dict['null_active_allele']][0:min(len(observed_per_seq_dict['observed_scores_effect']),len(observed_per_seq_dict['observed_scores_noneffect']))]
    null_per_seq_dict['null_scores_inactive'] = null_per_seq_dict['null_scores_'+null_per_seq_dict['null_inactive_allele']][0:min(len(observed_per_seq_dict['observed_scores_effect']),len(observed_per_seq_dict['observed_scores_noneffect']))]
    
    start = 99
    end = 100
    while True:
        if null_per_seq_dict['null_scores_active'][start - 1] <= null_upper_thresh:
            if null_per_seq_dict['null_scores_active'][start - 2] <= null_upper_thresh:
                break
            else:
                start -= 1
        else:
            start -= 1
    while True:
        if null_per_seq_dict['null_scores_active'][end] <= null_upper_thresh:
            if null_per_seq_dict['null_scores_active'][end + 1] <= null_upper_thresh:
                break
            else:
                end += 1
        else:
            end += 1
    if start != (end - 1):
        null_seqlet_len = end - start
        if (null_seqlet_len) < 7:
            if (end - 100) > (99 - start):
                null_seqlet_start = start - math.ceil((7 - null_seqlet_len) / 2)
                null_seqlet_end = end + math.floor((7 - null_seqlet_len) / 2)
            else:
                null_seqlet_start = start - math.floor((7 - null_seqlet_len) / 2)
                null_seqlet_end = end + math.ceil((7 - null_seqlet_len) / 2)
        else:
            null_seqlet_start = start
            null_seqlet_end = end
    else:
        null_seqlet_start = start - 3
        null_seqlet_end = end + 3
    
    null_confidence_dict['null_seqlet_start'].append(null_seqlet_start)
    null_confidence_dict['null_seqlet_end'].append(null_seqlet_end)
    null_confidence_dict['null_seqlet_len'].append(null_seqlet_end - null_seqlet_start)
    
    for i,j in enumerate(null_per_seq_dict['null_scores_active'][null_seqlet_start:null_seqlet_end]):
        if null_per_seq_dict['null_scores_effect'][i+null_seqlet_start] >= 0:
            null_seqlet_effect += null_per_seq_dict['null_scores_effect'][i+null_seqlet_start]
        if null_per_seq_dict['null_scores_noneffect'][i+null_seqlet_start] >= 0:
            null_seqlet_noneffect += null_per_seq_dict['null_scores_noneffect'][i+null_seqlet_start]
        assert j == null_per_seq_dict['null_scores_active'][i+null_seqlet_start]
    for i,j in enumerate(null_per_seq_dict['null_scores_active']):
        if null_per_seq_dict['null_scores_effect'][i] >= 0:
            null_bg_effect += null_per_seq_dict['null_scores_effect'][i]
        if null_per_seq_dict['null_scores_noneffect'][i] >= 0:
            null_bg_noneffect += null_per_seq_dict['null_scores_noneffect'][i]
        assert j == null_per_seq_dict['null_scores_active'][i]
    
    null_magnitude_score = null_seqlet_effect - null_seqlet_noneffect
    null_prominence_score = (null_seqlet_effect / null_bg_effect) - (null_seqlet_noneffect / null_bg_noneffect)
    
    null_confidence_dict['null_seqlet_effect'].append(null_seqlet_effect)
    null_confidence_dict['null_seqlet_noneffect'].append(null_seqlet_noneffect)
    null_confidence_dict['null_bg_effect'].append(null_bg_effect)
    null_confidence_dict['null_bg_noneffect'].append(null_bg_noneffect)
    null_confidence_dict['null_magnitude_score'].append(null_magnitude_score)
    null_confidence_dict['null_prominence_score'].append(null_prominence_score)

print('Mean Null seqlet length:', '\t', '\t', statistics.mean(null_confidence_dict['null_seqlet_len']))
print('Median Null seqlet length:', '\t', '\t', statistics.median(null_confidence_dict['null_seqlet_len']))
print('St. Dev Null seqlet length:', '\t', '\t', statistics.stdev(null_confidence_dict['null_seqlet_len']))
print('Mode Null seqlet length:', '\t', '\t', statistics.mode(null_confidence_dict['null_seqlet_len']))
print('Max Null seqlet length:', '\t', '\t', max(null_confidence_dict['null_seqlet_len']))
print('Min Null seqlet lenght:', '\t', '\t', min(null_confidence_dict['null_seqlet_len']))

# Plot Null prominence scores
distrib_name = 't'
distrib = getattr(scipy.stats, distrib_name)
prominence_params = distrib.fit(null_confidence_dict['null_prominence_score'])
x_prominence = np.linspace(distrib.ppf(0.0001, *prominence_params[:-2], prominence_params[-2], prominence_params[-1]),
                        distrib.ppf(0.9999, *prominence_params[:-2], prominence_params[-2], prominence_params[-1]), 10000)
y_prominence = distrib.pdf(x_prominence, *prominence_params[:-2], prominence_params[-2], prominence_params[-1])
distrib_prominence = pd.DataFrame(list(zip(x_prominence, y_prominence)), columns =['x', 'y'])

print('Fitted ' + distrib_name + ' Distribution: ')
print()
print('Mean:', '\t', prominence_params[-2])
print('Stdev:', '\t', prominence_params[-1])
print()
print(kstest(null_confidence_dict['null_prominence_score'], distrib_name, args=[*prominence_params[:-2], prominence_params[-2], prominence_params[-1]]))
print()
null_prominence_quantiles = np.quantile(null_confidence_dict['null_prominence_score'], [0.025, 0.975])


# Add P-value

df['prominence_pval'] = [(2 * min((percentileofscore(list(null_confidence_dict['null_prominence_score']), x) / 100),
                                          (1 - (percentileofscore(list(null_confidence_dict['null_prominence_score']), x) / 100)))) \
                                 for x in df['observed_prominence_score']]
df.shape

# Plot Magnitude scores

distrib_name = 't'
distrib = getattr(scipy.stats, distrib_name)
magnitude_params = distrib.fit(null_confidence_dict['null_magnitude_score'])
x_magnitude = np.linspace(distrib.ppf(0.0001, *magnitude_params[:-2], magnitude_params[-2], magnitude_params[-1]),
                        distrib.ppf(0.9999, *magnitude_params[:-2], magnitude_params[-2], magnitude_params[-1]), 10000)
y_magnitude = distrib.pdf(x_magnitude, *magnitude_params[:-2], magnitude_params[-2], magnitude_params[-1])
distrib_magnitude = pd.DataFrame(list(zip(x_magnitude, y_magnitude)), columns =['x', 'y'])

print('Fitted ' + distrib_name + ' Distribution: ')
print()
print('Mean:', '\t', magnitude_params[-2])
print('Stdev:', '\t', magnitude_params[-1])
print()
print(kstest(null_confidence_dict['null_magnitude_score'], distrib_name, args=[*magnitude_params[:-2], magnitude_params[-2], magnitude_params[-1]]))
print()

null_magnitude_quantiles = np.quantile(null_confidence_dict['null_magnitude_score'], [0.025, 0.975])



# Add P-value
df['magnitude_pval'] = [(2 * min((percentileofscore(list(null_confidence_dict['null_magnitude_score']), x) / 100),
                                          (1 - (percentileofscore(list(null_confidence_dict['null_magnitude_score']), x) / 100)))) \
                                 for x in df['observed_magnitude_score']]


#### Work on getting Negative scores

observed_confidence_dict = {'neg_observed_seqlet_start': [], 'neg_observed_seqlet_end': [],
                            'neg_observed_seqlet_effect': [], 'neg_observed_seqlet_noneffect': [],
                            'neg_observed_bg_effect': [], 'neg_observed_bg_noneffect': [],
                            'neg_observed_magnitude_score': [], 'neg_observed_prominence_score': [],
                            'neg_observed_active_allele': [], 'neg_observed_inactive_allele': [],
                            'neg_observed_seqlet_len': []}

score = []
for x in np.arange(len(shuf_impscores_ref)):
    # get tested element
    y = shuf_impscores_ref[x]
    # get score of the bps
    for z in np.arange(len(y)):
        score.append(y[z][np.nonzero(y[z])][0])

null_lower_thresh = np.quantile(score, 0.025)
print('Per-base Importance Score Threshold: ', '\t', null_lower_thresh)

for seq in range(len(impscores_ref)):
    observed_per_seq_dict = {'neg_observed_scores_active': [], 'neg_observed_scores_inactive': [],
                             'neg_observed_scores_effect': [], 'neg_observed_scores_noneffect': [],
                             'neg_observed_active_allele': ''}
    observed_seqlet_start = 0
    observed_seqlet_end = 0
    observed_seqlet_effect = 0
    observed_seqlet_noneffect = 0
    observed_bg_effect = 0
    observed_bg_noneffect = 0
    
    observed_effect_total = np.sum(np.array([max(0, np.sum(i)) for i in impscores_ref[seq][0:200]]))
    observed_noneffect_total = np.sum(np.array([max(0, np.sum(i)) for i in impscores_alt[seq][0:200]]))
    observed_per_seq_dict['neg_observed_scores_effect'] = [np.sum(j) for j in impscores_ref[seq]]
    observed_per_seq_dict['neg_observed_scores_noneffect'] = [np.sum(j) for j in impscores_alt[seq]]
    if observed_effect_total > observed_noneffect_total:
        observed_per_seq_dict['neg_observed_active_allele'] = 'effect'
        observed_per_seq_dict['neg_observed_inactive_allele'] = 'noneffect'
        observed_confidence_dict['neg_observed_active_allele'].append('effect')
        observed_confidence_dict['neg_observed_inactive_allele'].append('noneffect')
    else:
        observed_per_seq_dict['neg_observed_active_allele'] = 'noneffect'
        observed_per_seq_dict['neg_observed_inactive_allele'] = 'effect'
        observed_confidence_dict['neg_observed_active_allele'].append('noneffect')
        observed_confidence_dict['neg_observed_inactive_allele'].append('effect')
    observed_per_seq_dict['neg_observed_scores_active'] = observed_per_seq_dict['neg_observed_scores_'+observed_per_seq_dict['neg_observed_active_allele']][0:min(len(observed_per_seq_dict['neg_observed_scores_effect']),len(observed_per_seq_dict['neg_observed_scores_noneffect']))]
    observed_per_seq_dict['neg_observed_scores_inactive'] = observed_per_seq_dict['neg_observed_scores_'+observed_per_seq_dict['neg_observed_inactive_allele']][0:min(len(observed_per_seq_dict['neg_observed_scores_effect']),len(observed_per_seq_dict['neg_observed_scores_noneffect']))]
    
    start = 99
    end = 100
    while True:
        if observed_per_seq_dict['neg_observed_scores_active'][start - 1] >= null_lower_thresh:
            if observed_per_seq_dict['neg_observed_scores_active'][start - 2] >= null_lower_thresh:
                break
            else:
                start -= 1
        else:
            start -= 1
    while True:
        if observed_per_seq_dict['neg_observed_scores_active'][end] >= null_lower_thresh:
            if observed_per_seq_dict['neg_observed_scores_active'][end + 1] >= null_lower_thresh:
                break
            else:
                end += 1
        else:
            end += 1
    if start != (end - 1):
        observed_seqlet_len = end - start
        if (observed_seqlet_len) < 7:
            if (end - 100) > (99 - start):
                observed_seqlet_start = start - math.ceil((7 - observed_seqlet_len) / 2)
                observed_seqlet_end = end + math.floor((7 - observed_seqlet_len) / 2)
            else:           
                observed_seqlet_start = start - math.floor((7 - observed_seqlet_len) / 2)
                observed_seqlet_end = end + math.ceil((7 - observed_seqlet_len) / 2)
        else:
            observed_seqlet_start = start
            observed_seqlet_end = end
    else:
        observed_seqlet_start = start - 3
        observed_seqlet_end = end + 3
    
    observed_confidence_dict['neg_observed_seqlet_start'].append(observed_seqlet_start)
    observed_confidence_dict['neg_observed_seqlet_end'].append(observed_seqlet_end)
    observed_confidence_dict['neg_observed_seqlet_len'].append(observed_seqlet_end - observed_seqlet_start)
    
    for i,j in enumerate(observed_per_seq_dict['neg_observed_scores_active'][observed_seqlet_start:observed_seqlet_end]):
        if observed_per_seq_dict['neg_observed_scores_effect'][i+observed_seqlet_start] <= 0:
            observed_seqlet_effect += observed_per_seq_dict['neg_observed_scores_effect'][i+observed_seqlet_start]
        if observed_per_seq_dict['neg_observed_scores_noneffect'][i+observed_seqlet_start] <= 0:
            observed_seqlet_noneffect += observed_per_seq_dict['neg_observed_scores_noneffect'][i+observed_seqlet_start]
        assert j == observed_per_seq_dict['neg_observed_scores_active'][i+observed_seqlet_start]
    for i,j in enumerate(observed_per_seq_dict['neg_observed_scores_active']):
        if observed_per_seq_dict['neg_observed_scores_effect'][i] <= 0:
            observed_bg_effect += observed_per_seq_dict['neg_observed_scores_effect'][i]
        if observed_per_seq_dict['neg_observed_scores_noneffect'][i] <= 0:
            observed_bg_noneffect += observed_per_seq_dict['neg_observed_scores_noneffect'][i]
        assert j == observed_per_seq_dict['neg_observed_scores_active'][i]
    
        observed_magnitude_score = observed_seqlet_effect - observed_seqlet_noneffect
    
        if observed_bg_effect == 0:
            observed_prominence_score = 0
        elif observed_bg_noneffect == 0:
            observed_prominence_score = (observed_seqlet_effect / observed_bg_effect)
        else:
            observed_prominence_score = (observed_seqlet_effect / observed_bg_effect) - (observed_seqlet_noneffect / observed_bg_noneffect)
    
    
    observed_confidence_dict['neg_observed_seqlet_effect'].append(observed_seqlet_effect)
    observed_confidence_dict['neg_observed_seqlet_noneffect'].append(observed_seqlet_noneffect)
    observed_confidence_dict['neg_observed_bg_effect'].append(observed_bg_effect)
    observed_confidence_dict['neg_observed_bg_noneffect'].append(observed_bg_noneffect)
    observed_confidence_dict['neg_observed_magnitude_score'].append(observed_magnitude_score)
    observed_confidence_dict['neg_observed_prominence_score'].append(observed_prominence_score)

print('Mean observed seqlet length:', '\t', '\t', statistics.mean(observed_confidence_dict['neg_observed_seqlet_len']))
print('Median observed seqlet length:', '\t', '\t', statistics.median(observed_confidence_dict['neg_observed_seqlet_len']))
print('St. Dev observed seqlet length:', '\t', statistics.stdev(observed_confidence_dict['neg_observed_seqlet_len']))
print('Mode observed seqlet length:', '\t', '\t', statistics.mode(observed_confidence_dict['neg_observed_seqlet_len']))
print('Max observed seqlet length:', '\t', '\t', max(observed_confidence_dict['neg_observed_seqlet_len']))
print('Min observed seqlet lenght:', '\t', '\t', min(observed_confidence_dict['neg_observed_seqlet_len']))


# Add observed scores to df
df2 = pd.DataFrame(observed_confidence_dict)
df2['name'] = df['name']
df2['alts'] = df['alts']
df = pd.merge(df, df2)
df.shape
df2.shape

# Get Null Distrubiton for neg

# Code the same for the negative scores

null_confidence_dict = {'null_seqlet_start': [], 'null_seqlet_end': [],
                        'null_seqlet_effect': [], 'null_seqlet_noneffect': [],
                        'null_bg_effect': [], 'null_bg_noneffect': [],
                        'null_magnitude_score': [], 'null_prominence_score': [],
                        'null_active_allele': [], 'null_inactive_allele': [],
                        'null_seqlet_len': []}

score = []
for x in np.arange(len(shuf_impscores_ref)):
    # get tested element
    y = shuf_impscores_ref[x]
    # get score of the bps
    for z in np.arange(len(y)):
        score.append(y[z][np.nonzero(y[z])][0])


null_lower_thresh = np.quantile(score, 0.025)
print('Per-base Importance Score Threshold: ', '\t', null_lower_thresh)

    
for seq in range(len(shuf_impscores_ref)):
    null_per_seq_dict = {'null_scores_active': [], 'null_scores_inactive': [],
                         'null_scores_effect': [], 'null_scores_noneffect': [],
                         'null_active_allele': ''}
    null_seqlet_start = 0
    null_seqlet_end = 0
    null_seqlet_effect = 0
    null_seqlet_noneffect = 0
    null_bg_effect = 0
    null_bg_noneffect = 0
    
    
    null_effect_total = np.sum(np.array([max(0, np.sum(i)) for i in shuf_impscores_ref[seq][0:200]]))
    null_noneffect_total = np.sum(np.array([max(0, np.sum(i)) for i in shuf_impscores_alt[seq][0:200]]))
    null_per_seq_dict['null_scores_effect'] = [np.sum(j) for j in shuf_impscores_ref[seq]]
    null_per_seq_dict['null_scores_noneffect'] = [np.sum(j) for j in shuf_impscores_alt[seq]]
    if null_effect_total > null_noneffect_total:
        null_per_seq_dict['null_active_allele'] = 'effect'
        null_per_seq_dict['null_inactive_allele'] = 'noneffect'
        null_confidence_dict['null_active_allele'].append('effect')
        null_confidence_dict['null_inactive_allele'].append('noneffect')
    else:
        null_per_seq_dict['null_active_allele'] = 'noneffect'
        null_per_seq_dict['null_inactive_allele'] = 'effect'
        null_confidence_dict['null_active_allele'].append('noneffect')
        null_confidence_dict['null_inactive_allele'].append('effect')
    null_per_seq_dict['null_scores_active'] = null_per_seq_dict['null_scores_'+null_per_seq_dict['null_active_allele']][0:min(len(observed_per_seq_dict['neg_observed_scores_effect']),len(observed_per_seq_dict['neg_observed_scores_noneffect']))]
    null_per_seq_dict['null_scores_inactive'] = null_per_seq_dict['null_scores_'+null_per_seq_dict['null_inactive_allele']][0:min(len(observed_per_seq_dict['neg_observed_scores_effect']),len(observed_per_seq_dict['neg_observed_scores_noneffect']))]
    
    start = 99
    end = 100
    while True:
        if null_per_seq_dict['null_scores_active'][start - 1] >= null_lower_thresh:
            if null_per_seq_dict['null_scores_active'][start - 2] >= null_lower_thresh:
                break
            else:
                start -= 1
        else:
            start -= 1
    while True:
        if null_per_seq_dict['null_scores_active'][end] >= null_lower_thresh:
            if null_per_seq_dict['null_scores_active'][end + 1] >= null_lower_thresh:
                break
            else:
                end += 1
        else:
            end += 1
    if start != (end - 1):
        null_seqlet_len = end - start
        if (null_seqlet_len) < 7:
            if (end - 100) > (99 - start):
                null_seqlet_start = start - math.ceil((7 - null_seqlet_len) / 2)
                null_seqlet_end = end + math.floor((7 - null_seqlet_len) / 2)
            else:
                null_seqlet_start = start - math.floor((7 - null_seqlet_len) / 2)
                null_seqlet_end = end + math.ceil((7 - null_seqlet_len) / 2)
        else:
            null_seqlet_start = start
            null_seqlet_end = end
    else:
        null_seqlet_start = start - 3
        null_seqlet_end = end + 3
    
    null_confidence_dict['null_seqlet_start'].append(null_seqlet_start)
    null_confidence_dict['null_seqlet_end'].append(null_seqlet_end)
    null_confidence_dict['null_seqlet_len'].append(null_seqlet_end - null_seqlet_start)
    
    for i,j in enumerate(null_per_seq_dict['null_scores_active'][null_seqlet_start:null_seqlet_end]):
        if null_per_seq_dict['null_scores_effect'][i+null_seqlet_start] <= 0:
            null_seqlet_effect += null_per_seq_dict['null_scores_effect'][i+null_seqlet_start]
        if null_per_seq_dict['null_scores_noneffect'][i+null_seqlet_start] <= 0:
            null_seqlet_noneffect += null_per_seq_dict['null_scores_noneffect'][i+null_seqlet_start]
        assert j == null_per_seq_dict['null_scores_active'][i+null_seqlet_start]
    for i,j in enumerate(null_per_seq_dict['null_scores_active']):
        if null_per_seq_dict['null_scores_effect'][i] <= 0:
            null_bg_effect += null_per_seq_dict['null_scores_effect'][i]
        if null_per_seq_dict['null_scores_noneffect'][i] <= 0:
            null_bg_noneffect += null_per_seq_dict['null_scores_noneffect'][i]
        assert j == null_per_seq_dict['null_scores_active'][i]
    
    null_magnitude_score = null_seqlet_effect - null_seqlet_noneffect
    
    if null_bg_effect == 0:
        null_prominence_score = 0
    else:
        null_prominence_score = (null_seqlet_effect / null_bg_effect) - (null_seqlet_noneffect / null_bg_noneffect)
    
    null_confidence_dict['null_seqlet_effect'].append(null_seqlet_effect)
    null_confidence_dict['null_seqlet_noneffect'].append(null_seqlet_noneffect)
    null_confidence_dict['null_bg_effect'].append(null_bg_effect)
    null_confidence_dict['null_bg_noneffect'].append(null_bg_noneffect)
    null_confidence_dict['null_magnitude_score'].append(null_magnitude_score)
    null_confidence_dict['null_prominence_score'].append(null_prominence_score)

print('Mean Null seqlet length:', '\t', '\t', statistics.mean(null_confidence_dict['null_seqlet_len']))
print('Median Null seqlet length:', '\t', '\t', statistics.median(null_confidence_dict['null_seqlet_len']))
print('St. Dev Null seqlet length:', '\t', '\t', statistics.stdev(null_confidence_dict['null_seqlet_len']))
print('Mode Null seqlet length:', '\t', '\t', statistics.mode(null_confidence_dict['null_seqlet_len']))
print('Max Null seqlet length:', '\t', '\t', max(null_confidence_dict['null_seqlet_len']))
print('Min Null seqlet lenght:', '\t', '\t', min(null_confidence_dict['null_seqlet_len']))


# Plot Null prominence scores
distrib_name = 't'
distrib = getattr(scipy.stats, distrib_name)
prominence_params = distrib.fit(null_confidence_dict['null_prominence_score'])
x_prominence = np.linspace(distrib.ppf(0.0001, *prominence_params[:-2], prominence_params[-2], prominence_params[-1]),
                        distrib.ppf(0.9999, *prominence_params[:-2], prominence_params[-2], prominence_params[-1]), 10000)
y_prominence = distrib.pdf(x_prominence, *prominence_params[:-2], prominence_params[-2], prominence_params[-1])
distrib_prominence = pd.DataFrame(list(zip(x_prominence, y_prominence)), columns =['x', 'y'])

print('Fitted ' + distrib_name + ' Distribution: ')
print()
print('Mean:', '\t', prominence_params[-2])
print('Stdev:', '\t', prominence_params[-1])
print()
print(kstest(null_confidence_dict['null_prominence_score'], distrib_name, args=[*prominence_params[:-2], prominence_params[-2], prominence_params[-1]]))
print()
null_prominence_quantiles = np.quantile(null_confidence_dict['null_prominence_score'], [0.025, 0.975])


# Add P-value

df['neg_prominence_pval'] = [(2 * min((percentileofscore(list(null_confidence_dict['null_prominence_score']), x) / 100),
                                          (1 - (percentileofscore(list(null_confidence_dict['null_prominence_score']), x) / 100)))) \
                                 for x in df['neg_observed_prominence_score']]
df.shape

# Plot Magnitude scores

distrib_name = 't'
distrib = getattr(scipy.stats, distrib_name)
magnitude_params = distrib.fit(null_confidence_dict['null_magnitude_score'])
x_magnitude = np.linspace(distrib.ppf(0.0001, *magnitude_params[:-2], magnitude_params[-2], magnitude_params[-1]),
                        distrib.ppf(0.9999, *magnitude_params[:-2], magnitude_params[-2], magnitude_params[-1]), 10000)
y_magnitude = distrib.pdf(x_magnitude, *magnitude_params[:-2], magnitude_params[-2], magnitude_params[-1])
distrib_magnitude = pd.DataFrame(list(zip(x_magnitude, y_magnitude)), columns =['x', 'y'])

print('Fitted ' + distrib_name + ' Distribution: ')
print()
print('Mean:', '\t', magnitude_params[-2])
print('Stdev:', '\t', magnitude_params[-1])
print()
print(kstest(null_confidence_dict['null_magnitude_score'], distrib_name, args=[*magnitude_params[:-2], magnitude_params[-2], magnitude_params[-1]]))
print()

null_magnitude_quantiles = np.quantile(null_confidence_dict['null_magnitude_score'], [0.025, 0.975])



# Add P-value
df['neg_magnitude_pval'] = [(2 * min((percentileofscore(list(null_confidence_dict['null_magnitude_score']), x) / 100),
                                          (1 - (percentileofscore(list(null_confidence_dict['null_magnitude_score']), x) / 100)))) \
                                 for x in df['neg_observed_magnitude_score']]


#### Read in delta scores

DeltaSNPs = pd.read_table(delta, header=None, names=['loc','delta_score'])
print(DeltaSNPs.shape)
print(DeltaSNPs.head)

# Added delta score to df
df['delta_score'] = DeltaSNPs['delta_score']


DeltaShufSNPs = pd.read_table(delta_shuf, header=None, names=['loc','delta_score'])
print(DeltaShufSNPs.shape)
print(DeltaShufSNPs.head)

# Plot distrubtion
distrib_name = 't'
distrib = getattr(scipy.stats, distrib_name)
explain_params = distrib.fit(DeltaShufSNPs.delta_score)
x_explain = np.linspace(distrib.ppf(0.0001, *explain_params[:-2], explain_params[-2], explain_params[-1]),
                        distrib.ppf(0.9999, *explain_params[:-2], explain_params[-2], explain_params[-1]), 10000)
y_explain = distrib.pdf(x_explain, *explain_params[:-2], explain_params[-2], explain_params[-1])
distrib_explain = pd.DataFrame(list(zip(x_explain, y_explain)), columns =['x', 'y'])

print('Fitted ' + distrib_name + ' Distribution: ')
print()
print('Mean:', '\t', explain_params[-2])
print('Stdev:', '\t', explain_params[-1])
print()
print(kstest(DeltaShufSNPs.delta_score, distrib_name, args=[*explain_params[:-2], explain_params[-2], explain_params[-1]]))
print()


# Get pval
df['delta_p_val'] = [(2 * min(distrib.cdf(x, *explain_params[:-2], explain_params[-2], explain_params[-1]),
                                      1 - distrib.cdf(x, *explain_params[:-2], explain_params[-2], explain_params[-1]))) \
                              for x in df['delta_score']
                 ]


#### Get GKM Scores

snp_scores = []
for x in np.arange(len(impscores_ref)):
        # print(x)
        mid = 100
        # print(mid)
        low = mid - 25
        high = mid + 25
        
        
        score_ref = np.sum(impscores_ref[x][low:high])
        # print(score_ref)
        
        score_alt = np.sum(impscores_alt[x][low:high])
        # print(score_alt)
        
        score= score_alt - score_ref
        # print(score)
        
        snp_scores.append(score)
print(snp_scores[0:10])
plt.hist(snp_scores)

d = {'scores': snp_scores}
Explain_scores = pd.DataFrame(d)
Explain_scores

shuf_snp_scores = []
for x in np.arange(len(shuf_impscores_ref)):
        # print(x)
        mid = 100
        # print(mid)
        low = mid - 25
        high = mid + 25
        
        
        shuf_score_ref = np.sum(shuf_impscores_ref[x][low:high])
        # print(score_ref)
        
        shuf_score_alt = np.sum(shuf_impscores_alt[x][low:high])
        # print(score_alt)
        
        score = shuf_score_alt - shuf_score_ref
        # print(score)
        
        shuf_snp_scores.append(score)
print(shuf_snp_scores[0:10])
plt.hist(shuf_snp_scores)

d = {'shuffle_scores': shuf_snp_scores}
Explain_shuffle_scores = pd.DataFrame(d)
Explain_shuffle_scores

# Add scores to df
df['explain_score'] = snp_scores

# Plot distrubtiution
distrib_name = 't'
distrib = getattr(scipy.stats, distrib_name)
explain_params = distrib.fit(shuf_snp_scores)
x_explain = np.linspace(distrib.ppf(0.0001, *explain_params[:-2], explain_params[-2], explain_params[-1]),
                        distrib.ppf(0.9999, *explain_params[:-2], explain_params[-2], explain_params[-1]), 10000)
y_explain = distrib.pdf(x_explain, *explain_params[:-2], explain_params[-2], explain_params[-1])
distrib_explain = pd.DataFrame(list(zip(x_explain, y_explain)), columns =['x', 'y'])

print('Fitted ' + distrib_name + ' Distribution: ')
print()
print('Mean:', '\t', explain_params[-2])
print('Stdev:', '\t', explain_params[-1])
print()
print(kstest(shuf_snp_scores, distrib_name, args=[*explain_params[:-2], explain_params[-2], explain_params[-1]]))
print()


# Get pval
df['explain_pval'] = [(2 * min(distrib.cdf(x, *explain_params[:-2], explain_params[-2], explain_params[-1]),
                                      1 - distrib.cdf(x, *explain_params[:-2], explain_params[-2], explain_params[-1]))) \
                              for x in df['explain_score']]

#### Get ISM value
ISM_ref_df = pd.read_table(ISM_ref, header = None)
ISM_alt_df = pd.read_table(ISM_alt, header = None)
ISM_ref_shuf_df = pd.read_table(ISM_ref_shuf, header = None)
ISM_alt_shuf_df = pd.read_table(ISM_alt_shuf, header = None)


ISM_score = ISM_alt_df[1] - ISM_ref_df[1]
ISM_shuf_score = ISM_alt_shuf_df[1] - ISM_ref_shuf_df[1]

df['ISM_score'] = ISM_score

distrib_name = 't'
distrib = getattr(scipy.stats, distrib_name)
explain_params = distrib.fit(ISM_shuf_score)
x_explain = np.linspace(distrib.ppf(0.0001, *explain_params[:-2], explain_params[-2], explain_params[-1]),
                        distrib.ppf(0.9999, *explain_params[:-2], explain_params[-2], explain_params[-1]), 10000)
y_explain = distrib.pdf(x_explain, *explain_params[:-2], explain_params[-2], explain_params[-1])
distrib_explain = pd.DataFrame(list(zip(x_explain, y_explain)), columns =['x', 'y'])

print('Fitted ' + distrib_name + ' Distribution: ')
print()
print('Mean:', '\t', explain_params[-2])
print('Stdev:', '\t', explain_params[-1])
print()
print(kstest(ISM_shuf_score, distrib_name, args=[*explain_params[:-2], explain_params[-2], explain_params[-1]]))
print()

# Get pval
df['ISM_pval'] = [(2 * min(distrib.cdf(x, *explain_params[:-2], explain_params[-2], explain_params[-1]),
                                      1 - distrib.cdf(x, *explain_params[:-2], explain_params[-2], explain_params[-1]))) \
                              for x in df['ISM_score']]

#### Write the output
df.to_csv(output)