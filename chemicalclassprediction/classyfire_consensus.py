#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 16:19:11 2018

@author: madeleineernst
"""
import pandas as pd
from functools import reduce

import functools
import operator
from collections import Counter

def unique_smiles(matches):
    
    # combine SMILES for same feature into one string features
    for index, item in enumerate(matches):
        if 'Scan' in matches[index].columns:
            matches[index] = matches[index].groupby('Scan', as_index=False).agg(lambda x: ','.join(set(x.dropna()))) 
            matches[index] = matches[index].rename(columns = {'Scan':'cluster.index'})
            
    comb = reduce(lambda left,right: pd.merge(left,right,on='cluster.index', how = "outer"), matches)
    if 'FusionSMILES' in comb.columns:
        comb = comb.drop(['FusionSMILES', 'ConsensusSMILES'], axis=1)
    
    # concatenate all SMILES
    comb['AllSmiles'] = comb.filter(regex='^.*(Smiles|SMILES).*$').apply(lambda x: ','.join([y for y in x.tolist() if str(y) != 'nan']), axis=1)
    
    # create dictionary of SMILES
    comb_dic = comb.set_index('cluster.index')['AllSmiles'].to_dict()
    
    # split comma separated strings of SMILES 
    for k in comb_dic:
        comb_dic[k] = comb_dic[k].split(',')
    
    # retrieve unique SMILES per feature
    for i in comb_dic:
        comb_dic[i] = list(set(comb_dic[i]))
        
    # remove empty values
    comb_dic = {k: comb_dic[k] for k in comb_dic if not comb_dic[k] == ['']} # 2172
    
    # convert dictionary into list of unique SMILES
    l = list(set([item for sublist in list(comb_dic.values()) for item in sublist]))
    
    # convert list into dataframe
    df = pd.DataFrame({"SMILES": l})
    
    return {'df':df, 'dic':comb_dic}

# retrieve most predominant chemical group per componentindex at one hierarchical level
def highestscore(a, chem_dic, score):
    
    # creates list of lists, where each list item corresponds to one componentindex
    # each sublist item corresponds to the chemical classes of the SMILES matched to one node
    chem_ci = []
    for i in a:
        chem_ci.append([chem_dic[x] for x in i])

    # creates a list of dictionaries. Each list item corresponds to one componentindex and each
    # dictionary corresponds to the chemical categories found within this node with a corresponding score.
    # This score sums up to 1 per node and corresponds to the number of SMILES associated with this class 
    # divided by the total number of SMILES per node
    chem_scores = [] 
    for i in chem_ci:
        subl = []
        for j in i:
            dic_a = dict(Counter(j))
            out = {k: v / total for total in (sum(dic_a.values()),) for k, v in dic_a.items()}
            subl.append(out)
        chem_scores.append(subl)

    # List of dictionaries, each list item corresponds to one componentindex and dictionaries contain
    # summed values per chemical category per component index
    # E.g. Organic acids and derivatives': 2.1246753246753247 
    # means that 2.12 nodes in this molecular family were classified as Organic acids and derivatives
    chem_sums = []
    for i in chem_scores:
        if i:
            sums = functools.reduce(operator.add, map(Counter, i))
        if not i:
            sums = []
        chem_sums.append(sums)

    # list of items, each list item corresponds to one componentindex and name and score of the most highly abundant
    # chemical class are found in it 
    # E.g. Lipids and 0.75 means that 0.75 of the molecular family was associated to lipids, 
    # this could be the case if e.g. Lipids score is 2.25 (nodes) and the molecular family consits of a total of 3 nodes
    chem_finalscore = []
    for index, item in enumerate(chem_sums): 
        if chem_sums[index]:
            char = max(chem_sums[index].items(), key=operator.itemgetter(1))
            num = char[1]/score[index]
        if not chem_sums[index]:
            char = ["no matches","no matches"]
            num = "no matches"
        chem_finalscore.append([char[0],num])
    
    return chem_finalscore


def molfam_classes(net, df, smilesdict):
    
    # df is a pandas dataframe comprising all unique SMILES, inchikeys and corresponding ClassyFire ontology
    # net is the GNPS network data
    # smilesdict is a dictionary of nodes with corresponding unique SMILES
    
    # create dictionaries for each hierarchical level of Classyfire

    kingdom_dic = smilesdict.copy()
    d = {k: v for k, v in zip(df.SMILES, df.kingdom)}
    for k in kingdom_dic:
            kingdom_dic[k] = [d.get(k, k) for k in kingdom_dic[k]]

    # Superclass
    superclass_dic = smilesdict.copy()
    d = {k: v for k, v in zip(df.SMILES, df.superclass)}
    for k in superclass_dic:
            superclass_dic[k] = [d.get(k, k) for k in superclass_dic[k]]
            
    # Class
    class_dic = smilesdict.copy()
    d = {k: v for k, v in zip(df.SMILES, df.CF_class)}
    for k in class_dic:
            class_dic[k] = [d.get(k, k) for k in class_dic[k]]
    
    # Sublass
    subclass_dic = smilesdict.copy()
    d = {k: v for k, v in zip(df.SMILES, df.subclass)}
    for k in subclass_dic:
            subclass_dic[k] = [d.get(k, k) for k in subclass_dic[k]]
            
    # Direct parent
    Dparent_dic = smilesdict.copy()
    d = {k: v for k, v in zip(df.SMILES, df.direct_parent)}
    for k in Dparent_dic:
            Dparent_dic[k] = [d.get(k, k) for k in Dparent_dic[k]]
            
    # Molecular framework
    MFramework_dic = smilesdict.copy()
    d = {k: v for k, v in zip(df.SMILES, df.molecular_framework)}
    for k in MFramework_dic:
            MFramework_dic[k] = [d.get(k, k) for k in MFramework_dic[k]]
    
    # iterate over all cluster indexes
    ci = list(set(net['componentindex']))
    
    # get list of all cluster indexes per componentindex
    a = []
    for i in ci:
        a.append(list(net.loc[net['componentindex'] == i, 'cluster index']))
        
    # remove all items that do not have any SMILES (are not contained in keys of smilesdict)
    # score retrieves number of nodes per cluster index
    score = []
    for index, item in enumerate(a): 
        red = list(set(list(smilesdict.keys())).intersection(a[index]))
        score.append(len(a[index]))
        a[index] = red 
    
    kingdom_finalscore = highestscore(a, kingdom_dic, score)
    superclass_finalscore = highestscore(a, superclass_dic, score)
    class_finalscore = highestscore(a, class_dic, score)
    subclass_finalscore = highestscore(a, subclass_dic, score)
    Dparent_finalscore = highestscore(a, Dparent_dic, score)
    MFramework_finalscore = highestscore(a, MFramework_dic, score)
   
    kingdom = [item[0] for item in kingdom_finalscore]
    superclass = [item[0] for item in superclass_finalscore]
    CF_class = [item[0] for item in class_finalscore]
    subclass = [item[0] for item in subclass_finalscore]
    Dparent = [item[0] for item in Dparent_finalscore]
    MFramework = [item[0] for item in MFramework_finalscore]
    kingdom_score = [item[1] for item in kingdom_finalscore]
    superclass_score =[item[1] for item in superclass_finalscore]
    CF_class_score =[item[1] for item in class_finalscore]
    subclass_score =[item[1] for item in subclass_finalscore]
    Dparent_score =[item[1] for item in Dparent_finalscore]
    MFramework_score =[item[1] for item in MFramework_finalscore]
    
        
    data_tuples = list(zip(ci,score,kingdom,kingdom_score,superclass,superclass_score,CF_class,CF_class_score,subclass,subclass_score,Dparent,Dparent_score,MFramework,MFramework_score))
    sumary = pd.DataFrame(data_tuples, columns=['componentindex','CF_NrNodes','CF_kingdom','CF_kingdom_score','CF_superclass','CF_superclass_score','CF_class','CF_class_score','CF_subclass','CF_subclass_score','CF_Dparent','CF_Dparent_score','CF_MFramework','CF_MFramework_score'])
        
    # remove selfloops
    sumary = sumary[sumary.componentindex != -1]
        
    final = pd.merge(net[['componentindex','cluster index']], sumary, on='componentindex')
    # make cluster index first column
    final = final[['cluster index','componentindex','CF_NrNodes','CF_kingdom','CF_kingdom_score','CF_superclass','CF_superclass_score','CF_class','CF_class_score','CF_subclass','CF_subclass_score','CF_Dparent','CF_Dparent_score','CF_MFramework','CF_MFramework_score']]
    
    return final
        
    
    