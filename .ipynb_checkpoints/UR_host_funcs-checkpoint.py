
import pandas as pd
import re
from functools import reduce
from collections import Counter
import pald_funcs as pldfunc

def intersect_series_values(df, func):
    '''Returns the intersection of lists in a Series input. In the Series
    input, each row is a list. The function returns the intersection over 
    all lists. Currently not used'''
    all_sets = [set(df.iloc[i][:]) for i in range(df.shape[0])] # convert to a list of sets
    return reduce((lambda x, y: x & y), all_sets)
    
def reduce_series_values(df, func):
    '''Reduces the Series values using 'func'. In the Series
    input, each row is a list. The function returns a reduce operation over 
    all lists, where the operation is defined by the input "func"'''
    # convert to a list of sets (each set in a Series value (list)) and then reduce
    return reduce(func, [set(df.iloc[i][:]) for i in range(df.shape[0])])

def most_series_values(df, perc=90.0):
    '''df is a Series, where each value is a list. The function returns the elements that appear 
    at least perc% in all lists. Here, an element can appear only once in each list.'''
    thres = df.count()*perc/100.0
    # aggregate all URs (with multiplicity), count the multiplicity (Counter) and return only
    # elements with multiplicity >= thres
    return {k for k, v in Counter(df.aggregate('sum')).items() if v>=thres}

def rename_hdb(name):
    '''This provides a better host DB names'''
    return 'Vertebrate' if name=='Ensembl' else name.split('_')[-1]
    
    
def cluster_NTs(ur, mlen):
    '''Clusters list of nucleotide (each of length mlen) to several groups:
    1. 'eql': All with same NT,
    2. 'gc': All with just G and C,
    3. 'at': All with just A and T, 
    4. 'noneql': ur \ {'eql'}, 
    5. 'nongc': ur \ {'gc'},
    6. 'nonat': ur \ {'at'}.
    7. 'pald': Palindromic sequences (i.e. x==reverse(complement(x)))
    Input 'ur' must be a set.'''
    clst_nt = {}
    clst_nt['all'] = ur
    clst_nt['eql'] = {s for s in ur if len(set(s))==1}  # URs with the same nucleotides (e.g., AAAA)
    clst_nt['gc'] = {s for s in ur if re.search('[CG]{'+str(mlen)+'}',s.upper())}  # URs that contain only G and/or C
    clst_nt['at'] = {s for s in ur if re.search('[AT]{'+str(mlen)+'}',s.upper())}  # URs that contain only A and/or T
    clst_nt['neql'] = ur - clst_nt['eql']  # all except URs with equal NTs
    clst_nt['ngc'] = ur - clst_nt['gc']    # all except URs with only G and/or C
    clst_nt['nat'] = ur - clst_nt['at']    # all except URs with only A and/or T
    clst_nt['pald'] = {s for s, f in zip(ur, pldfunc.is_palindromic(list(ur))) if f}  # Palindromic URs (e.g. GATC)
    return clst_nt

nhdb = {  # simpler names for host DB 
'Ensembl' : 'Vertebrate',
'Ensembl_Bacteria' : 'Bacteria',
'Ensembl_Fungi' : 'Fungi',
'Ensembl_Metazoa' : 'Metazoa',
'Ensembl_Plants' : 'Plants', 
'Ensembl_Protists' : 'Protists'
}
