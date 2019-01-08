import itertools as itr
import csv
import re

def NT_base_seq_score(two_seqs, eql=+5, neql=-4):
    '''Returns the score of the distance between two sequences. Default scorings are nuc44-based.'''
    return(sum([eql if s1==s2 else neql for s1,s2 in zip(two_seqs[0], two_seqs[1])]))


def NT_pair_seqs_score(all_seqs, score_func, eql=+5, neql=-4):
    '''Returns a dictionary where the keys are all possible unique pairs of all_seqs, and the values are
    the corresponding scores based on score_func'''
    scores = {}
    for a in itr.combinations(range(len(all_seqs)), 2):
        #scores[all_seqs[a[0]]+'~'+all_seqs[a[1]]] = score_func((all_seqs[a[0]], all_seqs[a[1]]), eql, neql)
        scores[(all_seqs[a[0]], all_seqs[a[1]])] = score_func((all_seqs[a[0]], all_seqs[a[1]]), eql, neql)
    return scores


def parse_starcode_cls_out(file):
    '''This function parses the output of the tool starcode that 
    is used to cluster nucleotide sequences. The returned dictionary contains
    a centroid as a key and a list of the correponding sequences in that
    cluster as a value.'''
    with open(file, 'rt') as fin:
        clusts = [r[0].split('\t') for r in csv.reader(fin, delimiter=' ')]

    clusters = {}
    for info in clusts: clusters[info[0]] = info[-1].split(',')
    return clusters

def parse_meshclsut_cls_out(file):
    '''This function parses the output of the tool Meshclust that 
    is used to cluster nucleotide sequences. The returned dictionary contains
    a centroid as a key and a list of the correponding sequences in that
    cluster as a value.'''
    clusters = {}
    with open(file, 'rt') as fin:
        for line in fin:
            vals = line.split(' ')
            if vals[0]=='>Cluster': # a new cluster
                clst_num = int(vals[-1][0])
                if clst_num > 0: # save the previous cluster
                    clusters[centroid]=clst_seqs
                clst_seqs = []
                centroid = None
            else:
                seq = re.search('(\w+)', vals[2]).group(0)
                clst_seqs.append(seq)
                if(vals[-1][0]=='*'):
                    centroid=seq
        # save the last cluster
        clusters[centroid]=clst_seqs
    return clusters
