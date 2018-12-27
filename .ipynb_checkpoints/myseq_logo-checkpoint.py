from Bio import motifs, SeqIO
from Bio.Seq import Seq
import numpy as np
import math
import matplotlib as mpl
from matplotlib.text import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
import itertools as it

def compute_pssm(lseq, pscount={'A':0.25,'C':0.25,'G':0.25,'T':0.25}, \
                 background={'A':0.25,'C':0.25,'G':0.25,'T':0.25}):
    '''Computes the PFM, PWM, and PSSM of a list of nucleotide sequences'''
    m = motifs.create(list(map(Seq, lseq)))

    pfm = m.counts
    # pwm with no pscount is a stochastic matrix
    pwm = m.counts.normalize(pseudocounts=pscount)
    pssm = pwm.log_odds(background)
    return pfm, pwm, pssm, m, m.consensus


# The following code is based on:
# https://github.com/saketkc/motif-logos-matplotlib/blob/master/Sequence%20logos%20in%20Python.ipynb
# However, the code there has few errors. The code below follows the paper "Information Content
# of Binding Sites of Nucleotide Sequences"
def calc_IC_approx_err(motif):
    '''Approximate calculate of small-sample correction error'''
    print('Computing approximate correction error...')
    bases = list(motif.pwm.keys())
    n = len(motif.counts[bases[0]])  # sequence length 
    return (len(bases)-1)/(2 * n * np.log(2))

def calc_IC_exact_err(motif):
    '''Exact computation of small-sample correction error'''
    print('Computing exact correction error...')
    pwm = motif.pwm
    bases = list(pwm.keys())
    n = na = len(motif.counts['A'])  # sequence length 
    nc = ng = nt = exact_error = 0
    done = False
    while not done:
        #print (na,nc,ng,nt)
        pp = (0.25**na)*(0.25**nc)*(0.25**ng)*(0.25**nt)
        frac = pp*math.factorial(na+nc+ng+nt)/(math.factorial(na)*math.factorial(nc)*\
                                               math.factorial(ng)*math.factorial(nt))
        exact_error += frac*sum([-p*np.nan_to_num(np.log2(p)) for p in \
                                 [na/n, nc/n, ng/n, nt/n]])
        if nt<=0:
            ## iterate inner loop            
            if ng > 0:
                ## g => t
                ng = ng - 1
                nt = nt + 1
            elif nc > 0:
                ## c -> g 
                nc = nc - 1;
                ng = ng + 1;
            else:
                ## a->c
                na = na - 1
                nc = nc + 1
        else:
            if ng > 0:
                ## g => t
                ng = ng - 1 
                nt = nt + 1
            elif nc>0:
                ## c => g; all t -> g
                nc = nc - 1
                ng = nt + 1
                nt = 0
            elif na>0:
                ## a => c; all g,t -> c
                nc = nt + 1
                na = na - 1
                nt = 0
            else:
                done = True
    return exact_error

def calc_info_content(motif, corr_type = 'no'):
    '''Calculate information content with small sample correction.
    Note that for both corr_type=='approx' (should be used for
    sequences larger than 50 nt) and for corr_type=='exact' (should
    be used for sequences smaller than 50 nt), the output can attain
    both negative and positive values. See the paper "Information Content
    of Binding Sites of Nucleotide Sequences". Thus, for sequence logos, use
    the default (i.e. corr_type = 'no). Note, the PWM used is motif.pwm (i.e. PWM
    without pseudocounts)
    '''
    pwm = motif.pwm  # should not use relative information
    bases = list(pwm.keys())
    if corr_type=='no':
        Hg = np.log2(len(bases))
    elif corr_type=='approx':
        Hg = np.log2(len(bases)) - calc_IC_approx_err(motif)
    else:  # exact 
        Hg = calc_IC_exact_err(motif)
    #print('Hg = {}'.format(Hg))
    return [Hg+sum([pwm[b][l]*np.nan_to_num(np.log2(pwm[b][l])) for b in bases])\
            for l in range(0, len(motif))]

def calc_rel_info(motif, corr_type = 'no'):
    '''Calculate relative information, i.e. the information content
    of each base along the sequence. Note, the PWM used is motif.pwm (i.e. PWM
    without pseudocounts)'''
    info_cont = calc_info_content(motif, corr_type)
    return {b: [np.nan_to_num(p*i) for p, i in zip(motif.pwm[b], info_cont)] for b in list(motif.pwm.keys())}

def gen_nt_sequence_logo(ax, rel_info):
    '''Generates nucleotide sequence logo.
    rel_info is computed by: rel_info=calc_rel_info(motif, 'no')'''
    fp = FontProperties(family="Arial", weight="bold") 
    globscale = 1.35  # this, and the values here below were set for family="Arial"
    LETTERS = { 'T' : TextPath((-0.305, 0), "T", size=1, prop=fp),
                'G' : TextPath((-0.384, 0), "G", size=1, prop=fp),
                'A' : TextPath((-0.35, 0), "A", size=1, prop=fp),
                'C' : TextPath((-0.366, 0), "C", size=1, prop=fp) }
    COLOR_SCHEME = {'G': 'orange', 
                    'A': 'red', 
                    'C': 'blue', 
                    'T': 'darkgreen'}

    def letterAt(letter, x, y, yscale=1, ax=None):
        '''Plots a letter at a given position with a given scale.'''
        text = LETTERS[letter]

        t = mpl.transforms.Affine2D().scale(1*globscale, yscale*globscale) + \
            mpl.transforms.Affine2D().translate(x,y) + ax.transData
        p = PathPatch(text, lw=0, fc=COLOR_SCHEME[letter], transform=t)
        if ax != None: ax.add_artist(p)
        return p

    cnt = it.count(1)
    maxy = 0
    for i in range(0, len(list(rel_info.values())[0])):
        scores = [(b,rel_info[b][i]) for b in COLOR_SCHEME.keys()]
        scores.sort(key=lambda t: t[1])
        y, x = 0, next(cnt)
        for base, score in scores:
            letterAt(base, x, y, score, ax)
            y += score
        maxy = max(maxy, y)
    return maxy # this can be used to set the y-axis limit in the figure
