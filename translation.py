
import random as rnd

# Basic DNA and RNA bp parameters
DNA_bases = {'T', 'C', 'A', 'G'}
RNA_bases = {s for s in DNA_bases if s!='T'}.union('U')

def DNA2RNA(seq):
    return seq.upper().replace('T', 'U')

def RNA2DNA(seq):
    return seq.upper().replace('U', 'T')

def is_valid_DNA_seq(seq):
    '''Returns True [False] is seq is a valid [not a valid] DNA sequence'''
    return set(seq.upper()) <= DNA_bases

def is_valid_RNA_seq(seq):
    '''Returns True [False] is seq is a valid [not a valid] RNA sequence'''
    return set(seq.upper()) <= RNA_bases

def is_valid_bp_seq(seq, seq_type='DNA'):
    '''Returns True [False] is seq is a valid [not a valid] seq_type sequence'''
    return is_valid_DNA_seq(seq) if seq_type=='DNA' else is_valid_RNA_seq(seq)

# Codon translation table
RNA_codon_table = {
    #    U             C             A             G
# U
    'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys', # UxU
    'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys', # UxC
    'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---', # UxA
    'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Trp', # UxG
# C
    'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg', # CxU
    'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg', # CxC
    'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg', # CxA
    'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg', # CxG
# A
    'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser', # AxU
    'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser', # AxC
    'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg', # AxA
    'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg', # AxG
# G
    'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly', # GxU
    'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly', # GxC
    'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly', # GxA
    'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'  # GxG
}

# reverse translation ({aa : a codon of aa})
RNA_aa_table = {v : k for k, v in RNA_codon_table.items()}

# reverse translation with synonymous codons {aa : [all codons of aa]}
RNA_aa_synon_table = {v: [k for k in RNA_codon_table.keys() if RNA_codon_table[k]==v] \
                      for v in set(RNA_codon_table.values())}

# AA code (Three-letter code : [IUPAC, Name])
AA_code = {
    'Ala' : ['A', 'Alanine'],
    'Arg' : ['R', 'Arginine'], 
    'Asn' : ['N', 'Asparagine'],
    'Asp' : ['D', 'Aspartic-acid'],
    'Cys' : ['C', 'Cysteine'],
    'Glu' : ['E', 'Glutamic-acid'],
    'Gln' : ['Q', 'Glutamine'],
    'Gly' : ['G', 'Glycine'],
    'His' : ['H', 'Histidine'],
    'Ile' : ['I', 'Isoleucine'],
    'Leu' : ['L', 'Leucine'],
    'Lys' : ['K', 'Lysine'],
    'Met' : ['M', 'Methionine'],
    'Phe' : ['F', 'Phenylalanine'],
    'Pro' : ['P', 'Proline'],
    'Ser' : ['S', 'Serine'],
    'Thr' : ['T', 'Threonine'],
    'Trp' : ['W', 'Tryptophan'],
    'Tyr' : ['Y', 'Tyrosine'],
    'Val' : ['V', 'Valine']
}

# aa IUPAC to three-letter code (IUPAC AA: three letters AA)
AA_IUPAC23code = {v[0] : k for k, v in AA_code.items()}

def convert_aa_IUPAC23(aa_seq):
    '''Converts a IUPAC AA sequence to three-letter code AA sequence'''
    return ''.join([AA_IUPAC23code.get(i, '???') for i in list(aa_seq)])

def convert_aa_32IUPAC(aa_seq):
    '''Converts a three-letter AA code sequence to IUPAC AA sequence'''
    return ''.join([AA_code.get(aa_seq[i:i+3], '?')[0] for i in range(0,len(aa_seq),3)])

def translate_RNA_codon(codon):
    '''Returns the animo acid corresponding to a codon'''
    return RNA_codon_table[DNA2RNA(codon)]

def translate_RNA_aa(aa):
    '''Returns a codon corresponding to an animo acid'''
    return(RNA_aa_table.get(aa, '???'))

def translate(seq):
    '''Returns the animo acid sequence corresponding to the RNA sequence seq'''
    translation = ''
    for n in range(0, len(seq) - (len(seq) % 3), 3): # every third base
        translation += translate_RNA_codon(seq[n:n+3]) 
    return translation

# this can be used for both translation and reverse translation
def gen_translate(seq, func=translate_RNA_codon):
    '''Performs translation based on a translation function'''
    translation = ''
    for n in range(0, len(seq) - (len(seq) % 3), 3): # every third base of three letters of AA
        translation += func(seq[n:n+3]) 
    return translation

def translate_in_frame(seq, framenum=1):
    '''Returns the translation of seq in a given reading frame'''
    #return translate(seq[framenum-1:])
    return gen_translate(seq[framenum-1:])

def print_translation_in_frame(seq, framenum=1, prefix=''):
    '''Prints the translation of seq in reading frame framenum preceded by prefix'''
    print(prefix, framenum, ' ' * framenum, translate_in_frame(seq, framenum), sep='')
    
def print_translations(seq, prefix=''):
    '''Prints the translation of seq in all three reading frames, each preceded by prefix'''
    print('\n' ,' ' * (len(prefix) + 2), seq, sep='', end=':\n')
    [print_translation_in_frame(seq, fnum, prefix) for fnum in range(1,4)] 
        

# Translation of ORF (start codon to stop codon)
def translate_with_open_reading_frames(seq, framenum=1, open=False):
    '''Returns the translation of seq in framenum (1, 2, or 3), with ---'s when not 
    within an open reading frame'''
    translation = ''
    seqlength = len(seq) - (framenum - 1)
    for n in range(framenum-1, seqlength - (seqlength % 3), 3):
        codon = translate_RNA_codon(seq[n:n+3])
        open = (open or codon == "Met") and not (codon == "---")
        translation += codon if open else "---" 
    return translation

def print_translation_with_open_reading_frame(seq, framenum=1, prefix=''): 
    '''Prints ORF translation of a seq'''
    print(prefix, framenum, ' ' * framenum, translate_with_open_reading_frames(seq, framenum), sep='')

def print_translations_with_open_reading_frames(seq, prefix=''):
    '''Prints ORF translation of a seq for the three reading frame cases'''
    print('\n', ' ' * (len(prefix) + 2), seq, sep='')
    [print_translation_with_open_reading_frame(seq, frame, prefix) for frame in range(1,4)]
    
    
def gen_random_bp_seq(len, DNA=True):
    '''Generates a random sequence of DNA (default) or RNA nucleotides'''
    return ''.join(rnd.choices(list(DNA_bases),k=len)) if DNA else ''.join(rnd.choices(list(RNA_bases),k=len))

def gen_random_aa(len, IUPAC=True):
    '''Generates a random sequence of AA IUPAC letters (default) or a sequence of three letter-code AA'''
    return ''.join(rnd.choices(list(AA_IUPAC23code.keys()),k=len)) if IUPAC \
        else ''.join(rnd.choices(list(AA_code),k=len))
# ===============================================================================================================
