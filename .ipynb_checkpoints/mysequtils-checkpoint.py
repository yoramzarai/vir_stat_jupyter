from Bio import SeqIO

def seq_create_fasta(seq, labels=False):
    '''Given a sequence of strings and correponding labels, return a string containing 
    the FASTA data in the form of:
    > labels[0]
    seq[0]
    > labels[1]
    seq[1]
    ...
    '''
    if not labels: labels=['s'+str(s) for s in range(len(seq))]
    assert len(seq)==len(labels), "seq and labels length missmatch !"
    return ''.join([s+'\n' if f else '> '+l+'\n' for s,l in zip(seq, labels) for f in range(2)])

    
def seq2plain_text_file(seq, fname):
    '''This function writes the input list of sequences to a file, where
    each sequence is written in a separate line.'''
    with open(fname, 'wt') as fout: print(*seq, sep='\n', file=fout)

def create_fasta_file(fname, seq, labels=False):
    '''This function creates a FASTA file.'''
    with open(fname, 'wt') as fout: print(seq_create_fasta(seq,labels), file=fout)

def my_fasta_read(fname):
    '''Reads fasta file. 
    The functions returns a list of the sequences, a list of the sequences IDs 
    and a list of the sequences descriptions (content of the text following >>)'''
    seq = []
    iD = []
    desc = []  # description 
    for record in SeqIO.parse(fname, "fasta"):
        desc.append(record.description)
        seq.append(str(record.seq))
        iD.append(record.id)
    return seq, iD, desc

def get_seq_count(seqs):
    '''This function returns a dictionary, where the keys are the (unique)
    sequences and the values are the number of occurances. Input seqs is a list
    of strings.'''
    from collections import defaultdict
    d = defaultdict(int)
    for s in seqs: d[s] += 1
    #sorted(d, key=d.get, reverse=True) # this returns a sorted list in descending order
    return d 

def my_filter_df(df, field, val, cols=''):
    '''Given a Pandas dataframe, this function returns a sub-set dataframe where the value of
    the column "field" is equal to "val". "cols" (if give) is a list of column fields (strings) to return
    (otherwise, all columns are returned).'''
    return df.loc[df[field]==val, cols] if cols else df.loc[df[field]==val, :]