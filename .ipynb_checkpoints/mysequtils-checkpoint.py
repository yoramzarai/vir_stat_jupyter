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