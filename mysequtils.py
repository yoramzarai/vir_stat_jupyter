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

    
    
    
