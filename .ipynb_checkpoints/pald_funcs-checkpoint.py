from Bio.Seq import Seq

def is_palindromic(lseq):
    '''Given a list of strings returns a list of 1 [0] if the
    corresponding string is [is not] palindromic'''
    assert isinstance(lseq, list), 'Input to is_palindromic() must be a list of strings'
    return [i == i.reverse_complement() for i in map(Seq,lseq)]

def get_palindromic(lseq):
    '''Given a list of strings returns a list of the palindromic strings'''
    assert isinstance(lseq, list), 'Input to is_palindromic() must be a list of strings'
    return [str(i) for i in map(Seq,lseq) if i == i.reverse_complement()]
