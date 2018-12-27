import sys
import getopt
import csv
import numpy
import weblogolib


#_alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

class Scales:
    Information = 1
    Probability = 2

def create_plot( data, basename, length, scale=Scales.Information, ):

    assert(length>1)
    assert(len(basename)>2)
    filename = ""
    options = weblogolib.LogoOptions()
    #options.title = "A Logo Title"
    #options.size = weblogolib.std_sizes['large']
    options.size = weblogolib.LogoSize( stack_height=200, stack_width=20 )
    options.first_index = -length
    options.show_errorbars = False
    options.title = "Stereotype for purified model 1001 (HSMM22s35c)"
    options.alphabet = weblogolib.unambiguous_protein_alphabet
    options.xaxis_label = "Position relative to C-terminal end"
    options.logo_title = "Amino-acid distribution in stereotype sequence (for purified model HSMM22/35c)"
    options.color_scheme = weblogolib.std_color_schemes["chemistry"]


    if( scale == Scales.Information ):
        options.yaxis_scale = 1.2
        options.yaxis_label = "AA emission probability and information (bits)"
        filename = "{}_information_scale".format(basename)

    elif (scale == Scales.Probability):
        options.unit_name = 'probability' # request probability scale (instead of information scale)
        options.yaxis_label = "Amino acid emission probability"
        filename = "{}_probability_scale".format(basename)
    else:
        raise "Invalid scale argument to create_plot()"


    format = weblogolib.LogoFormat(data, options)

    with open("{}.png".format(filename), 'w') as fout:
        weblogolib.png_formatter( data, format, fout)

    with open("{}.eps".format(filename), 'w') as fout:
        weblogolib.eps_formatter( data, format, fout)

    with open("{}.pdf".format(filename), 'w') as fout:
        weblogolib.pdf_formatter( data, format, fout)


def main_test(model_path, args = None):
    counts = []
    reader = csv.reader( open(model_path, 'rb'), delimiter=' ')
    for row in reader:
        vals = map(lambda x: float(x)*1000, filter(lambda x: len(x)>0, row) )
        counts.append(vals);

    counts_arr = numpy.array(counts).transpose()

    assert(counts_arr.shape[1] == 20 )

    #data = weblogolib.LogoData.from_counts(weblogolib.unambiguous_protein_alphabet, counts_arr, weblogolib.equiprobable_distribution(20) )
    data = weblogolib.LogoData.from_counts(weblogolib.unambiguous_protein_alphabet, counts_arr, None )

    for scale in [Scales.Probability, Scales.Information]:
        create_plot(data, "effectors.283.clean", counts_arr.shape[0], scale )



class UsageError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __init__(self):
        self.msg = "Usage: score_sequences <model.pickle> <sequences.fa>";


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.GetoptError as msg:
            raise UsageError(msg)

        if( len(argv) < 1 ):
            raise UsageError();

        main_test(argv[1], argv)

    except UsageError as err:
        print( err.msg )
        return 2;



if __name__ == "__main__":
    sys.exit(main())





