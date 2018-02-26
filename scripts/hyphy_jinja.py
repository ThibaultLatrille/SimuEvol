'''  
T. Latrille (Thibault.latrille@ens-lyon.org), adapted from SJS (stephanie.spielman@gmail.com, https://github.com/clauswilke/Omega_MutSel).
Script to generate Hyphy batchfile, for a variety of model parameterizations, from an alignment in .fasta and a tree.

USAGE: Usage: python prefs_to_freqs.py <fasta> <newick>. The first argument is the alignment file. The second argument is the newick tree file

Dependencies - numpy, jinja2


What does the code do?
First, we find the F61, F1x4, F3x4 frequencies. We use the global alignment frequencies.
Second, we set up the hyphy batch file which makes use of these frequencies.
Third, we generate the MG1 and MG3 matrix files.
'''

import sys
import numpy as np
import jinja2

codon_dict = {"AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N", "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
              "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S", "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",
              "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H", "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
              "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R", "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
              "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D", "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
              "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G", "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
              "TAC": "Y", "TAT": "Y", "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S", "TGC": "C", "TGG": "W",
              "TGT": "C", "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F"}

codons = sorted(codon_dict.keys())
nucindex = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def get_nuc_diff(source, target, grab_position=False):
    diff = ''
    position = 5
    for i in range(3):
        if source[i] != target[i]:
            diff += source[i] + target[i]
            position = i
    if grab_position:
        return diff, position
    else:
        return diff


def extract_codon_frequencies(fasta_path):
    ''' Compute codon frequencies from fasta file '''
    f61_freqs = np.zeros(61)

    fasta_file = open(fasta_path, 'r')

    codon_to_index = {codon: index for index, codon in enumerate(codons)}
    for seq_unstriped in fasta_file:
        seq = seq_unstriped.strip()
        if seq[0] != ">":
            assert len(seq) % 3 == 0
            for site in range(int(len(seq) / 3)):
                f61_freqs[codon_to_index[seq[3 * site:3 * site + 3]]] += 1

    f61_freqs /= sum(f61_freqs)
    return f61_freqs


def codon_to_nuc_freqs(codon_freqs):
    ''' Calculate nucleotide and positional nucleotide frequencies from codon frequencies. '''
    nuc_freqs = {k: 0 for k in nucindex.keys()}

    for index, codon in enumerate(codons):
        assert len(codon) == 3
        for nuc in codon:
            nuc_freqs[nuc] += codon_freqs[index] / 3

    return nuc_freqs


def extract_nuc_pct(fasta_path):
    ''' Compute codon frequencies from fasta file '''
    return codon_to_nuc_freqs(extract_codon_frequencies(fasta_path))


def build_emp_codon_freqs(nuc_freqs):
    ''' Compute F1x4 codon frequencies from nucleotide frequencies. '''
    codon_freqs = np.ones(61)
    for i, codon in enumerate(codons):
        for j in range(3):
            codon_freqs[i] *= nuc_freqs[codon[j]]
            
    pi_stop = nuc_freqs['T'] * nuc_freqs['A'] * nuc_freqs['G'] + \
              nuc_freqs['T'] * nuc_freqs['G'] * nuc_freqs['A'] + \
              nuc_freqs['T'] * nuc_freqs['A'] * nuc_freqs['A']
    codon_freqs /= (1. - pi_stop)
    assert (abs(np.sum(codon_freqs) - 1.) < 1e-8), "Could not properly calculate F1x4 frequencies."
    return codon_freqs


def build_codon_freqs(nuc_freqs):
    ''' Compute codon frequencies from GTR nucleotide frequencies. '''
    codon_freqs = [""] * 61
    for i, codon in enumerate(codons):
        codon_freqs[i] = "*".join([nuc_freqs[codon[j]] for j in range(3)])
        codon_freqs[i] += "/(1-p_stop)"
    return codon_freqs


def array_to_hyphy_freq(f):
    ''' Convert array of codon frequencies to a hyphy frequency string to directly write into a batchfile. '''
    hyphy_f = "{"
    for freq in f:
        hyphy_f += "{" + str(freq) + "},"
    hyphy_f = hyphy_f[:-1]
    hyphy_f += "};"
    return hyphy_f


def build_nuc_vars(nuc_freqs):
    ''' Compute codon frequencies from GTR nucleotide frequencies. '''
    const_freq = nuc_freqs['T']
    const_freq_sum = "+".join([var for var in nuc_freqs.values() if var != const_freq])
    nuc_globals = " ".join(["global {0}={1};".format(var, 1/len(nuc_freqs)) for var in nuc_freqs.values() if var != const_freq])
    nuc_globals += " global {0}:=1-({1});".format(const_freq, const_freq_sum)
    nuc_globals += " global p_stop:={0}; ".format(
        "+".join(["*".join([nuc_freqs[n] for n in stop]) for stop in ["TGA", "TAA", "TAG"]]))

    nuc_globals += " ".join(["{0}:>0; {0}:<1;".format(p) for p in nuc_freqs.values()])
    return nuc_globals


def is_TI(source, target):
    purines = ["A", "G"]
    pyrims = ["C", "T"]
    check1 = source in purines and target in purines
    check2 = source in pyrims and target in pyrims
    if check1 or check2:
        return True
    else:
        return False


def build_rates(param):
    assert param in [0, 1, 5]
    gtr_vars = {}
    for n_1 in nucindex.keys():
        for n_2 in nucindex.keys():
            if n_1 != n_2:
                key = n_1 + n_2
                if param == 1 and is_TI(n_1, n_2):
                    gtr_vars[key] = "k"
                if param == 5:
                    value = "r_" + "".join(sorted(n_1 + n_2))
                    if value != 'r_AC':
                        gtr_vars[key] = value
    return gtr_vars


def build_matrices(nuc_freqs, sym_vars):
    ''' Create MG94-style matrices (use target nucleotide frequencies).  '''
    matrix = '{61, 61, \n'  # MG94

    for i in range(61):
        source = codons[i]
        for j in range(61):
            target = codons[j]

            diff, x = get_nuc_diff(source, target, grab_position=True)
            if len(diff) == 2:
                assert (len(str(x)) == 1), "Problem with determining nucleotide difference between codons"
                freq = str(nuc_freqs[diff[1]])
                # Create string for matrix element
                element = '{' + str(i) + ',' + str(j) + ',t'
                if diff in sym_vars:
                    element += '*' + sym_vars[diff]
                if codon_dict[source] != codon_dict[target]:
                    element += '*w'
                matrix += element + '*' + freq + '} '

    # And save to file.
    return matrix + '}\n'


def build_hyphy_batchfile(raw_batch_path, batch_outfile,  fasta_infile, tree_infile,
                          rates=list([0, 5]), freqs=list([0, 1, 3])):
    # Parse input arguments and set up input/outfile files accordingly
    name = batch_outfile.split('/')[-1]
    raw_batchfile = raw_batch_path + "/globalDNDS_raw_jinja.bf"

    # Calculate frequency parameterizations
    print("Calculating frequency parameterizations.")
    nuc_freqs_dict, nuc_vars_dict, codon_freqs_dict = dict(), dict(), dict()
    nuc_freqs_dict[1] = {'A': 'p_at', 'C': 'p_cg', 'G': 'p_cg', 'T': 'p_at'}
    nuc_freqs_dict[3] = {'A': 'p_a', 'C': 'p_c', 'G': 'p_g', 'T': 'p_t'}

    for freq_param in freqs:
        assert freq_param in [0, 1, 3], "The number of free parameters for nucleotide frequencies is invalid"

        if freq_param == 0:
            nuc_freqs_dict[freq_param] = codon_to_nuc_freqs(extract_codon_frequencies(fasta_infile))
            codon_freqs = build_emp_codon_freqs(nuc_freqs_dict[freq_param])
        else:
            codon_freqs = build_codon_freqs(nuc_freqs_dict[freq_param])
            nuc_vars_dict[freq_param] = build_nuc_vars(nuc_freqs_dict[freq_param])

        codon_freqs_dict[freq_param] = array_to_hyphy_freq(codon_freqs)

    matrices_dict, rate_vars_dict = dict(), dict()
    for rate_param in rates:
        assert rate_param in [0, 1, 5], "The number of free parameters for symetric rates frequencies is invalid"
        rates = build_rates(rate_param)
        rate_vars_dict[rate_param] = " ".join(["global {0}=1.0; {0}:>0;".format(var) for var in set(rates.values())])

        matrices_dict[rate_param] = dict()
        for freq_param, freqs in nuc_freqs_dict.items():
            matrices_dict[rate_param][freq_param] = build_matrices(freqs, rates)

    # Create the hyphy batchfile to include the frequencies calculated here. Note that we need to do this since no
    # actual alignment exists which includes all protein positions, so cannot be read in by hyphy file.
    tree_file = open(tree_infile, 'r')
    tree = tree_file.readline()[:-2]
    tree_file.close()

    env = jinja2.Environment(loader=jinja2.FileSystemLoader("/"))
    template = env.get_template(raw_batchfile)
    t = template.stream(matrices_dict=matrices_dict,
                        codon_freqs_dict=codon_freqs_dict,
                        nuc_vars_dict=nuc_vars_dict,
                        rate_vars_dict=rate_vars_dict,
                        tree=tree, name=name,
                        fasta_infile=fasta_infile
                        ).dump(batch_outfile)

    # Use nucleotide and positional nucleotide frequencies to construct MG-style matrices
    print("Building and saving the MG-style matrices")
    return batch_outfile


def parse_input(arguments):
    usage_error = "\n\n Usage: python prefs_to_freqs.py <raw_batch> <batch> <fasta> <newick>.\n" \
                  "The first argument is the path to the dir containing the raw batch. " \
                  "The second argument is the path to the output batch. " \
                  "The third argument is the alignment file. " \
                  "The fourth argument is the newick tree file"
    assert (len(arguments) == 5), usage_error

    directory = arguments[1]
    batch = arguments[2]
    fasta = arguments[3]
    tree = arguments[4]

    return directory, batch, fasta, tree


if __name__ == '__main__':
    if len(sys.argv) == 5:
        raw, output, fasta, newick = parse_input(sys.argv)
    else:
        raw, output, fasta, newick = "/home/thibault/SimuEvol/", \
                                   "/home/thibault/SimuEvol/data_hyphy/npcst.bf", \
                                   "/home/thibault/SimuEvol/data_alignment/npcst.fasta", \
                                   "/home/thibault/SimuEvol/data_trees/np.newick"
    build_hyphy_batchfile(raw, output, fasta, newick, [0, 5], [0, 1, 3])
