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

import numpy as np
import jinja2
from collections import defaultdict
import argparse
from ete3 import Tree

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


def nested_dict_init():
    return defaultdict(nested_dict_init)


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


def p_stop_var(nuc_freqs):
    return "p_stop_{0}".format(len(set(nuc_freqs.values())))


def build_codon_freqs(nuc_freqs):
    ''' Compute codon frequencies from GTR nucleotide frequencies. '''
    codon_freqs = [""] * 61
    for i, codon in enumerate(codons):
        codon_freqs[i] = "*".join([nuc_freqs[codon[j]] for j in range(3)])
        codon_freqs[i] += "/(1-{0})".format(p_stop_var(nuc_freqs))
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
    values_set = set(nuc_freqs.values())
    assert len(values_set) == 2 or len(values_set) == 4, "There must either 2 or 4 frequencies parameters"
    nuc_globals = " ".join(["global {0}=.25;".format(var) for var in values_set if var != const_freq])

    const_freq_sum = "+".join([var for var in values_set if var != const_freq])
    ratio = len(values_set) / len(nuc_freqs.values())
    nuc_globals += " global {0}:={1}-({2});".format(const_freq, ratio, const_freq_sum)

    p_stop = "+".join(["*".join([nuc_freqs[n] for n in stop]) for stop in ["TGA", "TAA", "TAG"]])
    nuc_globals += " global {0}:={1}; ".format(p_stop_var(nuc_freqs), p_stop)

    nuc_globals += " ".join(["{0}:>0; {0}:<{1};".format(p, ratio) for p in values_set])
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


def weak_strong(source, target):
    weak = ["A", "T"]
    out_str = ""
    for nuc in [source, target]:
        if nuc in weak:
            out_str += "W"
        else:
            out_str += "S"
    return out_str


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


def build_omega_rates(param):
    assert param in [1, 4, 6]
    gtr_omega_vars = {}
    for n_1 in nucindex.keys():
        for n_2 in nucindex.keys():
            if n_1 != n_2:
                key = n_1 + n_2
                if param == 1:
                    gtr_omega_vars[key] = "w"
                if param == 4:
                    gtr_omega_vars[key] = "w_" + weak_strong(n_1, n_2)
                elif param == 6:
                    gtr_omega_vars[key] = "w_" + "".join(sorted(n_1 + n_2))
    return gtr_omega_vars


def build_matrices(nuc_freqs, exchan_vars, omega_vars):
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
                if diff in exchan_vars:
                    element += '*' + exchan_vars[diff]
                if codon_dict[source] != codon_dict[target]:
                    element += '*' + omega_vars[diff]
                matrix += element + '*' + freq + '} '

    # And save to file.
    return matrix + '}\n'


def build_hyphy_batchfile(raw_batch_path, batch_outfile,  fasta_infile, tree_infile,
                          rates=list([0, 1, 5]), freqs=list([0, 1, 3]), omega=list([1, 4])):
    # Parse input arguments and set up input/outfile files accordingly
    name = batch_outfile.split('/')[-1]
    raw_batchfile = raw_batch_path + "/globalDNDS_raw_jinja.bf"

    # Calculate frequency parameterizations
    print("Calculating frequency parametrization.")
    nuc_freqs_dict, codon_freqs_dict = dict(), dict()
    nuc_freqs_dict[1] = {'A': 'p_at', 'C': 'p_cg', 'G': 'p_cg', 'T': 'p_at'}
    nuc_freqs_dict[3] = {'A': 'p_a', 'C': 'p_c', 'G': 'p_g', 'T': 'p_t'}

    nuc_vars_dict, rate_vars_dict, omega_vars_dict = dict(), dict(), dict()
    for freq_param in freqs:
        assert freq_param in [0, 1, 3], "The number of free parameters for nucleotide frequencies is invalid"

        if freq_param == 0:
            nuc_freqs_dict[freq_param] = codon_to_nuc_freqs(extract_codon_frequencies(fasta_infile))
            codon_freqs = build_emp_codon_freqs(nuc_freqs_dict[freq_param])
        else:
            codon_freqs = build_codon_freqs(nuc_freqs_dict[freq_param])
            nuc_vars_dict[freq_param] = build_nuc_vars(nuc_freqs_dict[freq_param])

        codon_freqs_dict[freq_param] = array_to_hyphy_freq(codon_freqs)

    omega_rates_dict = dict()
    for omega_param in omega:
        assert omega_param in [1, 4, 6], "The number of free parameters for omega is invalid"
        omega_rates_dict[omega_param] = build_omega_rates(omega_param)
        omega_rates_set = set(omega_rates_dict[omega_param].values())
        omega_vars_dict[omega_param] = " ".join(["global {0}=1.0; {0}:>0;".format(var) for var in omega_rates_set])

    matrices_dict = nested_dict_init()
    for rate_param in rates:
        assert rate_param in [0, 1, 5], "The number of free parameters for symmetric rates frequencies is invalid"
        rates = build_rates(rate_param)
        rate_vars_dict[rate_param] = " ".join(["global {0}=1.0; {0}:>0;".format(var) for var in set(rates.values())])

        for freq_param, freqs in nuc_freqs_dict.items():
            for omega_param, omega_rates in omega_rates_dict.items():
                matrices_dict[rate_param][freq_param][omega_param] = build_matrices(freqs, rates, omega_rates)

    # Create the hyphy batchfile to include the frequencies calculated here. Note that we need to do this since no
    # actual alignment exists which includes all protein positions, so cannot be read in by hyphy file.
    tree_file = Tree(tree_infile)
    tree = tree_file.write(format=9)

    env = jinja2.Environment(loader=jinja2.FileSystemLoader("/"))
    template = env.get_template(raw_batchfile)
    t = template.stream(matrices_dict=matrices_dict,
                        codon_freqs_dict=codon_freqs_dict,
                        nuc_vars_dict=nuc_vars_dict,
                        rate_vars_dict=rate_vars_dict,
                        omega_vars_dict=omega_vars_dict,
                        tree=tree, name=name,
                        fasta_infile=fasta_infile
                        ).dump(batch_outfile)

    # Use nucleotide and positional nucleotide frequencies to construct MG-style matrices
    print("Building and saving the MG-style matrices")
    return batch_outfile


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--directory', required=False, type=str,
                        default='/home/thibault/panhome/SimuEvol/scripts',
                        dest="d", metavar="<dir_raw_batch>",
                        help="The path to the directory containing the raw batch (globalDNDS_raw_jinja.bf)")
    parser.add_argument('-b', '--batch', required=False, type=str,
                        default='/home/thibault/SimuEvol/data_hyphy/npcst.bf',
                        dest="b", metavar="<batch_output_path.bf>",
                        help="The path to the output batch.")
    parser.add_argument('-f', '--fasta', required=False, type=str,
                        default="/home/thibault/SimuEvol/data_alignment/npcst.fasta",
                        dest="f", metavar="<.fasta>",
                        help="The path to fasta alignment file")
    parser.add_argument('-t', '--tree', required=False, type=str,
                        default="/home/thibault/SimuEvol/data_trees/np.newick",
                        dest="t", metavar="<.newick>",
                        help="The path to the newick tree file")
    parser.add_argument('-p', '--parameters', required=False, type=str,
                        default="0_5-0_1_3-1_4_6",
                        dest="p", metavar="<parameters>",
                        help="The parameters of the GTR-Matrix")
    args = parser.parse_args()
    params = [list(map(int, p.split("_"))) for p in args.p.split("-")]
    build_hyphy_batchfile(args.d, args.b, args.f, args.t, params[0], params[1], params[2])
