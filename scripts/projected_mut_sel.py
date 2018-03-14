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


def array_to_hyphy_freq(f):
    ''' Convert array of codon frequencies to a hyphy frequency string to directly write into a batchfile. '''
    hyphy_f = "{"
    for freq in f:
        hyphy_f += "{" + str(freq) + "},"
    hyphy_f = hyphy_f[:-1]
    hyphy_f += "};"
    return hyphy_f


def build_nuc_vars(nuc_freqs, vars_dict):
    ''' Compute codon frequencies from GTR nucleotide frequencies. '''
    const_freq = nuc_freqs['T']
    values_set = set(nuc_freqs.values())
    assert len(values_set) == 2 or len(values_set) == 4, "There must be either 2 or 4 frequencies parameters"
    const_freq_sum = "+".join([var for var in values_set if var != const_freq])
    ratio = len(values_set) / len(nuc_freqs.values())

    for var in values_set:
        if var != const_freq:
            vars_dict[var] = "global {0}=.25; {0}:>0; {0}:<{1};".format(var, ratio)

    vars_dict[const_freq] = "global {0}:={1}-({2}); {0}:>0; {0}:<{1};".format(const_freq, ratio, const_freq_sum)
    return 1


def is_TI(source, target):
    purines = ["A", "G"]
    pyrims = ["C", "T"]
    check1 = source in purines and target in purines
    check2 = source in pyrims and target in pyrims
    if check1 or check2:
        return True
    else:
        return False


def build_rates(param, vars_dict):
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
                        vars_dict[value] = "global {0}=1.0; {0}:>0;".format(value)
    return gtr_vars


def build_matrices(nuc_freqs, exchan_vars, omega_param, vars_dict):
    ''' Create MG94-style matrices (use target nucleotide frequencies).  '''
    matrix = '{61, 61, \n'  # MG94
    codon_freqs = [""] * 61
    beta_set = set()
    for i, source in enumerate(codons):
        codon_freqs[i] = "*".join([nuc_freqs[source[j]] for j in range(3)])
        if omega_param == 95 or omega_param == 20:
            epsilon = "e_" + codon_dict[source]
            codon_freqs[i] += "*" + epsilon
            vars_dict[epsilon] = "global {0}=1.0; {0}:>0;".format(epsilon)

        for j, target in enumerate(codons):

            diff, x = get_nuc_diff(source, target, grab_position=True)
            if len(diff) == 2:
                assert (len(str(x)) == 1), "Problem with determining nucleotide difference between codons"
                freq = str(nuc_freqs[diff[1]])
                # Create string for matrix element
                element = '{' + str(i) + ',' + str(j) + ',mu*t'
                if diff in exchan_vars:
                    element += '*' + exchan_vars[diff]
                if codon_dict[source] != codon_dict[target]:
                    if omega_param == 1:
                        element += '*w'
                        vars_dict["w"] = "global w=1.0; w:>0;"
                    elif omega_param == 95 or omega_param == 20:
                        if omega_param == 95:
                            beta = 'w_' + "".join(sorted(codon_dict[source] + codon_dict[target]))
                            vars_dict[beta] = "global {0}=1.0; {0}:>0;".format(beta)
                            beta_set.add(beta)
                            element += '*' + beta
                        epsilon = 'e_' + codon_dict[target]
                        element += '*' + epsilon

                matrix += element + '*' + freq + '} '
    matrix += '}\n'

    print("{0} beta parameters out of {1:.0f} possible".format(len(beta_set), 20 * 19 / 2))
    z_var = "z_{0}".format(omega_param)
    vars_dict[z_var] = " global {0}:={1};".format(z_var, "+".join(codon_freqs))

    for i, codon in enumerate(codon_freqs):
        codon_freqs[i] = codon + "/" + z_var
    # And save to file.
    return matrix, array_to_hyphy_freq(codon_freqs)


def build_hyphy_batchfile(raw_batch_path, batch_outfile, fasta_infile, tree_infile,
                          rates=list([0, 1, 5]), freqs=list([0, 1, 3]), omega=list([1, 4])):
    # Parse input arguments and set up input/outfile files accordingly
    name = batch_outfile.split('/')[-1]
    raw_batchfile = raw_batch_path + "/projected_mut_sel.bf"

    # Calculate frequency parameterizations
    print("Calculating frequency parametrization.")
    nuc_freqs_dict = dict()
    nuc_freqs_dict[1] = {'A': 'p_at', 'C': 'p_cg', 'G': 'p_cg', 'T': 'p_at'}
    nuc_freqs_dict[3] = {'A': 'p_a', 'C': 'p_c', 'G': 'p_g', 'T': 'p_t'}

    matrices_nested_dict = nested_dict_init()
    freqs_nested_dict = nested_dict_init()
    vars_nested_dict = nested_dict_init()
    vars_all = dict()

    for rate_param in rates:
        for freq_param in freqs:
            for omega_param in omega:
                vars_nested_dict[rate_param][freq_param][omega_param] = dict()
                vars_dict = vars_nested_dict[rate_param][freq_param][omega_param]
                vars_dict["mu"] = "global mu=1; mu:>0;"

                nuc_freqs = nuc_freqs_dict[freq_param]
                build_nuc_vars(nuc_freqs_dict[freq_param], vars_dict)

                exchan_vars = build_rates(rate_param, vars_dict)

                matrix, codon_freqs = build_matrices(nuc_freqs, exchan_vars, omega_param, vars_dict)
                matrices_nested_dict[rate_param][freq_param][omega_param] = matrix

                if omega_param in freqs_nested_dict[freq_param]:
                    assert freqs_nested_dict[freq_param][omega_param] == codon_freqs
                else:
                    freqs_nested_dict[freq_param][omega_param] = codon_freqs

                for var in vars_dict.keys():
                    if var in vars_all:
                        assert vars_all[var] == vars_dict[var]
                    else:
                        vars_all[var] = vars_dict[var]

    # Create the hyphy batchfile to include the frequencies calculated here. Note that we need to do this since no
    # actual alignment exists which includes all protein positions, so cannot be read in by hyphy file.
    tree_file = Tree(tree_infile)
    tree = tree_file.write(format=0)

    env = jinja2.Environment(loader=jinja2.FileSystemLoader("/"))
    template = env.get_template(raw_batchfile)
    template.stream(matrices_nested_dict=matrices_nested_dict,
                    freqs_nested_dict=freqs_nested_dict,
                    vars_nested_dict=vars_nested_dict,
                    vars_all=vars_all,
                    tree=tree, name=name,
                    fasta_infile=fasta_infile).dump(batch_outfile)

    print("Building and saving Hyphy Batchfile ({0})".format(batch_outfile))
    return batch_outfile


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--directory', required=False, type=str,
                        default='/home/thibault/SimuEvol/scripts',
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
                        default="0-3-1_20_95",
                        dest="p", metavar="<parameters>",
                        help="The parameters of the GTR-Matrix")
    args = parser.parse_args()
    params = [list(map(int, p.split("_"))) for p in args.p.split("-")]
    build_hyphy_batchfile(args.d, args.b, args.f, args.t, params[0], params[1], params[2])
