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
amino_acids = sorted(set(codon_dict.values()))
nucindex = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
nuc_vars = {'A': 'p_a', 'C': 'p_c', 'G': 'p_g', 'T': 'p_t'}
const_freq = 'p_t'
const_freq_sum = "+".join([var for var in nuc_vars.values() if var != const_freq])

nuc_globals = " ".join(["global {0}=0.25;".format(var) for var in nuc_vars.values() if var != const_freq])
nuc_globals += " global {0}:=1-({1});".format(const_freq, const_freq_sum)
nuc_globals += " global p_stop:={0}; ".format(
    "+".join(["*".join([nuc_vars[n] for n in stop]) for stop in ["TGA", "TAA", "TAG"]]))

nuc_constrains = " ".join(["{0}:>0; {0}:<1;".format(p) for p in nuc_vars.values()])
var_mutational = nuc_globals + nuc_constrains

gtr_vars = {}
for n_1 in nucindex.keys():
    for n_2 in nucindex.keys():
        if n_1 != n_2:
            key = n_1 + n_2
            value = "s_" + "".join(sorted(n_1 + n_2))
            if value != 's_AC':
                gtr_vars[key] = value
var_empirical = " ".join(["global {0}=1.0; {0}:>0;".format(var) for var in set(gtr_vars.values())])

purines = ["A", "G"]
pyrims = ["C", "T"]
weak = ["A", "T"]
strong = ["A", "T"]


def extract_nuc_pct(fasta_path):
    ''' Compute codon frequencies from fasta file '''
    nuc_freqs = {k: 0.0 for k in nucindex.keys()}

    fasta_file = open(fasta_path, 'r')
    total = 0.
    for seq in fasta_file:
        seq = seq.strip()
        if seq[0] != ">":
            for n in seq:
                nuc_freqs[n] += 1.0
                total += 1.0

    for k in nuc_freqs.keys():
        nuc_freqs[k] /= total
    return nuc_freqs


def extract_frequencies(fasta_path):
    ''' Compute codon frequencies from fasta file '''
    f61_freqs = np.zeros(61)

    fasta_file = open(fasta_path, 'r')

    codon_to_index = {codon: index for index, codon in enumerate(codons)}
    for seq in fasta_file:
        if seq[0] != ">":
            for site in range(int(len(seq) / 3)):
                f61_freqs[codon_to_index[seq[3 * site:3 * site + 3]]] += 1

    f61_freqs /= sum(f61_freqs)
    return f61_freqs


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


def codon_to_nuc(codon_freqs):
    ''' Calculate nucleotide and positional nucleotide frequencies from codon frequencies. '''
    pos_nuc_freqs = np.zeros([4, 3])
    for i in range(3):
        for count, codon in enumerate(codons):
            pos_nuc_freqs[nucindex[codon[i]]][i] += codon_freqs[count]
    assert (np.allclose(np.sum(pos_nuc_freqs, axis=0), np.ones(3))), "bad positional nucleotide frequencies"
    return np.mean(pos_nuc_freqs, axis=1)


def calc_f1x4_freqs(nuc_freqs):
    ''' Compute F1x4 codon frequencies from nucleotide frequencies. '''
    f1x4 = np.ones(61)
    pi_stop = nuc_freqs[nucindex['T']] * nuc_freqs[nucindex['A']] * nuc_freqs[nucindex['G']] + \
              nuc_freqs[nucindex['T']] * nuc_freqs[nucindex['G']] * nuc_freqs[nucindex['A']] + \
              nuc_freqs[nucindex['T']] * nuc_freqs[nucindex['A']] * nuc_freqs[nucindex['A']]
    for i, codon in enumerate(codons):
        for j in range(3):
            f1x4[i] *= nuc_freqs[nucindex[codon[j]]]
    f1x4 /= (1. - pi_stop)
    assert (abs(np.sum(f1x4) - 1.) < 1e-8), "Could not properly calculate F1x4 frequencies."
    return f1x4


def calc_gtr_freqs():
    ''' Compute codon frequencies from GTR nucleotide frequencies. '''
    gtr_freqs = [""] * 61
    for i, codon in enumerate(codons):
        gtr_freqs[i] = "*".join([nuc_vars[codon[j]] for j in range(3)])
        gtr_freqs[i] += "/(1-p_stop)"
    return gtr_freqs


def build_mg_matrices(nuc_freqs):
    ''' Create MG94-style matrices (use target nucleotide frequencies).  '''

    matrix_emp = 'EmpMatrix = {61, 61, \n'  # MG94
    matrix_mut = 'MutMatrix = {61, 61, \n'  # MG94 with 4 omega parameters weak-to-strong

    for i in range(61):
        source = codons[i]
        for j in range(61):
            target = codons[j]

            diff, x = get_nuc_diff(source, target, grab_position=True)
            if len(diff) == 2:
                assert (len(str(x)) == 1), "Problem with determining nucleotide difference between codons"
                emp = str(nuc_freqs[nucindex[diff[1]]])
                mut = nuc_vars[diff[1]]
                # Create string for matrix element
                element = '{' + str(i) + ',' + str(j) + ',t'
                if diff in gtr_vars:
                    element += '*' + gtr_vars[diff]
                if codon_dict[source] != codon_dict[target]:
                    element += '*w'
                matrix_emp += element + '*' + emp + '} '
                matrix_mut += element + '*' + mut + '} '

    # And save to file.
    return matrix_emp + '};\n\n' + matrix_mut + '};\n\n'


def array_to_hyphy_freq(f):
    ''' Convert array of codon frequencies to a hyphy frequency string to directly write into a batchfile. '''
    hyphy_f = "{"
    for freq in f:
        hyphy_f += "{" + str(freq) + "},"
    hyphy_f = hyphy_f[:-1]
    hyphy_f += "};"
    return hyphy_f


def parse_input(arguments):
    usage_error = "\n\n Usage: python prefs_to_freqs.py <dir> <fasta> <newick>.\n" \
                  "The first argument is the output directory. " \
                  "The second argument is the alignment file. " \
                  "The third argument is the newick tree file"
    assert (len(arguments) == 4), usage_error

    directory = arguments[1]
    fasta = arguments[2]
    tree = arguments[3]

    return directory, fasta, tree


def build_hyphy_batchfile(hyphy_dir, fasta_infile, tree_infile):
    # Parse input arguments and set up input/outfile files accordingly
    name = fasta_infile.split('/')[-1]
    raw_batchfile = hyphy_dir + "/scripts/globalDNDS_raw_jinja.bf"
    batch_outfile = hyphy_dir + "/data_hyphy/globalDNDS_" + name + '.bf'

    # Calculate frequency parameterizations
    print("Calculating frequency parameterizations.")

    nuc_freqs = codon_to_nuc(extract_frequencies(fasta_infile))
    f1x4_freqs = calc_f1x4_freqs(nuc_freqs)

    gtr_freqs = calc_gtr_freqs()

    # Convert frequency arrays to hyphy strings
    freqs_mutational = array_to_hyphy_freq(gtr_freqs)
    freqs_empirical = array_to_hyphy_freq(f1x4_freqs)

    # Create the hyphy batchfile to include the frequencies calculated here. Note that we need to do this since no
    # actual alignment exists which includes all protein positions, so cannot be read in by hyphy file.
    tree_file = open(tree_infile, 'r')
    tree = tree_file.readline()[:-2]
    tree_file.close()

    matrices = build_mg_matrices(nuc_freqs)
    env = jinja2.Environment(loader=jinja2.FileSystemLoader("/"))
    template = env.get_template(raw_batchfile)
    t = template.stream(matrices=matrices,
                        freqs_empirical=freqs_empirical,
                        freqs_mutational=freqs_mutational,
                        var_empirical=var_empirical,
                        var_mutational=var_mutational,
                        tree=tree, name=name,
                        fasta_infile=fasta_infile
                        ).dump(batch_outfile)

    # Use nucleotide and positional nucleotide frequencies to construct MG-style matrices
    print("Building and saving the MG-style matrices")
    return batch_outfile


if __name__ == '__main__':
    if len(sys.argv) == 4:
        directory, fasta, newick = parse_input(sys.argv)
    else:
        directory, fasta, newick = "/home/thibault/SimuEvol/", \
                                   "/home/thibault/SimuEvol/data_alignment/npcst.fasta", \
                                   "/home/thibault/SimuEvol/data_trees/np.newick"
    build_hyphy_batchfile(directory, fasta, newick)
