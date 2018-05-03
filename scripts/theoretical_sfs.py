# GLOBAL IMPORTS
import numpy as np
from scipy.special import binom
from scipy import integrate
import matplotlib.pyplot as plt

alpha = 0.1
at_gc_ratio = 1.0
theta = 0.1
n = 40
nbr_sites = 100

codon_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W'}

nucleotides = "ACGT"
assert len(codon_table.keys()) == 64, "There is not 64 codons in the codon table"

codons = list([k for k, v in codon_table.items() if v != 'X'])
assert len(codons) == 61, "There is not 3 stop codons in the codon table"
assert not [n for n in "".join(codons) if (n not in nucleotides)], "There is a codon with an unrecognized nucleotide"

amino_acids_set = set(codon_table.values())
amino_acids_set.remove('X')

aa_char_to_int = {v: k for k, v in enumerate(sorted(amino_acids_set))}
amino_acids = "".join(amino_acids_set)
assert len(amino_acids) == 20, "There is not 20 amino-acids in the codon table"

aa_table = {}
codon_to_aa = [None] * len(codons)
for index, codon in enumerate(codons):
    aa = codon_table[codon]
    if aa not in aa_table:
        aa_table[aa] = []
    aa_table[aa].append(codon)
    codon_to_aa[index] = aa_char_to_int[aa]
assert len(aa_table) == 20, "There is not 20 amino-acids in the aa table"


def generate_neighbors(synonymous=True):
    neighbor_dict = dict()
    for i, codon_origin in enumerate(codons):
        neighbor_dict[i] = []
        for j, codon_target in enumerate(codons):
            mutations = [(n_from, n_to) for n_from, n_to in zip(codon_origin, codon_target) if n_from != n_to]
            if len(mutations) == 1:
                if synonymous == (codon_table[codon_origin] == codon_table[codon_target]):
                    _, nuc_target = mutations[0]
                    neighbor_dict[i].append((j, nuc_target))
    return neighbor_dict


def dfe(site_preferences, mut_bias, synonymous):
    neighbors = generate_neighbors(synonymous=synonymous)

    nbr_weak = np.array([c.count("A") + c.count("T") for c in codons])
    mut_frequencies = np.power(mut_bias, nbr_weak)

    coefficients = list()
    weights = list()

    for preferences in site_preferences:
        pref_codons = np.array([preferences[codon_to_aa[i]] for i in range(len(codons))])
        codon_frequencies = mut_frequencies * pref_codons
        codon_frequencies /= np.sum(codon_frequencies)
        log_fitness = np.log(preferences)

        for i, codon_origin in enumerate(codons):
            for j, nuc_target in neighbors[i]:

                if synonymous:
                    coefficients.append(0.0)
                else:
                    coefficients.append(log_fitness[codon_to_aa[j]] - log_fitness[codon_to_aa[i]])

                weight = codon_frequencies[i]
                if nuc_target == 'A' or nuc_target == 'T':
                    weight *= mut_bias

                weights.append(weight)

    return np.array(coefficients), np.array(weights) / np.sum(weights)


def sfs(nbr_sample, site_preferences, mut_bias, theta, synonymous):
    neighbors = generate_neighbors(synonymous=synonymous)
    
    nbr_weak = np.array([c.count("A") + c.count("T") for c in codons])
    mut_frequencies = np.power(mut_bias, nbr_weak)

    sample_range = np.array(range(1, nbr_sample - 1))
    sample_sfs = np.zeros(len(sample_range))

    for preferences in site_preferences:
        pref_codons = np.array([preferences[codon_to_aa[i]] for i in range(len(codons))])
        codon_frequencies = mut_frequencies * pref_codons
        codon_frequencies /= np.sum(codon_frequencies)
        log_fitness = np.log(preferences)

        def residency(x):
            res = 0
            for i, codon_origin in enumerate(codons):
                tmp_res = 0
                for j, nuc_target in neighbors[i]:
                    pij = theta

                    if synonymous:
                        pij *= 1 - x
                    else:
                        s = log_fitness[codon_to_aa[j]] - log_fitness[codon_to_aa[i]]
                        if s == 0:
                            pij *= 1 - x
                        else:
                            pij *= (1 - np.exp(- s * (1 - x))) / (1 - np.exp(-s))

                    if nuc_target == 'A' or nuc_target == 'T':
                        pij *= mut_bias

                    tmp_res += pij

                res += codon_frequencies[i] * tmp_res
            return res

        precision = 8
        x_array = np.linspace(0, 1, 2 ** precision + 1)
        res_array = np.array([residency(x) for x in x_array])

        for index, a in enumerate(sample_range):
            y_array = res_array * np.power(x_array, a - 1) * np.power(1 - x_array, nbr_sample - a - 1)
            sample_sfs[index] += binom(nbr_sample, a) * integrate.simps(y_array, x_array)

    return sample_range, sample_sfs


bar_width = 0.35
opacity = 1
RED = "#EB6231"
YELLOW = "#E29D26"
BLUE = "#5D80B4"
LIGHTGREEN = "#6ABD9B"
GREEN = "#8FB03E"

propensities = np.random.dirichlet(alpha * np.ones(20), nbr_sites)

coeffs, weights = dfe(propensities, at_gc_ratio, synonymous=False)
plt.hist(coeffs, bins=50, density=True, weights=weights,
         label='Non-synonymous ($\\sum = {0:.3g}$)'.format(sum(weights)))

plt.legend()
plt.xlabel("Selection coefficient")
plt.ylabel("Density")
plt.title('DFE for $\\lambda = {0}$, $n = {1}$ '
          'and $\\alpha = {2}$'.format(at_gc_ratio, n, alpha))
plt.show()

sampling, non_syn_sfs = sfs(n, propensities, at_gc_ratio, theta, synonymous=False)
plt.bar(sampling, non_syn_sfs, bar_width,
        label='Non-synonymous ($\\sum = {0:.3g}$)'.format(sum(non_syn_sfs)),
        alpha=opacity, color=BLUE)
sampling, syn_sfs = sfs(n, propensities, at_gc_ratio, theta, synonymous=True)
plt.bar(sampling + bar_width, syn_sfs, bar_width,
        label='Synonymous ($\\sum = {0:.3g}$)'.format(sum(syn_sfs)),
        alpha=opacity, color=LIGHTGREEN)

plt.legend()
plt.xlabel("At derived frequency")
plt.ylabel("Expected number of SNPs")
plt.title('SFS for $\\lambda = {0}$, $n = {1}$ '
          'and $\\alpha = {2}$'.format(at_gc_ratio, n, alpha))
plt.show()
