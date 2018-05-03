# GLOBAL IMPORTS
import os
import glob
from collections import defaultdict
import numpy as np


subset_list = ["WS", "WW", "SS", "SW"]

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


def dico_from_file(filename):
    tmp_dico = {}
    tmp_file = open(filename, "r")
    for line in tmp_file:
        split_line = line.split("=")
        if len(split_line) > 1:
            value = split_line[1].strip()
            try:
                tmp_dico[split_line[0]] = float(value)
            except:
                pass
    tmp_file.close()
    return tmp_dico


def extract_nuc_pct(fasta_path):
    ''' Compute nucleotides frequencies from fasta file '''
    nucindex = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}
    total = 0.0
    nbr_sites = 0
    nbr_species = 0
    fasta_file = open(fasta_path, 'r')
    for seq_unstriped in fasta_file:
        if seq_unstriped[0] != ">":
            nbr_species += 1
            seq = seq_unstriped.strip()
            assert len(seq) % 3 == 0
            nbr_sites = int(len(seq) / 3)
            for site in seq:
                if site in nucleotides:
                    nucindex[site] += 1.
                    total += 1.
    fasta_file.close()
    print("\tNumber of species: {0}".format(nbr_species))
    print("\tNumber of sites: {0}".format(nbr_sites))

    return {k: v / total for k, v in nucindex.items()}, nbr_sites, nbr_species


def append_or_init(dico, key, value):
    if key in dico:
        dico[key].append(value)
    else:
        dico[key] = [value]


def nested_dict_init():
    return defaultdict(nested_dict_init)


def weak_strong(nuc):
    if nuc == "A" or nuc == "T":
        return "W"
    else:
        return "S"


def format_hyphy_dico(hyphy_dico):
    if "pnCG" in hyphy_dico:
        hyphy_dico["pnC"] = hyphy_dico["pnCG"]
        hyphy_dico["pnG"] = hyphy_dico["pnCG"]
        hyphy_dico["pnA"] = 0.5 - hyphy_dico["pnCG"]
    hyphy_dico["pnT"] = 1.0 - (hyphy_dico["pnA"] + hyphy_dico["pnC"] + hyphy_dico["pnG"])

    if "epsA" in hyphy_dico and "epsM" not in hyphy_dico:
        hyphy_dico["epsM"] = 20 - sum([hyphy_dico["eps" + aa] for aa in amino_acids_set if aa != "M"])

    if 'w' not in hyphy_dico:
        codon_frequencies = np.ones(len(codons))
        for codon_index, codon in enumerate(codons):
            for nuc in codon:
                codon_frequencies[codon_index] *= hyphy_dico["pn" + nuc]
            epsilon = "eps" + codon_table[codon]
            if epsilon in hyphy_dico:
                codon_frequencies[codon_index] *= hyphy_dico[epsilon]

        codon_frequencies /= np.sum(codon_frequencies)

        d_dict, d0_dict = defaultdict(float), defaultdict(float)
        d, d0 = 0.0, 0.0

        for c_origin_id, codon_origin in enumerate(codons):
            for c_target_id, codon_target in enumerate(codons):
                mutations = [(n_from, n_to) for n_from, n_to in zip(codon_origin, codon_target) if n_from != n_to]
                if len(mutations) == 1 and codon_table[codon_origin] != codon_table[codon_target]:
                    nuc_origin, nuc_target = mutations[0]

                    d0_partiel = codon_frequencies[c_origin_id] * hyphy_dico["pn" + nuc_target]

                    exchan = "exch" + "".join(sorted((nuc_target, nuc_origin)))
                    if exchan in hyphy_dico:
                        d0_partiel *= hyphy_dico[exchan]

                    d_partiel = d0_partiel
                    beta = "b_" + "".join(sorted((codon_table[codon_origin], codon_table[codon_target])))
                    if beta in hyphy_dico:
                        if hyphy_dico[beta] > 100.0:
                            print("{0}={1} ({2} sites)".format(beta, hyphy_dico[beta], hyphy_dico["n"]))
                            hyphy_dico[beta] = 1.0
                        d_partiel *= hyphy_dico[beta]

                    omega_subset = "w_" + weak_strong(nuc_origin) + weak_strong(nuc_target)
                    if omega_subset in hyphy_dico:
                        d_partiel *= hyphy_dico[omega_subset]

                    epsilon = "eps" + codon_table[codon_target]
                    if epsilon in hyphy_dico:
                        d_partiel *= hyphy_dico[epsilon]

                    d += d_partiel
                    d0 += d0_partiel

                    if omega_subset not in hyphy_dico:
                        d_dict[omega_subset] += d_partiel
                        d0_dict[omega_subset] += d0_partiel

        hyphy_dico["w"] = d / d0
        for omega_subset in d_dict.keys():
            if d_dict[omega_subset] != 0.0 and (omega_subset not in hyphy_dico):
                hyphy_dico[omega_subset] = d_dict[omega_subset] / d0_dict[omega_subset]


def equilibrium_lambda(hyphy_dico):
    at_pct = equilibrium_at_pct(hyphy_dico)
    return at_pct / (1 - at_pct)


def equilibrium_at_pct(hyphy_dico):
    nbr_weak = np.array([c.count("A") + c.count("T") for c in codons])
    codon_frequencies = np.ones(len(codons))

    for index, codon in enumerate(codons):
        for nuc in codon:
            codon_frequencies[index] *= hyphy_dico["pn" + nuc]
        epsilon = "eps" + codon_table[codon]
        if epsilon in hyphy_dico:
            codon_frequencies[index] *= hyphy_dico[epsilon]

    codon_frequencies /= np.sum(codon_frequencies)

    return np.sum(codon_frequencies * nbr_weak) / 3

def prefs_path_to_list(preferences_path):
    preferences_list = []
    preferences_file = open(preferences_path, "r")
    preferences_file.readline()
    for line in preferences_file:
        preferences = list(map(float, line.strip().split(" ")[3:]))
        assert len(preferences) == 20
        preferences_list.append(preferences)
    preferences_file.close()
    return preferences_list


def theoretical_at_gc_pct(preferences_list, mut_bias):
    at_pct = theoretical_at_pct(preferences_list, mut_bias)
    return at_pct / (1 - at_pct)


def theoretical_at_pct(preferences_list, mut_bias):
    nbr_weak = np.array([c.count("A") + c.count("T") for c in codons])
    at_pct_list = []

    for preferences in preferences_list:
        codon_frequencies = np.power(mut_bias, nbr_weak)
        pref_codons = [preferences[codon_to_aa[i]] for i in range(len(codons))]
        codon_frequencies *= np.power(pref_codons, 1)
        codon_frequencies /= np.sum(codon_frequencies)

        at_pct_list.append(np.sum(codon_frequencies * nbr_weak) / 3)

    return np.mean(at_pct_list)


def projected_omega(preferences_list, mut_bias):
    nbr_weak = np.array([c.count("A") + c.count("T") for c in codons])

    codon_frequencies = np.ones((len(codons), len(preferences_list)))
    for site, preferences in enumerate(preferences_list):
        site_frequencies = np.power(mut_bias, nbr_weak)
        pref_codons = [preferences[codon_to_aa[i]] for i in range(len(codons))]
        site_frequencies *= np.power(pref_codons, 1)
        codon_frequencies[:, site] = site_frequencies / np.sum(site_frequencies)

    d, d0 = 0, 0
    for c_origin_id, codon_origin in enumerate(codons):
        mean_codon_frequency = np.mean(codon_frequencies[c_origin_id, :])

        for c_target_id, codon_target in enumerate(codons):
            mutations = [(n_from, n_to) for n_from, n_to in zip(codon_origin, codon_target) if n_from != n_to]
            if len(mutations) == 1 and codon_table[codon_origin] != codon_table[codon_target]:
                _, nuc_target = mutations[0]

                pij = 0
                for site, preferences in enumerate(preferences_list):
                    fi = np.log(preferences[codon_to_aa[c_origin_id]])
                    fj = np.log(preferences[codon_to_aa[c_target_id]])
                    site_pij = (fj - fi) / (1 - np.exp(fi - fj))
                    pij += codon_frequencies[c_origin_id, site] * site_pij

                pij /= np.sum(codon_frequencies[c_origin_id, :])

                d0_partiel = mean_codon_frequency

                if nuc_target == 'A' or nuc_target == 'T':
                    d0_partiel *= mut_bias

                d0 += d0_partiel
                d += d0_partiel * pij

    return d / d0