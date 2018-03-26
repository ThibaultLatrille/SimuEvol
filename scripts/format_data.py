# GLOBAL IMPORTS
import os
import glob
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import lines
from matplotlib import markers

colors_array = list(colors.cnames.keys())
lines_array = list(lines.lineStyles.keys())
markers_array = list(markers.MarkerStyle.markers.keys())


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
    
    fasta_file = open(fasta_path, 'r')
    for seq_unstriped in fasta_file:
        if seq_unstriped[0] != ">":
            seq = seq_unstriped.strip()
            assert len(seq) % 3 == 0
            for site in seq:
                nucindex[site] += 1.
                total += 1.
    fasta_file.close()
    
    return {k: v/total for k, v in nucindex.items()}


def append_or_init(dico, key, value):
    if key in dico:
        dico[key].append(value)
    else:
        dico[key] = [value]


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


def nested_dict_init():
    return defaultdict(nested_dict_init)


def format_hyphy_dico(hyphy_dico):
    if "p_CG" in hyphy_dico:
        hyphy_dico["p_C"] = hyphy_dico["p_CG"]
        hyphy_dico["p_G"] = hyphy_dico["p_CG"]
        hyphy_dico["p_A"] = 0.5 - hyphy_dico["p_CG"]
    hyphy_dico["p_T"] = 1.0 - (hyphy_dico["p_A"] + hyphy_dico["p_C"] + hyphy_dico["p_G"])

    if "e_A" in hyphy_dico and "e_M" not in hyphy_dico:
        hyphy_dico["e_M"] = 20 - sum([hyphy_dico["e_" + aa] for aa in amino_acids_set if aa != "M"])

    if 'w' not in hyphy_dico:
        d, d0 = 0, 0
        codon_frequencies = np.ones(len(codons))
        for codon_index, codon in enumerate(codons):
            for nuc in codon:
                codon_frequencies[codon_index] *= hyphy_dico["p_" + nuc]
            epsilon = "e_" + codon_table[codon]
            if epsilon in hyphy_dico:
                codon_frequencies[codon_index] *= hyphy_dico[epsilon]

        codon_frequencies /= np.sum(codon_frequencies)

        for c_origin_id, codon_origin in enumerate(codons):
            for c_target_id, codon_target in enumerate(codons):
                mutations = [(n_from, n_to) for n_from, n_to in zip(codon_origin, codon_target) if n_from != n_to]
                if len(mutations) == 1 and codon_table[codon_origin] != codon_table[codon_target]:
                    nuc_origin, nuc_target = mutations[0]

                    d0_partiel = codon_frequencies[c_origin_id] * hyphy_dico["p_" + nuc_target]

                    exchan = "r_" + "".join(sorted((nuc_target, nuc_origin)))
                    if exchan in hyphy_dico:
                        d0_partiel *= hyphy_dico[exchan]

                    d_partiel = d0_partiel
                    beta = "b_" + "".join(sorted((codon_table[codon_origin], codon_table[codon_target])))
                    if beta in hyphy_dico:
                        d_partiel *= hyphy_dico[beta]
                    epsilon = "e_" + codon_table[codon_target]
                    if epsilon in hyphy_dico:
                        d_partiel *= hyphy_dico[epsilon]

                    d += d_partiel
                    d0 += d0_partiel

        hyphy_dico["w"] = d / d0


def equilibrium_lambda(hyphy_dico):
    at_pct = equilibrium_at_pct(hyphy_dico)
    return at_pct / (1 - at_pct)


def equilibrium_at_pct(hyphy_dico):
    nbr_weak = np.array([c.count("A") + c.count("T") for c in codons])
    codon_frequencies = np.ones(len(codons))

    for index, codon in enumerate(codons):
        for nuc in codon:
            codon_frequencies[index] *= hyphy_dico["p_"+nuc]
        epsilon = "e_" + codon_table[codon]
        if epsilon in hyphy_dico:
            codon_frequencies[index] *= hyphy_dico[epsilon]

    codon_frequencies /= np.sum(codon_frequencies)

    return np.sum(codon_frequencies * nbr_weak) / 3


nested_dict = nested_dict_init()
current_dir = "/home/thibault/SimuEvol/"
for file in sorted(os.listdir(current_dir + "/data_prefs")):
    prefix = file.strip().replace(".txt", "")
    hyphy_path = "{0}/data_hyphy/{1}".format(current_dir, prefix)
    if os.path.isdir(hyphy_path):
        param = prefix.split("_")
        nbr_sites = int(param[0])
        protein = param[1]
        mixture = param[2] == "True"
        alpha = float(param[3])

        protein_prefs_path = "{0}/data_prefs/{1}".format(current_dir, file)

        batchfiles = [batch[:-3] for batch in os.listdir(hyphy_path) if batch[-3:] == ".bf"]

        for batch in batchfiles:
            at_gc_pct = float(batch.split("_m")[-1])
            simu_evol_path = "{0}/data_alignment/{1}_{2}".format(current_dir, prefix, batch)
            simuevol_dico = dico_from_file(simu_evol_path + ".txt")

            nested_dict["$\omega$"]["$\omega$"][at_gc_pct] = simuevol_dico["w"]
            nested_dict["$\omega$"]["$\omega^0$"][at_gc_pct] = simuevol_dico["w0"]
            nested_dict["$\omega$"]["$<\omega^0>$"][at_gc_pct] = simuevol_dico["<w0>"]

            nuc_freqs = extract_nuc_pct(simu_evol_path + ".fasta")
            at_gc_pct_obs = (nuc_freqs['A'] + nuc_freqs['T']) / (nuc_freqs['G'] + nuc_freqs['C'])
            nested_dict["%$_{AT}$"]["Alignment"][at_gc_pct] = nuc_freqs['A'] + nuc_freqs['T']
            nested_dict["$\lambda$"]["Alignment"][at_gc_pct] = at_gc_pct_obs

            prefs_list = prefs_path_to_list(protein_prefs_path)
            nested_dict["%$_{AT}$"]["Theoretical MutSel"][at_gc_pct] = theoretical_at_pct(prefs_list, at_gc_pct)
            nested_dict["$\lambda$"]["Theoretical MutSel"][at_gc_pct] = theoretical_at_gc_pct(prefs_list, at_gc_pct)

            for hyphy_result in glob.glob("{0}/{1}*_hyout.txt".format(hyphy_path, batch)):
                experiment = hyphy_result.replace("_hyout.txt", "").split(".bf_")[-1]
                hyphy_dico = dico_from_file(hyphy_result)

                format_hyphy_dico(hyphy_dico)

                for param in hyphy_dico:
                    nested_dict[param][experiment][at_gc_pct] = hyphy_dico[param]

                nested_dict["$\omega$"][experiment][at_gc_pct] = hyphy_dico["w"]

                if ("p_G" in hyphy_dico) and ("p_G" in hyphy_dico):
                    gc_pct = hyphy_dico["p_G"] + hyphy_dico["p_C"]
                    nested_dict["%$_{AT}$"][experiment + " Mutation"][at_gc_pct] = 1 - gc_pct
                    nested_dict["$\lambda$"][experiment + " Mutation"][at_gc_pct] = (1 - gc_pct) / gc_pct

                nested_dict["%$_{AT}$"][experiment + " MutSel"][at_gc_pct] = equilibrium_at_pct(hyphy_dico)
                nested_dict["$\lambda$"][experiment + " MutSel"][at_gc_pct] = equilibrium_lambda(hyphy_dico)

        for param in ["%$_{AT}$", "$\lambda$", "$\omega$"]:
            nested_dict_l1 = nested_dict[param]
            my_dpi = 96
            fig = plt.figure(figsize=(880 / my_dpi, 400 / my_dpi), dpi=my_dpi)
            index = 0
            for experiment in sorted(nested_dict_l1.keys(), key=lambda x: x[::-1]):
                x_list = sorted(nested_dict_l1[experiment].keys())
                y_list = [nested_dict_l1[experiment][k] for k in x_list]
                plt.plot(x_list, y_list,
                         linestyle=':',
                         marker=markers_array[index % len(markers_array)],
                         label=experiment)
                index += 1
            plt.xscale('log')
            plt.xlabel('$\lambda$')
            if 'AT' in param:
                plt.plot(x_list, np.array(x_list) / (1 + np.array(x_list)), color="black", linewidth=0.5)
            if 'lambda' in param:
                plt.plot(x_list, x_list, color="black", linewidth=0.5)
                plt.yscale('log')
            elif 'omega' in param:
                plt.yscale('linear')
            plt.ylabel(param)
            plt.legend()
            plt.title(prefix)
            plt.tight_layout()
            plt.show()

        '''
        nbr_gc = np.zeros(20)
        for aa_index, aa in enumerate(amino_acids):
            for codon in aa_table[aa]:
                nbr_gc[aa_index] += codon.count('G') + codon.count('C')

            nbr_gc[aa_index] /= len(aa_table[aa]) * 3

        epsilon_dict = nested_dict_init()
        for experiment in nested_dict["e_A"]:
            for at_gc_pct in nested_dict["e_A"][experiment]:
                epsilon_dict[at_gc_pct][experiment] = np.zeros(20)
                for aa_index, aa in enumerate(amino_acids):
                    epsilon_dict[at_gc_pct][experiment][aa_index] = nested_dict["e_" + aa][experiment][at_gc_pct]

        for at_gc_pct in sorted(epsilon_dict):
            fig = plt.figure(figsize=(880 / my_dpi, 400 / my_dpi), dpi=my_dpi)
            plt.xlabel('%GC of amino-acids')
            plt.ylabel('$\epsilon$')
            for experiment in sorted(epsilon_dict[at_gc_pct]):
                plt.scatter(nbr_gc, epsilon_dict[at_gc_pct][experiment], label=str(experiment))
            plt.legend()
            plt.title(at_gc_pct)
            plt.tight_layout()
            plt.show()
        '''