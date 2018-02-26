# GLOBAL IMPORTS
import os
import glob
from collections import defaultdict
import numpy as np


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


def theoretical_at_gc_pct(preferences_list, mut_bias):
    nucleotides = "ACGT"
    codontable = {
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
    assert len(codontable.keys()) == 64, "There is not 64 codons in the codon table"

    codons = list([k for k, v in codontable.items() if v != 'X'])
    assert len(codons) == 61, "There is not 3 stop codons in the codon table"
    assert not [n for n in "".join(codons) if (n not in nucleotides)], "There is a codon with an unrecognized nucleotide"

    amino_acids_set = set(codontable.values())
    amino_acids_set.remove('X')

    nbr_weak = np.array([c.count("A") + c.count("T") for c in codons])
    nbr_strong = 3. - nbr_weak

    aa_char_to_int = {v: k for k, v in enumerate(sorted(amino_acids_set))}
    amino_acids = "".join(amino_acids_set)
    assert len(amino_acids) == 20, "There is not 21 amino-acids in the codon table"

    aa_table = {}
    codon_to_aa = [None] * len(codons)
    for index, codon in enumerate(codons):
        aa = codontable[codon]
        if aa not in aa_table:
            aa_table[aa] = []
        aa_table[aa].append(index)
        codon_to_aa[index] = aa_char_to_int[aa]
    assert len(aa_table) == 20, "There is not 21 amino-acids in the aa table"

    at_pct_list = []
    gc_pct_list = []

    for preferences in preferences_list:

        codon_frequencies = np.power(mut_bias, nbr_weak)
        pref_codons = [preferences[codon_to_aa[i]] for i in range(len(codons))]
        codon_frequencies *= np.power(pref_codons, 1)
        codon_frequencies /= np.sum(codon_frequencies)

        at_pct_list.append(np.sum(codon_frequencies * nbr_weak))
        gc_pct_list.append(np.sum(codon_frequencies * nbr_strong))

    return np.sum(at_pct_list) / np.sum(gc_pct_list)


def nested_dict_init():
    return defaultdict(nested_dict_init)


nested_dict = nested_dict_init()

current_dir = "/home/thibault/SimuEvol/"
for file in os.listdir(current_dir + "/data_prefs"):
    prefix = file.strip().replace(".txt", "")
    hyphy_path = "{0}/data_hyphy/{1}".format(current_dir, prefix)
    if os.path.isdir(hyphy_path):
        param = prefix.split("_")
        nbr_sites = int(param[0])
        protein = param[1]
        mixture = param[2] == "True"
        alpha = float(param[3])
        chain = int(param[4])

        protein_prefs_path = "{0}/data_prefs/{1}".format(current_dir, file)

        batchfiles = [batch[:-3] for batch in os.listdir(hyphy_path) if batch[-3:] == ".bf"]

        for batch in batchfiles:
            at_gc_pct = float(batch.split("_m")[-1])
            simu_evol_path = "{0}/data_alignment/{1}_{2}".format(current_dir, prefix, batch)
            simuevol_dico = dico_from_file(simu_evol_path + ".txt")

            nested_dict["Simulated"][at_gc_pct]["w"] = simuevol_dico["w"]
            nested_dict["Theoretical"][at_gc_pct]["w"] = simuevol_dico["w0"]

            nuc_freqs = extract_nuc_pct(simu_evol_path + ".fasta")
            at_gc_pct_obs = (nuc_freqs['A'] + nuc_freqs['T']) / (nuc_freqs['G'] + nuc_freqs['C'])
            nested_dict["Simulated"][at_gc_pct]["%AT/%GC"] = at_gc_pct_obs

            prefs_list = prefs_path_to_list(protein_prefs_path)
            nested_dict["Theoretical"][at_gc_pct]["%AT/%GC"] = theoretical_at_gc_pct(prefs_list, at_gc_pct)

            for hyphy_result in glob.glob("{0}/{1}*_hyout.txt".format(hyphy_path, batch)):
                experiment = hyphy_result.replace("_hyout.txt", "").split(".bf_")[-1]
                hyphy_dico = dico_from_file(hyphy_result)
                nested_dict[experiment][at_gc_pct]["w"] = hyphy_dico["w"]

                if "f3" in experiment:
                    gc_pct = hyphy_dico["p_g"] + hyphy_dico["p_c"]
                    nested_dict[experiment][at_gc_pct]["%AT/%GC"] = (1 - gc_pct) / gc_pct
                elif "f1" in experiment:
                    gc_pct = hyphy_dico["p_cg"]
                    nested_dict[experiment][at_gc_pct]["%AT/%GC"] = (1 - gc_pct) / gc_pct

                if "r5" in experiment:
                    for ex_param in ["r_AT", "r_CT", "r_CG", "r_GT", "r_AG"]:
                        nested_dict[experiment][at_gc_pct][ex_param] = hyphy_dico[ex_param]

        header = ["Experiment", "Lambda"]
        experiment_path = protein_prefs_path = "{0}/data_figures/{1}.tsv".format(current_dir, prefix)
        experiment_file = open(experiment_path, "w")
        parameter_set = set()
        for nested_dict_l1 in nested_dict.values():
            for nested_dict_l2 in nested_dict_l1.values():
                parameter_set = parameter_set.union(set(list(nested_dict_l2.keys())))

        parameters = sorted(list(parameter_set))
        header.extend(parameters)
        experiment_file.write("#" + "\t".join(header) + "\n")

        for experiment, nested_dict_l1 in nested_dict.items():
            for at_gc_pct in sorted(nested_dict_l1.keys()):
                val_list = [experiment, at_gc_pct]
                nested_dict_l2 = nested_dict_l1[at_gc_pct]
                for param in parameters:
                    if param in nested_dict_l2:
                        val_list.append(nested_dict_l2[param])
                    else:
                        val_list.append(1.)
                experiment_file.write("\t".join(map(str, val_list)) + "\n")
        experiment_file.close()
