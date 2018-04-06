# GLOBAL IMPORTS
import os
import glob
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import lines
from matplotlib import markers
import statsmodels.api as sm
from math import ceil


def interactive_legend(ax=None):
    if ax is None:
        ax = plt.gca()
    if ax.legend_ is None:
        ax.legend()

    return InteractiveLegend(ax.legend_)


class InteractiveLegend(object):
    def __init__(self, legend):
        self.legend = legend
        self.fig = legend.axes.figure

        self.lookup_artist, self.lookup_handle = self._build_lookups(legend)
        self._setup_connections()

        self.update()

    def _setup_connections(self):
        for artist in self.legend.texts + self.legend.legendHandles:
            artist.set_picker(10)  # 10 points tolerance

        self.fig.canvas.mpl_connect('pick_event', self.on_pick)
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)

    def _build_lookups(self, legend):
        labels = [t.get_text() for t in legend.texts]
        handles = legend.legendHandles
        label2handle = dict(zip(labels, handles))
        handle2text = dict(zip(handles, legend.texts))

        lookup_artist = {}
        lookup_handle = {}
        for artist in legend.axes.get_children():
            if artist.get_label() in labels:
                handle = label2handle[artist.get_label()]
                lookup_handle[artist] = handle
                lookup_artist[handle] = artist
                lookup_artist[handle2text[handle]] = artist

        lookup_handle.update(zip(handles, handles))
        lookup_handle.update(zip(legend.texts, handles))

        return lookup_artist, lookup_handle

    def on_pick(self, event):
        handle = event.artist
        if handle in self.lookup_artist:
            artist = self.lookup_artist[handle]
            artist.set_visible(not artist.get_visible())
            self.update()

    def on_click(self, event):
        if event.button == 3:
            visible = False
        elif event.button == 2:
            visible = True
        else:
            return

        for artist in self.lookup_artist.values():
            artist.set_visible(visible)
        self.update()

    def update(self):
        for artist in self.lookup_artist.values():
            handle = self.lookup_handle[artist]
            if artist.get_visible():
                handle.set_visible(True)
            else:
                handle.set_visible(False)
        self.fig.canvas.draw()

    def show(self):
        plt.show()


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

    return {k: v / total for k, v in nucindex.items()}


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

        d, d0 = {"All": .0}, {"All": .0}

        for key in subset_list:
            d[key] = 0.0
            d0[key] = 0.0

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
                        d_partiel *= hyphy_dico[beta]

                    epsilon = "eps" + codon_table[codon_target]
                    if epsilon in hyphy_dico:
                        d_partiel *= hyphy_dico[epsilon]

                    d['All'] += d_partiel
                    d0['All'] += d0_partiel

                    if beta in hyphy_dico:
                        subset = weak_strong(nuc_origin) + weak_strong(nuc_target)
                        d[subset] += d_partiel
                        d0[subset] += d0_partiel

        hyphy_dico["w"] = d['All'] / d0['All']
        for subset in subset_list:
            if d[subset] != 0.0:
                hyphy_dico["w_" + subset] = d[subset] / d0[subset]


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


current_dir = "/home/thibault/SimuEvol/"
for file in sorted(os.listdir(current_dir + "/data_prefs")):
    prefix = file.strip().replace(".txt", "")
    hyphy_path = "{0}/data_hyphy/{1}".format(current_dir, prefix)
    if os.path.isdir(hyphy_path):
        nested_dict = nested_dict_init()
        protein_prefs_path = "{0}/data_prefs/{1}".format(current_dir, file)
        batchfiles = [batch[:-3] for batch in os.listdir(hyphy_path) if batch[-3:] == ".bf"]


        for batch in batchfiles:
            at_gc_pct = float(batch.split("_m")[-1].split("_")[0])
            simu_evol_path = "{0}/data_alignment/{1}_{2}".format(current_dir, prefix, "_".join(batch.split("_")[:-1]))
            simuevol_dico = dico_from_file(simu_evol_path + ".txt")

            nested_dict["$\omega$"]["$\omega$"][at_gc_pct] = simuevol_dico["w"]
            nested_dict["$\omega$"]["$\omega^0$"][at_gc_pct] = simuevol_dico["w0"]
            nested_dict["$\omega$"]["$<\omega^0>$"][at_gc_pct] = simuevol_dico["<w0>"]

            for subset in subset_list:
                nested_dict["$\omega_{subset}$"]["$\omega_{" + subset + "}$"][at_gc_pct] = simuevol_dico["w_" + subset]
                nested_dict["$\omega_{subset}$"]["$\omega^0_{" + subset + "}$"][at_gc_pct] = simuevol_dico["w0_" + subset]

            nuc_freqs = extract_nuc_pct(simu_evol_path + ".fasta")
            at_gc_pct_obs = (nuc_freqs['A'] + nuc_freqs['T']) / (nuc_freqs['G'] + nuc_freqs['C'])
            nested_dict["%$_{AT}$"]["Alignment"][at_gc_pct] = nuc_freqs['A'] + nuc_freqs['T']
            nested_dict["$\lambda$"]["Alignment"][at_gc_pct] = at_gc_pct_obs

            prefs_list = prefs_path_to_list(protein_prefs_path)
            # nested_dict["$\omega$"]["$\omega_P$"][at_gc_pct] = projected_omega(prefs_list, at_gc_pct)
            nested_dict["%$_{AT}$"]["Theoretical Mut-Sel"][at_gc_pct] = theoretical_at_pct(prefs_list, at_gc_pct)
            nested_dict["$\lambda$"]["Theoretical Mut-Sel"][at_gc_pct] = theoretical_at_gc_pct(prefs_list, at_gc_pct)

            for hyphy_result in glob.glob("{0}/{1}*_hyout.txt".format(hyphy_path, batch)):
                experiment = hyphy_result.replace(".bf_hyout.txt", "").split("_")[-1]
                hyphy_dico = dico_from_file(hyphy_result)
                format_hyphy_dico(hyphy_dico)

                for param in hyphy_dico:
                    nested_dict[param][experiment][at_gc_pct] = hyphy_dico[param]

                nested_dict["$\omega$"][experiment][at_gc_pct] = hyphy_dico["w"]

                for subset in subset_list:
                    omega = "w_" + subset
                    if omega in hyphy_dico:
                        nested_dict["$\omega_{subset}$"][experiment + subset][at_gc_pct] = hyphy_dico[omega]

                if ("pnG" in hyphy_dico) and ("pnG" in hyphy_dico):
                    gc_pct = hyphy_dico["pnG"] + hyphy_dico["pnC"]
                    nested_dict["%$_{AT}$"][experiment + " Mutation"][at_gc_pct] = 1 - gc_pct
                    nested_dict["$\lambda$"][experiment + " Mutation"][at_gc_pct] = (1 - gc_pct) / gc_pct

                nested_dict["%$_{AT}$"][experiment + " Mut-Sel"][at_gc_pct] = equilibrium_at_pct(hyphy_dico)
                nested_dict["$\lambda$"][experiment + " Mut-Sel"][at_gc_pct] = equilibrium_lambda(hyphy_dico)

        for sub, param in enumerate(["$\lambda$", "$\omega$", "$\omega_{subset}$"]):
            nested_dict_l1 = nested_dict[param]
            my_dpi = 96
            fig, ax = plt.subplots()
            index = 0
            for experiment in sorted(nested_dict_l1.keys(), key=lambda x: x):
                x_list = sorted(nested_dict_l1[experiment].keys())
                y_list = [nested_dict_l1[experiment][k] for k in x_list]
                ax.plot(x_list, y_list,
                        linestyle=['-', '--'][index % 2], label=experiment, linewidth=2,
                        marker=markers_array[index % len(markers_array)])
                index += 1
            ax.set_xscale('log')
            ax.set_xlabel('$\lambda$')
            if 'AT' in param:
                ax.plot(x_list, np.array(x_list) / (1 + np.array(x_list)), color="black", linewidth=0.5)
            if 'lambda' in param:
                ax.plot(x_list, x_list, color="black", linewidth=0.5)
                ax.set_yscale('log')
            elif 'omega' in param:
                ax.set_yscale('linear')
            ax.set_ylabel(param)
            ax.legend()
            ax.set_title(prefix)
            interactive_legend().show()

        nbr_gc = np.zeros(20)
        for aa_index, aa in enumerate(amino_acids):
            for codon in aa_table[aa]:
                nbr_gc[aa_index] += codon.count('G') + codon.count('C')

            nbr_gc[aa_index] /= len(aa_table[aa]) * 3

        epsilon_dict = nested_dict_init()
        for experiment in nested_dict["epsA"]:
            for at_gc_pct in nested_dict["epsA"][experiment]:
                epsilon_dict[at_gc_pct][experiment] = np.zeros(20)
                for aa_index, aa in enumerate(amino_acids):
                    epsilon_dict[at_gc_pct][experiment][aa_index] = nested_dict["eps" + aa][experiment][at_gc_pct]

        RED = "#EB6231"
        YELLOW = "#E29D26"
        BLUE = "#5D80B4"
        LIGHTGREEN = "#6ABD9B"
        GREEN = "#8FB03E"
        colors = [YELLOW, BLUE, GREEN, RED, LIGHTGREEN]

        fontsize = 6
        dim = int(ceil(np.sqrt(len(epsilon_dict))))
        fig, axes = plt.subplots(dim, dim, sharex='all', sharey='all')
        for sub, at_gc_pct in enumerate(sorted(epsilon_dict)):
            ax = axes[int(sub / dim)][sub % dim]
            ax.set_xlabel('%GC of amino-acids', fontsize=fontsize)
            ax.set_ylabel('$\epsilon$', fontsize=fontsize)
            for index, experiment in enumerate(sorted(epsilon_dict[at_gc_pct])):
                ax.scatter(nbr_gc, epsilon_dict[at_gc_pct][experiment],
                           label=str(experiment),
                           color=colors[index % len(colors)])
                model = sm.OLS(epsilon_dict[at_gc_pct][experiment], sm.add_constant(nbr_gc))
                results = model.fit()
                b, a = results.params[0:2]
                idf = np.linspace(min(nbr_gc), max(nbr_gc), 30)
                ax.plot(idf, a * idf + b, 'r-',
                        label=r"$y={0:.3g}x {3} {1:.3g}$ ($r^2={2:.3g})$".format(
                            float(a), abs(float(b)), results.rsquared, "+" if float(b) > 0 else "-"),
                        color=colors[index % len(colors)])

            ax.legend(fontsize=fontsize)
            ax.set_title("$\lambda=$" + str(at_gc_pct), fontsize=fontsize)
        plt.show()
