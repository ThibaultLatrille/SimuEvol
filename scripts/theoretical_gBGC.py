# GLOBAL IMPORTS
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage

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
weak_nucleotides = "AT"
strong_nucleotides = "CG"
assert len(codon_table.keys()) == 64, "There is not 64 codons in the codon table"

codons = list([k for k, v in codon_table.items() if v != 'X'])
assert len(codons) == 61, "There is not 3 stop codons in the codon table"
assert not [n for n in "".join(codons) if (n not in nucleotides)], "There is a codon with an unrecognized nucleotide"
nbr_weak = np.array([sum([c.count(weak) for weak in weak_nucleotides]) for c in codons])
nbr_strong = 3 - nbr_weak

amino_acids_set = set(codon_table.values())
amino_acids_set.remove('X')

aa_char_to_int = {v: k for k, v in enumerate(sorted(amino_acids_set))}
amino_acids = "".join(amino_acids_set)
assert len(amino_acids) == 20, "There is not 20 amino-acids in the codon table"


def generate_codon_to_aa():
    codon_to_aa_list = [None] * len(codons)
    for codon_index, codon in enumerate(codons):
        codon_to_aa_list[codon_index] = aa_char_to_int[codon_table[codon]]
    return codon_to_aa_list


def generate_neighbors(synonymous=True):
    neighbor_dict = dict()
    for x, codon_origin in enumerate(codons):
        neighbor_dict[x] = []
        for y, codon_target in enumerate(codons):
            mutations = [(n_from, n_to) for n_from, n_to in zip(codon_origin, codon_target) if n_from != n_to]
            if len(mutations) == 1:
                if synonymous == (codon_table[codon_origin] == codon_table[codon_target]):
                    nuc_origin, nuc_target = mutations[0]
                    neighbor_dict[x].append((y, nucleotides.index(nuc_origin), nucleotides.index(nuc_target)))
    return neighbor_dict


codon_to_aa = generate_codon_to_aa()
non_syn_neighbors = generate_neighbors(synonymous=False)
syn_neighbors = generate_neighbors(synonymous=True)


def simpson_diversity(frequencies):
    return 1.0 / np.sum(frequencies * frequencies)


def shanon_diversity(frequencies):
    return np.exp(-np.sum(frequencies * np.log(frequencies)))


def hill_diversity(frequencies, q):
    if q == 1:
        return shanon_diversity(frequencies)
    else:
        return np.power(np.sum(np.power(frequencies, q)), 1.0 / (1.0 - q))


def compute_diversity(frequencies_seq):
    n = frequencies_seq.shape[0]
    diversity_seq = np.zeros(n)
    for site in range(n):
        diversity_seq[site] = simpson_diversity(frequencies_seq[site, :])

    diversity_mean = simpson_diversity(np.mean(frequencies_seq, axis=0))
    return np.mean(diversity_seq), diversity_mean


def compute_diversity_aa(codon_frequencies_seq):
    n = codon_frequencies_seq.shape[0]
    aa_frequencies_seq = np.zeros((n, len(amino_acids)))
    for site in range(n):
        for codon_index, freq in enumerate(codon_frequencies_seq[site, :]):
            aa_frequencies_seq[site, codon_to_aa[codon_index]] += freq

    return compute_diversity(aa_frequencies_seq)


def compute_lambda_obs(codon_frequencies_seq):
    n = codon_frequencies_seq.shape[0]
    at_pct_seq = np.zeros(n)
    for site in range(n):
        at_pct_seq[site] = np.sum(codon_frequencies_seq[site, :] * nbr_weak) / 3

    at_pct = np.mean(at_pct_seq, axis=0)
    return at_pct / (1 - at_pct)


def compute_omega(codon_frequencies_seq, aa_fitness_seq, mutation_matrix, gBGC_param):
    sub_flow_non_syn, mut_flow_non_syn = 0, 0
    sub_flow_syn, mut_flow_syn = 0, 0

    n = codon_frequencies_seq.shape[0]
    for site in range(n):
        aa_fitness_site = aa_fitness_seq[site, :]
        codon_frequencies_site = codon_frequencies_seq[site, :]
        for x in range(len(codons)):
            sub_flow_non_syn_tmp, mut_flow_non_syn_tmp = 0, 0
            sub_flow_syn_tmp, mut_flow_syn_tmp = 0, 0

            for y, a, b in non_syn_neighbors[x]:
                sel_coef = aa_fitness_site[codon_to_aa[y]] - aa_fitness_site[codon_to_aa[x]]
                if nucleotides[b] in strong_nucleotides:
                    sel_coef += gBGC_param
                if abs(sel_coef) < 1e-10:
                    p_fix = 1.0
                else:
                    p_fix = sel_coef / (1. - np.exp(-sel_coef))
                sub_flow_non_syn_tmp += mutation_matrix[a][b] * p_fix
                mut_flow_non_syn_tmp += mutation_matrix[a][b]

            for y, a, b in syn_neighbors[x]:
                if (nucleotides[b] in strong_nucleotides) and abs(gBGC_param) > 1e-10:
                    p_fix = gBGC_param / (1. - np.exp(-gBGC_param))
                else:
                    p_fix = 1.0

                sub_flow_syn_tmp += mutation_matrix[a][b] * p_fix
                mut_flow_syn_tmp += mutation_matrix[a][b]

            sub_flow_non_syn += codon_frequencies_site[x] * sub_flow_non_syn_tmp
            mut_flow_non_syn += codon_frequencies_site[x] * mut_flow_non_syn_tmp
            sub_flow_syn += codon_frequencies_site[x] * sub_flow_syn_tmp
            mut_flow_syn += codon_frequencies_site[x] * mut_flow_syn_tmp

    dn = sub_flow_non_syn / mut_flow_non_syn
    ds = sub_flow_syn / mut_flow_syn
    return dn / ds


def generate_frequencies(alpha_param, gBGC_param, lambda_param, n):
    mutation_matrix = np.ones((len(nucleotides), len(nucleotides)))
    for a in range(len(nucleotides)):
        for b in range(len(nucleotides)):
            if nucleotides[b] in weak_nucleotides:
                mutation_matrix[a, b] = lambda_param

    aa_propensities_seq = np.random.dirichlet(alpha_param * np.ones(len(amino_acids)), n)
    aa_fitness_seq = np.log(aa_propensities_seq)
    codon_frequencies_seq = np.zeros((n, len(codons)))
    mut_frequencies = np.power(lambda_param, nbr_weak)
    mut_frequencies *= np.exp(gBGC_param * nbr_strong)

    for site in range(n):
        aa_propensities_site = aa_propensities_seq[site, :]
        codon_propensities_site = np.array([aa_propensities_site[codon_to_aa[i]] for i in range(len(codons))])
        codon_frequencies_site = mut_frequencies * codon_propensities_site
        codon_frequencies_seq[site, :] = codon_frequencies_site / np.sum(codon_frequencies_site)

    omega = compute_omega(codon_frequencies_seq, aa_fitness_seq, mutation_matrix, gBGC_param)
    lambda_obs = compute_lambda_obs(codon_frequencies_seq)
    diversity_codon, diversity_mean_codon = compute_diversity(codon_frequencies_seq)
    diversity_aa, diversity_mean_aa = compute_diversity_aa(codon_frequencies_seq)

    return omega, lambda_obs / lambda_param, \
           diversity_codon, diversity_mean_codon, diversity_aa, diversity_mean_aa


def plot_heatmap(nbr_steps, n, alpha):
    zoom = 3
    x_axis = np.logspace(-1, 1, nbr_steps)
    y_axis = np.linspace(0, 5, nbr_steps)
    zommed_x_axis = np.logspace(-1, 1, nbr_steps * zoom)
    zommed_y_axis = np.linspace(0, 5, nbr_steps * zoom)
    y, x = np.meshgrid(zommed_y_axis, zommed_x_axis)

    heatmap_summary_stats = list()
    stat_names = ['$\\omega$',
                  '$\\lambda_{obs} / \\lambda$',
                  '$D_{\\mathrm{C}}^{\\mathrm{site}}$', '$D_{\\mathrm{C}}^{\\mathrm{seq}}$',
                  '$D_{\\mathrm{A}}^{\\mathrm{site}}$', '$D_{\\mathrm{A}}^{\\mathrm{seq}}$']
    stat_title = ['omega',
                  'at_over_gc',
                  'diversity_site_codon', 'diversity_seq_codon',
                  'diversity_site_aa', 'diversity_seq_aa']
    for _ in range(len(stat_names)):
        heatmap_summary_stats.append(np.zeros((nbr_steps, nbr_steps)))

    for x_i, lambda_param in enumerate(x_axis):
        for y_i, gBGC_param in enumerate(y_axis):
            stat_list = generate_frequencies(alpha, gBGC_param, lambda_param, n)
            for stat_id, stat in enumerate(stat_list):
                heatmap_summary_stats[stat_id][x_i][y_i] = stat
            print("{0:.4g}% computed".format(100 * (x_i * nbr_steps + (y_i + 1)) / (nbr_steps ** 2)))

    for index, z in enumerate(heatmap_summary_stats):
        zommed_z = scipy.ndimage.zoom(z, zoom)
        my_dpi = 196
        fig = plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)

        contourf = plt.contourf(x, y, zommed_z, 16, cmap='RdBu',
                                vmin=np.nanmin(zommed_z),
                                vmax=np.nanmax(zommed_z))
        contour = plt.contour(contourf,
                              levels=contourf.levels[::2],
                              linewidths=(1,), colors='black')
        plt.axis([x.min(), x.max(), y.min(), y.max()])
        cbar = plt.colorbar(contourf)
        cbar.ax.set_ylabel(stat_names[index], fontsize=16)
        # Add the contour line levels to the colorbar
        cbar.add_lines(contour)
        plt.clabel(contour, fmt='%2.1f', colors='black', fontsize=8)
        plt.xscale("log")
        plt.xlabel('$\\lambda$', fontsize=16)
        plt.ylabel('$B$', fontsize=16)
        plt.tight_layout()
        plt.title('$n = {0}$'.format(n))
        # set the limits of the plot to the limits of the data

        # plt.savefig("../figures/gBGC_{0}.svg".format(stat_title[index]), format="svg")
        plt.savefig("../figures/gBGC_{0}.png".format(stat_title[index]), format="png")
        plt.clf()
        plt.close('all')


plot_heatmap(30, 6000, 0.2)
