# GLOBAL IMPORTS
from analysis import *
from interactive_plot import *
import statsmodels.api as sm
from math import ceil

current_dir = "/home/thibault/SimuEvol/simulated"
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

            nuc_freqs, nbr_sites, nbr_species = extract_nuc_pct(simu_evol_path + ".fasta")
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
                hyphy_dico["n"] = nbr_sites
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
