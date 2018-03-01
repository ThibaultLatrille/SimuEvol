/* T. Latrille.
Hyphy inference for an "experimental" dataset using GTR mutation matrix
*/

t; t:>0;
{% for omega_vars in omega_vars_dict.values() %}
{{ omega_vars }}{% endfor %}{% for nuc_vars in nuc_vars_dict.values() %}
{{ nuc_vars }}{% endfor %}{% for rate_vars in rate_vars_dict.values() %}
{{ rate_vars }}{% endfor %}


LIKELIHOOD_FUNCTION_OUTPUT = 1;
RANDOM_STARTING_PERTURBATIONS = 1;
OPTIMIZATION_PRECSION = 0.00000001;

{% for rate_param, matrices_dict_l1 in matrices_dict.items() %}
{% for freq_param, matrices_dict_l2 in matrices_dict_l1.items() %}
{% for omega_param, matrix in matrices_dict_l2.items() %}
Matrix_w{{ omega_param }}_r{{ rate_param }}_f{{ freq_param }}={{ matrix }};
{% endfor %}
{% endfor %}
{% endfor %}

/* Read in the data */
DataSet	raw_data = ReadDataFile("{{ fasta_infile }}");

/* Filter the data to find and remove any stop codons*/
DataSetFilter   filt_data = CreateFilter(raw_data,3,"", "","TAA,TAG,TGA");


/* Set up frequencies. Note that these were all hard-coded in when the file was created via the script SimuEvol/scripts/prefs_to_freqs.py */

{% for freq_param, freqs in codon_freqs_dict.items() %}
Freq_{{ freq_param }}={{ freqs }};
{% endfor %}

/* Optimize likelihoods for each frequency specification */

{% for rate_param, matrices_dict_l1 in matrices_dict.items() %}
{% for freq_param, matrices_dict_l2 in matrices_dict_l1.items() %}
{% for omega_param, matrix in matrices_dict_l2.items() %}
////////////// w{{ omega_param }}_r{{ rate_param }}_f{{ freq_param }} //////////////
t; t:>0;
{{ omega_vars_dict[omega_param] }}
{{ nuc_vars_dict[freq_param] }}
{{ rate_vars_dict[rate_param] }}
Model Model_w{{ omega_param }}_r{{ rate_param }}_f{{ freq_param }} = (Matrix_w{{ omega_param }}_r{{ rate_param }}_f{{ freq_param }}, Freq_{{ freq_param }}, 0);
UseModel (USE_NO_MODEL);
UseModel(Model_w{{ omega_param }}_r{{ rate_param }}_f{{ freq_param }});
Tree    Tree01 = {{ tree }};
LikelihoodFunction  LikModel_w{{ omega_param }}_r{{ rate_param }}_f{{ freq_param }} = (filt_data, Tree01);
Optimize (paramValues, LikModel_w{{ omega_param }}_r{{ rate_param }}_f{{ freq_param }});
fprintf ("{{ name }}_w{{ omega_param }}_r{{ rate_param }}_f{{ freq_param }}_hyout.txt", LikModel_w{{ omega_param }}_r{{ rate_param }}_f{{ freq_param }});
{% endfor %}{% endfor %}{% endfor %}
