/* T. Latrille.
Hyphy inference for an "experimental" dataset using GTR mutation matrix
*/

global w=1.0; w:> 0; t; t:>0; // note that we use "global t" (instead of locally estimated for each branch) since all branch lengths are the same in the simulation tree.
{% for nuc_vars in nuc_vars_dict.values() %}
{{ nuc_vars }}{% endfor %}{% for rate_vars in rate_vars_dict.values() %}
{{ rate_vars }}{% endfor %}

LIKELIHOOD_FUNCTION_OUTPUT = 1;
RANDOM_STARTING_PERTURBATIONS = 1;
OPTIMIZATION_PRECSION = 0.00000001;

{% for rate_param, freq_to_matrix_dict in matrices_dict.items() %}
{% for freq_param, matrix in freq_to_matrix_dict.items() %}
Matrix_r{{ rate_param }}_f{{ freq_param }}={{ matrix }};
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

{% for rate_param, freq_to_matrix_dict in matrices_dict.items() %}
{% for freq_param, matrix in freq_to_matrix_dict.items() %}
////////////// r{{ rate_param }}_f{{ freq_param }} //////////////
global w=1.0; w:> 0; t; t:>0;
{{ nuc_vars_dict[freq_param] }}
{{ rate_vars_dict[rate_param] }}
Model Modelr{{ rate_param }}_f{{ freq_param }} = (Matrix_r{{ rate_param }}_f{{ freq_param }}, Freq_{{ freq_param }}, 0);
UseModel (USE_NO_MODEL);
UseModel(Modelr{{ rate_param }}_f{{ freq_param }});
Tree    Tree01 = {{ tree }};
LikelihoodFunction  LikModelr{{ rate_param }}_f{{ freq_param }} = (filt_data, Tree01);
Optimize (paramValues, LikModelr{{ rate_param }}_f{{ freq_param }});
fprintf ("{{ name }}_r{{ rate_param }}_f{{ freq_param }}_hyout.txt", LikModelr{{ rate_param }}_f{{ freq_param }});
{% endfor %}
{% endfor %}