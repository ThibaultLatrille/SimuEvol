/* T. Latrille.
Hyphy inference for an "experimental" dataset using GTR mutation matrix
*/
{% for var in vars_all.values() %}{{ var }} {% endfor %}

LIKELIHOOD_FUNCTION_OUTPUT = 1;
RANDOM_STARTING_PERTURBATIONS = 1;
OPTIMIZATION_PRECSION = 0.00000001;

{% for rate_param, matrices_dict_l1 in matrices_nested_dict.items() %}
{% for freq_param, matrices_dict_l2 in matrices_dict_l1.items() %}
{% for omega_param, matrix in matrices_dict_l2.items() %}
Matrix_r{{ rate_param  }}_f{{ freq_param }}_w{{ omega_param }}={{ matrix }};
{% endfor %}
{% endfor %}
{% endfor %}

/* Read in the data */
DataSet	raw_data = ReadDataFile("{{ fasta_infile }}");

/* Filter the data to find and remove any stop codons*/
DataSetFilter   filt_data = CreateFilter(raw_data,3,"", "","TAA,TAG,TGA");


/* Set up frequencies. Note that these were all hard-coded in when the file was created via the script SimuEvol/scripts/prefs_to_freqs.py */

{% for freq_param, freqs_dict_l1 in freqs_nested_dict.items() %}
{% for omega_param, freqs in freqs_dict_l1.items() %}
Freq_f{{ freq_param }}_w{{ omega_param }}={{ freqs }};
{% endfor %}
{% endfor %}

/* Optimize likelihoods for each frequency specification */

{% for rate_param, matrices_dict_l1 in matrices_nested_dict.items() %}
{% for freq_param, matrices_dict_l2 in matrices_dict_l1.items() %}
{% for omega_param, matrix in matrices_dict_l2.items() %}
////////////// w{{ omega_param }}_r{{ rate_param }}_f{{ freq_param }} //////////////
{% for var in vars_nested_dict[rate_param][freq_param][omega_param].values() %}{{ var }} {% endfor %}

Model Model_r{{ rate_param  }}_f{{ freq_param }}_w{{ omega_param }} = (Matrix_r{{ rate_param  }}_f{{ freq_param }}_w{{ omega_param }}, Freq_f{{ freq_param }}_w{{ omega_param }}, 0);
UseModel (USE_NO_MODEL);
UseModel(Model_r{{ rate_param  }}_f{{ freq_param }}_w{{ omega_param }});
Tree    Tree_r{{ rate_param  }}_f{{ freq_param }}_w{{ omega_param }} = {{ tree }};

branchNames 	= BranchName (Tree_r{{ rate_param  }}_f{{ freq_param }}_w{{ omega_param }}, -1);
branchLengths	= BranchLength (Tree_r{{ rate_param  }}_f{{ freq_param }}_w{{ omega_param }}, -1);

for (k = 0; k < Columns(branchNames)-1; k=k+1)
{
	ExecuteCommands("Tree_r{{ rate_param  }}_f{{ freq_param }}_w{{ omega_param }}." + branchNames[k] + ".t:=" + branchLengths[k] + ";");
}

LikelihoodFunction  LikModel_r{{ rate_param  }}_f{{ freq_param }}_w{{ omega_param }} = (filt_data, Tree_r{{ rate_param  }}_f{{ freq_param }}_w{{ omega_param }});
Optimize (paramValues, LikModel_r{{ rate_param  }}_f{{ freq_param }}_w{{ omega_param }});
fprintf (stdout, LikModel_r{{ rate_param  }}_f{{ freq_param }}_w{{ omega_param }});
fprintf ("{{ name }}_r{{ rate_param  }}_f{{ freq_param }}_w{{ omega_param }}_hyout.txt", LikModel_r{{ rate_param  }}_f{{ freq_param }}_w{{ omega_param }});
{% endfor %}{% endfor %}{% endfor %}
