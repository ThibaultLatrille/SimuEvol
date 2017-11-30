"How to calculate the non-synonymous to synonymous substitution rate ratio of protein-coding genes under the Fisher-Wright mutation-selection framework"
Mario dos Reis (2015) Biology Letters, 20141031, DOI: 10.1098/rsbl.2014.1031

The alignment (pb2.phyl, rbcl.phyl) and tree (pb2.tree, rbc.tree) files can be analysed with the program swMutSel (https://github.com/tamuri/swmutsel) to estimate the fitnesses of amino acids at codon sites and other evolutionary parameters for the PB2 and RBCL proteins. The estimated fitness values (e.g. see file rbcl_dir001_MLE.txt) can then be used to calculate the relative non-synonymous (and synonymous) substitution rates as described in the paper (Eq. 3.1, 3.2). R scripts to perform the calculations are available from the author (mariodosreis@gmail.com).

A detailed description of the datasets and of the statistical method to estimate the fitnesses of amino acids at codon sites can be found in

Tamuri AU, dos Reis M and Goldstein RA. (2012)  Estimating the distribution of selection coefficients from phylogenetic data using sitewise mutation-selection models. Genetics, 190: 1101-1115. DOI: 10.1534/genetics.111.136432

and

Tamuri AU, Goldman N and dos Reis M. (2014) A penalized-likelihood method to estimate the distribution of selection coefficients from phylogenetic data. Genetics, 197: 257-271. DOI: 10.1534/genetics.114.162263


FILES:

pb2.phyl: Alignment (in phylip format) of 401 sequences of the pb2 gene of influenza viruses isolated from human and avian hosts. The alignment is 2,277 nucleotides (759 codons) long. The first two characters of the sequence name indicate the host (Av: Avian, Hu: Human), and the characters after the underscore '_' are the GenBank accession number.

pb2.tree: Phylogenetic tree (in Newick format) of 401 sequences of the PB2 gene of influenza viruses. The branch lengths are in number of substitutions per codon site estimated under the FMutSel0 model with the program CODEML.

rbcl.phyl: Alignment (in phylip format) of 3,490 sequences of the rbcL chloroplast gene of Monocots (a group of flowering plants). The alignment is 1,371 nucleotides (457 codons) long. In the sequence name, the characters after the underscore '_' are the GenBank accession numbers.

rbcl.tree: Phylogenetic tree (in Newick format) of 3,490 sequences of the choloplast gene of Monocots. The branch lengths are in number of substitutions per codon site estimated under the FMutSel0 model with the program CODEML.

pb2_dir001_MLE.txt: Fitnesses and evolutionary parameters for the 401 pb2 sequences of influenza estimated by penalised likelihood with the program swMutSel. A Dirichlet-based penalty with alpha=0.01 was used (see Tamuri et al. 2014).

pb2_fitness_av.txt: Fitness estiamtes for 25 adaptive locations for viruses evolving in avian hosts. The estimates were obtained without penalty (see Tamuri et al. 2012).

pb2_fitness_hu.txt: Fitnesses estimates for 25 adaptive locations for viruses evolving in human hosts. The estimates were obtained without penalty (see Tamuri et al. 2012).


ACKNOWLEDGEMENTS:

Thanks to Alexandros Stamatakis and Guido Grimm for providing the rbcL sequence alignment.