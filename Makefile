
# ==============================================================================================================
#  COMPILATION
# ==============================================================================================================
.PHONY: all # Requires: cmake 3.1.0 or better
all: build
	@cd build ; make --no-print-directory -j8

build: CMakeLists.txt
	@rm -rf build
	@mkdir build
	@cd build ; cmake ..

.PHONY: coverage
coverage:
	@rm -rf build
	@mkdir build
	@cd build ; cmake -DCOVERAGE_MODE=ON ..
	@make --no-print-directory test

.PHONY: debug
debug:
	@rm -rf build
	@mkdir build
	@cd build ; cmake -DDEBUG_MODE=ON ..
	@make --no-print-directory

.PHONY: release
release:
	@rm -rf build
	@mkdir build
	@cd build ; cmake ..
	@make --no-print-directory

.PHONY: tiny
tiny:
	@rm -rf build
	@mkdir build
	@cd build ; cmake -DTINY=ON ..
	@make --no-print-directory

.PHONY: clean
clean:
	@rm -rf build
	@rm -rf test

# ==============================================================================================================
#  CODE QUALITY
# ==============================================================================================================
.PHONY: format # Requires: clang-format
format:
	@clang-format -i `find src/ -name *.*pp`

# ==============================================================================================================
#  TESTING
# ==============================================================================================================

.PHONY: gBGC
gBGC: tiny
	@rm -rf gBGC
	@mkdir gBGC
	@echo "\n\e[35m\e[1m== SimuProfile run ===============================================================\e[0m"
	build/SimuProfile --preferences data/preferences/gal4.prefs --newick data/trees/mammal_subtree.tre.annotated --nuc_matrix data/matrices/nucleotide_GTR.tsv --precision_matrix data/matrices/precision4x4.tsv --exon_size 25 --output gBGC/SimuProfile_gal4 --gBGC 1.0 --bias_gBGC 1.0 --branch_wise_correlation True


.PHONY: test
test: debug
	@rm -rf test
	@mkdir test
	@echo "\n\e[35m\e[1m== ToyGeo run =========================================================\e[0m"
	build/ToyGeo --output test/ToyGeo_run
	@echo "\n\e[35m\e[1m== ToyStab run =========================================================\e[0m"
	build/ToyStab --output test/ToyStab_run
	@echo "\n\e[35m\e[1m== SimuProfile run ===============================================================\e[0m"
	build/SimuProfile --preferences data/preferences/gal4.prefs --newick data/trees/mammal_subtree.tre.annotated --nuc_matrix data/matrices/nucleotide_GTR.tsv --precision_matrix data/matrices/precision4x4.tsv --exon_size 25 --output test/SimuProfile_gal4
	@echo "\n\e[35m\e[1m== SimuDfe run ===============================================================\e[0m"
	build/SimuDfe --newick data/trees/mammal_subtree.tre.annotated --nuc_matrix data/matrices/nucleotide_GTR.tsv --precision_matrix data/matrices/precision4x4.tsv --exon_size 50 --nbr_exons 2 --population_size 1e4 --gamma_mean 1e-5 --output test/SimuDfe_gal4
	@echo "\n\e[35m\e[1m== SimuGeo run ===============================================================\e[0m"
	build/SimuGeo --complexity 3 --newick data/trees/mammal_subtree.tre.annotated --nuc_matrix data/matrices/nucleotide_GTR.tsv --precision_matrix data/matrices/precision4x4.tsv --exon_size 200 --nbr_exons 5 --output test/SimuGeo
	@echo "\n\e[35m\e[1m== SimuFold run =========================================================\e[0m"
	build/SimuFold --pdb_folder data/pdbfiles --newick data/trees/mammal_subtree.tre.annotated --nuc_matrix data/matrices/nucleotide_GTR.tsv --precision_matrix data/matrices/precision4x4.tsv --exon_size 300 --nbr_exons 2 --population_size 1e5 --mutation_rate_per_generation 1e-9 --output test/SimuFold
	@echo "\n\e[35m\e[1m== SimuStab run =========================================================\e[0m"
	build/SimuStab --newick data/trees/mammal_subtree.tre.annotated --nuc_matrix data/matrices/nucleotide_GTR.tsv --precision_matrix data/matrices/precision4x4.tsv --exon_size 300 --nbr_exons 2 --population_size 1e5 --mutation_rate_per_generation 1e-9 --output test/SimuStab
	@echo "\n\e[35m\e[1m== SimuPoly run ==============================================================\e[0m"
	build/PolyProfile --preferences data/preferences/gal4.prefs --newick data/trees/mammal_subtree.tre.annotated --nuc_matrix data/matrices/nucleotide_GTR.tsv --precision_matrix data/matrices/precision4x4.tsv --exon_size 25 --noise_sigma 0.1 --output test/SimuPoly_gal4
	@echo "\n\e[35m\e[1m== PolyDfe run ==============================================================\e[0m"
	build/PolyDfe --newick data/trees/mammal_subtree.tre.annotated --nuc_matrix data/matrices/nucleotide_GTR.tsv --precision_matrix data/matrices/precision4x4.tsv --exon_size 50 --nbr_exons 2 --population_size 200 --gamma_mean 1e-3 --noise_sigma 0.1 --output test/SimuDfe_gal4
	@echo "\n\e[35m\e[1m== SimuRelax run =============================================================\e[0m"
	build/Relaxation --output test/SimuRelax_run
