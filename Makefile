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

.PHONY: test
test: all
	@rm -rf test
	@mkdir test
	@echo "\n\e[35m\e[1m== SimuDiv run ===============================================================\e[0m"
	build/SimuDiv --preferences data/preferences/gal4.prefs --newick data/trees/mammal_subtree.tre --nuc_matrix data/matrices/nucleotide_GTR.tsv --correlation_matrix data/matrices/correlation3x3_zero.tsv --output test/SimuDiv_gal4
	@echo "\n\e[35m\e[1m== SimuPoly run ===============================================================\e[0m"
	build/SimuPoly --preferences data/preferences/gal4.prefs --newick data/trees/mammal_subtree.tre --nuc_matrix data/matrices/nucleotide_GTR.tsv --correlation_matrix data/matrices/correlation3x3_zero.tsv --output test/SimuPoly_gal4
	@echo "\n\e[35m\e[1m== SimuRelax run ===============================================================\e[0m"
	build/SimuRelax --output test/SimuRelax_run
