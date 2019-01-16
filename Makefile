# ==============================================================================================================
#  COMPILATION
# ==============================================================================================================
.PHONY: all # Requires: cmake 3.1.0 or better
all: _build
	@cd bin ; make --no-print-directory -j8

_build: CMakeLists.txt
	@rm -rf bin
	@mkdir bin
	@cd bin ; cmake ..

.PHONY: coverage
coverage:
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake -DCOVERAGE_MODE=ON ..
	@make --no-print-directory test
	# find _build_coverage -type f -name "*.gcno" -exec mv -t src/ {} +
	# find _build_coverage -type f -name "*.gcda" -exec mv -t src/ {} +

.PHONY: clean
clean:
	@rm -rf bin
	@rm -rf _build
	@rm -rf _test

# ==============================================================================================================
#  CODE QUALITY
# ==============================================================================================================
.PHONY: format # Requires: clang-format
format:
	@clang-format -i `find -name *.*pp`

# ==============================================================================================================
#  TESTING
# ==============================================================================================================

.PHONY: test
test: all
	@rm -rf _test
	@mkdir _test
	@echo "\n\e[35m\e[1m== SimuRelax run ===============================================================\e[0m"
	bin/SimuRelax --output=_test/SimuRelax_run
	@echo "\n\e[35m\e[1m== SimuEvol run ===============================================================\e[0m"
	bin/SimuEvol --preferences=data/preferences/gal4.txt --newick=data/trees/gal4.newick --output=_test/SimuEvol_gal4
	@echo "\n\e[35m\e[1m== SimuPoly run ===============================================================\e[0m"
	bin/SimuPoly --preferences=data/preferences/gal4.txt --newick=data/trees/gal4.newick --output=_test/SimuPoly_gal4