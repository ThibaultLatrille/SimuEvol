#pragma once

#include <cassert>
#include <cmath>
#include <fstream>
#include "io.hpp"
#include "tclap/CmdLine.h"

class OutputArgParse {
  protected:
    TCLAP::CmdLine& cmd;

  public:
    explicit OutputArgParse(TCLAP::CmdLine& cmd) : cmd{cmd} {}

    TCLAP::ValueArg<std::string> output_path{"o", "output", "output path", true, "", "string", cmd};
    TCLAP::ValueArg<u_long> seed{
        "", "seed", "Random number generation seed", false, 0, "u_long", cmd};
};

class SimuArgParse : public OutputArgParse {
  public:
    explicit SimuArgParse(TCLAP::CmdLine& cmd) : OutputArgParse(cmd) {}

    TCLAP::ValueArg<std::string> precision_path{
        "c", "precision_matrix", "input precision matrix path", false, "", "string", cmd};
    TCLAP::ValueArg<double> mutation_rate_per_generation{"m", "mutation_rate_per_generation",
        "Mutation rate (at the root)", false, 1e-8, "double", cmd};
    TCLAP::ValueArg<double> root_age{
        "a", "root_age", "Age of the root", false, 50e6, "double", cmd};
    TCLAP::ValueArg<double> generation_time{"g", "generation_time",
        "Number of year between generations (at the root)", false, 40, "double", cmd};
    TCLAP::ValueArg<std::string> nuc_matrix_path{
        "q", "nuc_matrix", "input nucleotide matrix preferences path", false, "", "string", cmd};
    TCLAP::ValueArg<std::string> newick_path{
        "t", "newick", "input newick tree path", true, "", "string", cmd};
    TCLAP::ValueArg<u_long> exons{"s", "exon_size",
        "Number of codon sites per exon (default 0 means the size of the fitness profiles "
        "provided, thus assuming complete linkage between sites)",
        false, 0, "u_long", cmd};
    TCLAP::SwitchArg branch_wise_correlation{"", "branch_wise_correlation",
        "The Correlated parameters are determined branch-wise", cmd, false};

    TCLAP::SwitchArg fix_pop_size{"", "fix_pop_size", "Population size is fixed", cmd, false};
    TCLAP::SwitchArg fix_mut_rate{"", "fix_mut_rate", "Mutation rate is fixed", cmd, false};
    TCLAP::SwitchArg fix_gen_time{"", "fix_gen_time", "Generation time is fixed", cmd, false};

    void add_to_trace(Trace& trace) {
        std::cout << "Random generator seed: " << seed.getValue() << std::endl;
        assert(mutation_rate_per_generation.getValue() > 0.0);
        assert(root_age.getValue() > 0.0);
        assert(generation_time.getValue() > 0.0);
        assert(generation_time.getValue() < root_age.getValue());
        trace.add("fix_pop_size", fix_pop_size.getValue());
        trace.add("fix_mut_rate", fix_mut_rate.getValue());
        trace.add("fix_gen_time", fix_gen_time.getValue());
        trace.add("seed", seed.getValue());
        trace.add("output_path", output_path.getValue());
        trace.add("tree_path", newick_path.getValue());
        trace.add("generation_time", generation_time.getValue());
        trace.add("exon_size", exons.getValue());
        trace.add("nucleotide_matrix_path", output_path.getValue());
        trace.add("mutation_rate_per_generation", mutation_rate_per_generation.getValue());
    }
};