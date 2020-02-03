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

    TCLAP::ValueArg<double> mutation_rate_per_generation{"", "mutation_rate_per_generation",
        "Mutation rate (at the root)", false, 1e-8, "double", cmd};
    TCLAP::ValueArg<double> generation_time{"", "generation_time",
        "Number of year between generations (at the root)", false, 40, "double", cmd};
    TCLAP::ValueArg<std::string> nuc_matrix_path{
        "", "nuc_matrix", "input nucleotide matrix preferences path", false, "", "string", cmd};
    TCLAP::ValueArg<u_long> exons{"", "exon_size",
        "Number of codon sites per exon (default 0 means the size of the fitness profiles "
        "provided, thus assuming complete linkage between sites)",
        false, 0, "u_long", cmd};

    virtual void add_to_trace(Trace& trace) {
        std::cout << "Random generator seed: " << seed.getValue() << std::endl;
        assert(mutation_rate_per_generation.getValue() > 0.0);
        assert(generation_time.getValue() > 0.0);
        trace.add("seed", seed.getValue());
        trace.add("output_path", output_path.getValue());
        trace.add("generation_time", generation_time.getValue());
        trace.add("exon_size", exons.getValue());
        trace.add("nucleotide_matrix_path", output_path.getValue());
        trace.add("mutation_rate_per_generation", mutation_rate_per_generation.getValue());
    }
};

class SimuSubArgParse : public SimuArgParse {
  public:
    explicit SimuSubArgParse(TCLAP::CmdLine& cmd) : SimuArgParse(cmd) {}
    TCLAP::ValueArg<u_long> nbr_grid_step{"", "nbr_grid_step",
        "Number of intervals in which to discretize the brownian motion", false, 100, "u_long",
        cmd};
    TCLAP::ValueArg<double> pop_size{
        "", "population_size", "Effective population size", false, 1.0, "double", cmd};

    void add_to_trace(Trace& trace) override {
        SimuArgParse::add_to_trace(trace);
        trace.add("pop_size", pop_size.getValue());
        trace.add("nbr_grid_step", nbr_grid_step.getValue());
    }
};

class SimuPolyArgParse : public SimuArgParse {
  public:
    explicit SimuPolyArgParse(TCLAP::CmdLine& cmd) : SimuArgParse(cmd) {}
    TCLAP::ValueArg<u_long> sample_size{
        "", "sample_size", "Sample size (at the leaves)", false, 20, "u_long", cmd};
    TCLAP::ValueArg<double> noise_sigma{"", "noise_sigma",
        "The Ornstein–Uhlenbeck sigma (0<sigma) applied to Ne at each generation", false, 0.0,
        "double", cmd};
    TCLAP::ValueArg<double> noise_theta{"", "noise_theta",
        "The Ornstein–Uhlenbeck theta (0<=theta<1) applied to Ne at each generation", false, 0.9,
        "double", cmd};
    TCLAP::ValueArg<u_long> pop_size{
        "", "population_size", "Effective population size", false, 5000, "u_long", cmd};

    void add_to_trace(Trace& trace) override {
        SimuArgParse::add_to_trace(trace);
        trace.add("sample_size", sample_size.getValue());
        trace.add("noise_sigma", noise_sigma.getValue());
        trace.add("noise_theta", noise_theta.getValue());
        trace.add("population_size", pop_size.getValue());
    }
};

class TreeArgParse {
  protected:
    TCLAP::CmdLine& cmd;

  public:
    explicit TreeArgParse(TCLAP::CmdLine& cmd) : cmd{cmd} {}
    TCLAP::ValueArg<double> root_age{"", "root_age", "Age of the root", false, 50e6, "double", cmd};
    TCLAP::ValueArg<std::string> newick_path{
        "", "newick", "input newick tree path", false, "", "string", cmd};
    TCLAP::ValueArg<double> nbr_branches{"", "nbr_branches",
        "Nbr of sucessive branches if no tree is inputed", false, 2, "unsigned", cmd};
    TCLAP::ValueArg<std::string> precision_path{
        "", "precision_matrix", "input precision matrix path", false, "", "string", cmd};
    TCLAP::SwitchArg branch_wise_correlation{"", "branch_wise_correlation",
        "The Correlated parameters are determined branch-wise", cmd, false};
    TCLAP::SwitchArg fix_pop_size{
        "", "fix_pop_size", "Log-Brownian on population size is null", cmd, false};
    TCLAP::SwitchArg fix_mut_rate{
        "", "fix_mut_rate", "Log-Brownian on mutation rate is null", cmd, false};
    TCLAP::SwitchArg fix_gen_time{
        "", "fix_gen_time", "Log-Brownian on generation time is null", cmd, false};

    TCLAP::ValueArg<double> bias_pop_size{
        "", "bias_pop_size", "Log-Bias on population size", false, 0., "double", cmd};
    TCLAP::ValueArg<double> bias_mut_rate{
        "", "bias_mut_rate", "Log-Bias on mutation rate", false, 0., "double", cmd};
    TCLAP::ValueArg<double> bias_gen_time{
        "", "bias_gen_time", "Log-Bias on generation time", false, 0., "double", cmd};

    void add_to_trace(Trace& trace) {
        assert(root_age.getValue() > 0.0);
        trace.add("root_age", root_age.getValue());
        trace.add("fix_pop_size", fix_pop_size.getValue());
        trace.add("fix_mut_rate", fix_mut_rate.getValue());
        trace.add("fix_gen_time", fix_gen_time.getValue());
        trace.add("bias_pop_size", bias_pop_size.getValue());
        trace.add("bias_mut_rate", bias_mut_rate.getValue());
        trace.add("bias_gen_time", bias_gen_time.getValue());
        trace.add("tree_path", newick_path.getValue());
        trace.add("precision_matrix", precision_path.getValue());
        trace.add("branch_wise_correlation", branch_wise_correlation.getValue());
    }
};