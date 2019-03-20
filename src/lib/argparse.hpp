#pragma once

#include <cassert>
#include <cmath>
#include <fstream>
#include "tclap/CmdLine.h"

class OutputArgParse {
  protected:
    TCLAP::CmdLine& cmd;

  public:
    explicit OutputArgParse(TCLAP::CmdLine& cmd) : cmd{cmd} {}

    TCLAP::ValueArg<std::string> output_path{"o", "output", "output path", true, "", "string", cmd};
};

class SimuArgParse : public OutputArgParse {
  public:
    explicit SimuArgParse(TCLAP::CmdLine& cmd) : OutputArgParse(cmd) {}

    TCLAP::ValueArg<std::string> correlation_path{
        "c", "correlation_matrix", "input correlation matrix path", false, "", "string", cmd};
    TCLAP::ValueArg<double> mu{
        "m", "mu", "Mutation rate (at the root)", false, 1e-8, "double", cmd};
    TCLAP::ValueArg<double> root_age{
        "a", "root_age", "Age of the root", false, 50e6, "double", cmd};
    TCLAP::ValueArg<double> generation_time{"g", "generation_time",
        "Number of year between generations (at the root)", false, 40, "double", cmd};
    TCLAP::ValueArg<std::string> nuc_matrix_path{
        "q", "nuc_matrix", "input nucleotide matrix preferences path", false, "", "string", cmd};
    TCLAP::ValueArg<std::string> preferences_path{
        "f", "preferences", "input site-specific preferences path", true, "", "string", cmd};
    TCLAP::ValueArg<std::string> newick_path{
        "t", "newick", "input newick tree path", true, "", "string", cmd};
    TCLAP::ValueArg<double> beta{
        "b", "beta", "Stringency parameter of the fitness profiles", false, 0.0, "double", cmd};
    TCLAP::ValueArg<u_long> exons{"s", "exon_size",
        "Number of codon sites per exon (default 0 means the size of the fitness profiles "
        "provided, thus assuming complete linkage between sites)",
        false, 0, "u_long", cmd};
};

std::vector<std::array<double, 20>> open_preferences(
    std::string const& file_name, double const& beta) {
    std::vector<std::array<double, 20>> fitness_profiles{0};

    std::ifstream input_stream(file_name);
    if (!input_stream)
        std::cerr << "Preferences file " << file_name << " doesn't exist" << std::endl;

    std::string line;

    // skip the header of the file
    getline(input_stream, line);
    char sep{' '};
    u_long nbr_col = 0;
    for (char sep_test : std::vector<char>({' ', ',', '\t'})) {
        u_long n = static_cast<u_long>(std::count(line.begin(), line.end(), sep_test));
        if (n > nbr_col) {
            sep = sep_test;
            nbr_col = n + 1;
        }
    }
    nbr_col -= 20;

    while (getline(input_stream, line)) {
        std::array<double, 20> fitness_profil{0};
        std::string word;
        istringstream line_stream(line);
        u_long counter{0};

        while (getline(line_stream, word, sep)) {
            if (counter > nbr_col) {
                fitness_profil[counter - (nbr_col + 1)] = beta * std::log(stod(word));
            }
            counter++;
        }

        fitness_profiles.push_back(fitness_profil);
    }
    return fitness_profiles;
}

std::string char_to_str(char const& _char) {
    std::string _str(1, _char);
    return _str;
}