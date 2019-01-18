#pragma once

#include <cassert>
#include <fstream>
#include "tclap/CmdLine.h"

class OutputArgParse {
  protected:
    TCLAP::CmdLine& cmd;

  public:
    explicit OutputArgParse(TCLAP::CmdLine& cmd) : cmd{cmd} {}

    TCLAP::ValueArg<std::string> output_path{
        "o", "output", "output protein path", true, "chain", "string", cmd};
};

class SimuArgParse : public OutputArgParse {
  public:
    explicit SimuArgParse(TCLAP::CmdLine& cmd) : OutputArgParse(cmd) {}

    TCLAP::ValueArg<std::string> preferences_path{
        "f", "preferences", "input site-specific preferences path", true, "chain", "string", cmd};
    TCLAP::ValueArg<std::string> newick_path{
        "t", "newick", "input newick tree path", true, "chain", "string", cmd};
};

std::string open_newick(std::string const& file_name) {
    std::ifstream input_stream(file_name);
    if (!input_stream) std::cerr << "Can't open newick file!" << std::endl;

    std::string line;
    getline(input_stream, line);

    return line;
}

std::vector<std::array<double, 20>> open_preferences(
    std::string const& file_name, double const& beta) {
    std::vector<std::array<double, 20>> fitness_profiles{0};

    std::ifstream input_stream(file_name);
    if (!input_stream) std::cerr << "Can't open preferences file!" << std::endl;

    std::string line;

    // skip the header of the file
    getline(input_stream, line);

    while (getline(input_stream, line)) {
        std::array<double, 20> fitness_profil{0};
        std::string word;
        istringstream line_stream(line);
        unsigned counter{0};

        while (getline(line_stream, word, ' ')) {
            if (counter > 2) { fitness_profil[counter - 3] = beta * log(stod(word)); }
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