#pragma once

#include <array>
#include <vector>
#include "fitness.hpp"
#include "tclap/CmdLine.h"

class AdditiveLandscape final : public FitnessLandscape {
  public:
    // The fitness profiles of amino-acids.
    std::vector<std::array<double, 20>> profiles;

    explicit AdditiveLandscape(std::vector<std::array<double, 20>> &profiles) : profiles{profiles} {
        additive_phenotype = true;
        approx = true;
    }

    u_long nbr_sites() const override { return profiles.size(); }
};

class AdditiveState final : public FitnessState {
  private:
    AdditiveLandscape const &f;

  public:
    std::unique_ptr<FitnessState> clone() const override {
        return std::make_unique<AdditiveState>(*this);
    };

    explicit AdditiveState(AdditiveLandscape const &f) : FitnessState(f), f{f} {}

    bool operator==(FitnessState const &other) const override {
        return log_fitness == other.log_fitness;
    };

    u_long nbr_sites() const override { return f.nbr_sites(); }

    void update(std::vector<char> const &codon_seq, double const &pop_size) override {
        log_fitness = 0;
        for (u_long i = 0; i < nbr_sites(); ++i) {
            log_fitness += f.profiles.at(i).at(codonLexico.codon_to_aa.at(codon_seq.at(i)));
        }
    }

    void update(std::vector<char> const &codon_seq, u_long site, char codon_to, bool burn_in,
        double const &pop_size) override {
        log_fitness -= f.profiles.at(site).at(codonLexico.codon_to_aa.at(codon_seq.at(site)));
        log_fitness += f.profiles.at(site).at(codonLexico.codon_to_aa.at(codon_to));
        if (!burn_in) { summary_stats["sub-log-fitness"].add(log_fitness); }
    }

    double selection_coefficient(std::vector<char> const &codon_seq, u_long site, char codon_to,
        bool burn_in, double const &pop_size) const override {
        double s = f.profiles[site][codonLexico.codon_to_aa[codon_to]] -
                   f.profiles[site][codonLexico.codon_to_aa[codon_seq[site]]];
        if (!burn_in) {
            summary_stats["mut-s"].add(s);
            summary_stats["mut-log-fitness"].add(log_fitness + s);
        }
        return s;
    };

    std::array<double, 20> aa_selection_coefficients(
        std::vector<char> const &codon_seq, u_long site, double const &pop_size) const override {
        return f.profiles.at(site);
    }
};

std::unordered_map<std::string, SummaryStatistic> FitnessState::summary_stats = {};

class AdditiveArgParse {
  protected:
    TCLAP::CmdLine &cmd;

  public:
    explicit AdditiveArgParse(TCLAP::CmdLine &cmd) : cmd{cmd} {}
    TCLAP::ValueArg<std::string> preferences_path{
        "", "preferences", "input site-specific preferences path", true, "", "string", cmd};
    TCLAP::ValueArg<double> beta{
        "", "beta", "Stringency parameter of the fitness profiles", false, 1.0, "double", cmd};

    void add_to_trace(Trace &trace) {
        trace.add("site_preferences_path", preferences_path.getValue());
        trace.add("site_preferences_beta", beta.getValue());
    }
};

class SequenceAdditiveModel : public FitnessModel {
  public:
    SequenceAdditiveModel(double const &beta_prefactor, u_long &exon_size, AdditiveArgParse &args)
        : FitnessModel() {
        std::vector<std::array<double, 20>> site_aa_coeffs;

        std::ifstream input_stream(args.preferences_path.getValue());
        if (!input_stream)
            std::cerr << "Preferences file " << args.preferences_path.getValue() << " doesn't exist"
                      << std::endl;

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

            double min_val = 0;
            while (getline(line_stream, word, sep)) {
                if (counter >= nbr_col) {
                    double val = std::log(stod(word));
                    if (val < min_val) { min_val = val; }
                    fitness_profil[counter - nbr_col] = val;
                }
                counter++;
            }
            for (auto &val : fitness_profil) {
                val -= min_val;
                val *= beta_prefactor * args.beta.getValue();
            }

            site_aa_coeffs.push_back(fitness_profil);
        }
        if (exon_size == 0) { exon_size = site_aa_coeffs.size(); }
        auto dv = std::div(static_cast<long>(site_aa_coeffs.size()), exon_size);
        if (dv.rem != 0) { dv.quot++; }
        fitness_landscapes.reserve(dv.quot);
        fitness_states.reserve(dv.quot);
        for (int exon{0}; exon < dv.quot; exon++) {
            std::size_t begin_exon = exon * exon_size;
            std::size_t end_exon = std::min(begin_exon + exon_size, site_aa_coeffs.size());

            std::vector<std::array<double, 20>> exon_site_aa_coeffs(
                site_aa_coeffs.begin() + begin_exon, site_aa_coeffs.begin() + end_exon);
            fitness_landscapes.emplace_back(
                std::make_unique<AdditiveLandscape>(exon_site_aa_coeffs));

            fitness_states.emplace_back(std::make_unique<AdditiveState>(
                *dynamic_cast<AdditiveLandscape *>(fitness_landscapes.at(exon).get())));
        }
        assert(nbr_sites() == site_aa_coeffs.size());
        std::cout << fitness_landscapes.size() << " exons created." << std::endl;
        if (dv.rem != 0) { std::cout << "Last exon is " << dv.rem << " sites long." << std::endl; }
    }
};
