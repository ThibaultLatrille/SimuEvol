#pragma once

#include <array>
#include <vector>
#include "fitness.hpp"


class AdditiveLandscape final : public FitnessLandscape {
  public:
    // The fitness profiles of amino-acids.
    std::vector<std::array<double, 20>> profiles;

    explicit AdditiveLandscape(std::vector<std::array<double, 20>> &profiles)
        : profiles{profiles} {}

    u_long nbr_sites() const override { return profiles.size(); }
};

class AdditiveState final : public FitnessState {
  private:
    AdditiveLandscape const &f;
    double sum_selection_coeff;

  public:
    std::unique_ptr<FitnessState> clone() const override {
        return std::make_unique<AdditiveState>(*this);
    };

    explicit AdditiveState(AdditiveLandscape const &f) : f{f} { sum_selection_coeff = 0; }

    bool operator==(FitnessState const &other) const override {
        return sum_selection_coeff ==
               dynamic_cast<AdditiveState const *>(&other)->sum_selection_coeff;
    };

    u_long nbr_sites() const override { return f.nbr_sites(); }

    void update(std::vector<char> const &codon_seq) override {
        sum_selection_coeff = 0;
        for (u_long i = 0; i < nbr_sites(); ++i) {
            sum_selection_coeff += f.profiles.at(i).at(codonLexico.codon_to_aa.at(codon_seq.at(i)));
        }
    }


    void update(std::vector<char> const &codon_seq, u_long site, char codon_to) override {
        sum_selection_coeff -=
            f.profiles.at(site).at(codonLexico.codon_to_aa.at(codon_seq.at(site)));
        sum_selection_coeff += f.profiles.at(site).at(codonLexico.codon_to_aa.at(codon_to));
    }

    double selection_coefficient(
        std::vector<char> const &codon_seq, u_long site, char codon_to) const override {
        return f.profiles[site][codonLexico.codon_to_aa[codon_to]] -
               f.profiles[site][codonLexico.codon_to_aa[codon_seq[site]]];
    };

    std::array<double, 20> aa_selection_coefficients(
        std::vector<char> const &codon_seq, u_long site) const override {
        return f.profiles.at(site);
    }
};

class SequenceAdditiveModel : public FitnessModel {
  public:
    SequenceAdditiveModel(std::string const &file_name, double const &beta, u_long &exon_size)
        : FitnessModel() {
        std::vector<std::array<double, 20>> site_aa_coeffs;

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
                val *= beta;
            }

            site_aa_coeffs.push_back(fitness_profil);
        }
        if (exon_size == 0) { exon_size = site_aa_coeffs.size(); }
        auto dv = std::div(static_cast<long>(site_aa_coeffs.size()), exon_size);
        if (dv.rem != 0) { dv.quot++; }
        fitness_landscapes.reserve(dv.quot);
        fitness_states.reserve(dv.quot);
        for (int exon{0}; exon < dv.quot; exon++) {
            size_t begin_exon = exon * exon_size;
            size_t end_exon = std::min(begin_exon + exon_size, site_aa_coeffs.size());

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
