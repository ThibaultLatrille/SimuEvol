#pragma once

#include "matrices.hpp"
#include "tclap/CmdLine.h"

double Pfix(double const &pop_size, double const &selection_coefficient) {
    double S = 4 * selection_coefficient * pop_size;
    if ((abs(S)) < 1e-4) {
        return 1 + S / 2;
    } else {
        // return 2 * pop_size * (1.0 - exp(-2.0 * s)) / (1.0 - exp(-4 * pop_size * s));
        return S / (1.0 - exp(-S));
    }
}


// Method computing the equilibrium frequencies for one site.
std::array<double, 64> codon_frequencies(std::array<double, 20> const &aa_selection_coefficient,
    NucleotideRateMatrix const &nuc_matrix, double const &pop_size) {
    std::array<double, 64> codon_frequencies{0};
    // For each site of the vector of the site frequencies.
    for (char codon{0}; codon < 64; codon++) {
        double codon_freq = 1.0;

        // For all nucleotides in the codon
        for (auto const &nuc : codonLexico.codon_to_triplet[codon]) {
            codon_freq *= nuc_matrix.nuc_frequencies[nuc];
        }

        if (codonLexico.codon_to_aa[codon] != 20) {
            codon_frequencies[codon] =
                codon_freq *
                exp(4 * aa_selection_coefficient[codonLexico.codon_to_aa[codon]] * pop_size);
        } else {
            codon_frequencies[codon] = 0.;
        }
    }

    double sum_freq = std::accumulate(codon_frequencies.begin(), codon_frequencies.end(), 0.0);
    for (char codon{0}; codon < 64; codon++) { codon_frequencies[codon] /= sum_freq; }

    return codon_frequencies;
}

class FitnessLandscape {
  public:
    virtual ~FitnessLandscape() = default;

    virtual u_long nbr_sites() const = 0;
};

class FitnessState {
  public:
    static std::unordered_map<std::string, SummaryStatistic> summary_stats;

    virtual std::unique_ptr<FitnessState> clone() const = 0;

    virtual bool operator==(FitnessState const &other) const = 0;

    virtual u_long nbr_sites() const = 0;

    virtual void update(std::vector<char> const &codon_seq, double const &pop_size) = 0;

    virtual void update(
        std::vector<char> const &codon_seq, u_long site, char codon_to, bool burn_in, double const &pop_size) = 0;

    virtual double selection_coefficient(
        std::vector<char> const &codon_seq, u_long site, char codon_to, bool burn_in, double const &pop_size) const = 0;

    virtual std::array<double, 20> aa_selection_coefficients(
        std::vector<char> const &codon_seq, u_long site, double const &pop_size) const {
        std::array<double, 20> aa_sel_coeffs{};
        for (char aa = 0; aa < 20; aa++) {
            auto it = std::find(codonLexico.codon_to_aa.begin(), codonLexico.codon_to_aa.end(), aa);
            assert(it != codonLexico.codon_to_aa.end());
            char codon_to = std::distance(codonLexico.codon_to_aa.begin(), it);
            assert(codonLexico.codon_to_aa[codon_to] == aa);
            aa_sel_coeffs[aa] = selection_coefficient(codon_seq, site, codon_to, true, pop_size);
        }
        return aa_sel_coeffs;
    }

    // Theoretical computation of the predicted omega
    std::tuple<double, double> predicted_dn_dn0(std::vector<char> const &codon_seq,
        NucleotideRateMatrix const &mutation_rate_matrix, double pop_size) const {
        // For all site of the sequence.
        double dn{0.}, dn0{0.};
        for (u_long site = 0; site < nbr_sites(); ++site) {
            auto aa_sel_coeffs = aa_selection_coefficients(codon_seq, site, pop_size);
            // Codon original before substitution.
            std::array<double, 64> codon_freqs =
                codon_frequencies(aa_sel_coeffs, mutation_rate_matrix, pop_size);

            for (char codon_from{0}; codon_from < 64; codon_from++) {
                if (codonLexico.codon_to_aa[codon_from] != 20) {
                    // For all possible neighbors.
                    for (auto &neighbor : codonLexico.codon_to_neighbors[codon_from]) {
                        // Codon after mutation, Nucleotide original and Nucleotide after mutation.
                        char codon_to{0}, n_from{0}, n_to{0};
                        std::tie(codon_to, n_from, n_to) = neighbor;

                        // If the mutated amino-acid is a stop codon, the rate of fixation is 0.
                        // Else, if the mutated and original amino-acids are non-synonymous, we
                        // compute the rate of fixation. Note that, if the mutated and original
                        // amino-acids are synonymous, the rate of fixation is 1.
                        if (codonLexico.codon_to_aa[codon_to] != 20 and
                            codonLexico.codon_to_aa[codon_from] !=
                                codonLexico.codon_to_aa[codon_to]) {
                            double rate =
                                codon_freqs[codon_from] * mutation_rate_matrix(n_from, n_to);
                            dn0 += rate;
                            double s = aa_sel_coeffs[codonLexico.codon_to_aa[codon_to]] -
                                       aa_sel_coeffs[codonLexico.codon_to_aa[codon_from]];
                            rate *= Pfix(pop_size, s);
                            dn += rate;
                        }
                    }
                }
            }
        }
        return std::make_tuple(dn, dn0);
    }

    // Theoretical computation of the predicted omega
    std::tuple<double, double> flow_dn_dn0(std::vector<char> const &codon_seq,
        NucleotideRateMatrix const &mutation_rate_matrix, double pop_size) const {
        assert(nbr_sites() == codon_seq.size());
        double dn{0.}, dn0{0.};
        // For all site of the sequence.
        for (u_long site = 0; site < nbr_sites(); ++site) {
            auto aa_sel_coeffs = aa_selection_coefficients(codon_seq, site, pop_size);
            // For all possible neighbors.
            for (auto &neighbor : codonLexico.codon_to_neighbors[codon_seq[site]]) {
                // Codon after mutation, Nucleotide original and Nucleotide after mutation.
                char codon_to{0}, n_from{0}, n_to{0};
                std::tie(codon_to, n_from, n_to) = neighbor;

                // If the mutated amino-acid is a stop codon, the rate of fixation is 0.
                // Else, if the mutated and original amino-acids are non-synonymous, we
                // compute the rate of fixation. Note that, if the mutated and original
                // amino-acids are synonymous, the rate of fixation is 1.
                if (codonLexico.codon_to_aa[codon_to] != 20 and
                    codonLexico.codon_to_aa[codon_seq[site]] != codonLexico.codon_to_aa[codon_to]) {
                    double rate = mutation_rate_matrix(n_from, n_to);
                    dn0 += rate;
                    double s = aa_sel_coeffs[codonLexico.codon_to_aa[codon_to]] -
                               aa_sel_coeffs[codonLexico.codon_to_aa[codon_seq[site]]];
                    rate *= Pfix(pop_size, s);
                    dn += rate;
                }
            }
        }
        return std::make_tuple(dn, dn0);
    }

    virtual ~FitnessState() = default;
};

class FitnessModel {
  public:
    std::vector<std::unique_ptr<FitnessLandscape>> fitness_landscapes;
    std::vector<std::unique_ptr<FitnessState>> fitness_states;

    u_long nbr_sites() const {
        return std::accumulate(fitness_landscapes.begin(), fitness_landscapes.end(), 0,
            [](u_long a, auto const &b) { return a + b->nbr_sites(); });
    };

    u_long nbr_exons() const { return fitness_landscapes.size(); };
};