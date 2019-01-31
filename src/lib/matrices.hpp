#pragma once

#include <cassert>
#include <random>
#include "Eigen/Dense"
#include "statistic.hpp"

typedef Eigen::Matrix<double, 3, 3> Matrix3x3;
typedef Eigen::Matrix<double, 4, 4> Matrix4x4;
typedef Eigen::Matrix<double, 3, 1> Vector3x1;
typedef Eigen::Matrix<double, 4, 1> Vector4x1;
typedef Eigen::MatrixXd EMatrix;
typedef Eigen::VectorXd EVector;

int to_int(char c) { return c - '0'; }

class NucleotideRateMatrix : public Matrix4x4 {
  public:
    // Mutation
    Vector4x1 nuc_frequencies;
    double mutation_rate;

    Vector4x1 sum_mutation_rates;
    double max_sum_mutation_rates;
    mutable std::uniform_real_distribution<double> max_real_distr;
    mutable std::array<std::discrete_distribution<char>, 4> mutation_distr;

    explicit NucleotideRateMatrix(
        std::string const &input_filename, double mu = 1.0, bool normalize = true)
        : Matrix4x4(Matrix4x4::Ones()), mutation_rate{mu} {
        std::ifstream input_stream(input_filename);

        if (input_filename.empty()) {
            std::cerr << "No rate matrix file provided, use the default nucleotide rate matrix."
                      << std::endl;
        }

        if (!input_stream) {
            std::cerr << "Rate matrix file " << input_filename
                      << " doesn't exist, use the default rate matrix nucleotide instead."
                      << std::endl;
        }

        std::string line;
        getline(input_stream, line);
        std::stringstream header_stream(line);
        getline(input_stream, line);
        std::stringstream values_stream(line);
        std::string header_word, value;
        while (getline(header_stream, header_word, '\t')) {
            assert(header_word.size() == 4);
            assert(header_word.substr(0, 2) == "q_");
            getline(values_stream, value, '\t');
            (*this)(Codon::nuc_to_index.at(header_word[2]),
                Codon::nuc_to_index.at(header_word[3])) = stod(value);
        }

        for (int diag{0}; diag < 4; diag++) { (*this)(diag, diag) = 0.0; }
        (*this) -= this->rowwise().sum().asDiagonal();

        equilibrium_frequencies();
        std::cout << "The equilibrium nucleotides frequencies (" << Codon::nucleotides << ") are:\n"
                  << nuc_frequencies.transpose() << std::endl;

        if (normalize) { normalize_matrix(); }

        (*this) *= mutation_rate;

        std::cout << "The nucleotide rate matrix (" << Codon::nucleotides << ") is:\n"
                  << *this << std::endl;

        sum_mutation_rates = -(*this).diagonal();
        max_sum_mutation_rates = sum_mutation_rates.maxCoeff();

        max_real_distr = std::uniform_real_distribution<double>(0.0, max_sum_mutation_rates);
        for (char nuc_from{0}; nuc_from < 4; nuc_from++) {
            std::array<double, 4> rates{0.0};
            for (char nuc_to{0}; nuc_to < 4; nuc_to++) {
                if (nuc_from != nuc_to) { rates[nuc_to] = (*this)(nuc_from, nuc_to); }
            }
            mutation_distr[nuc_from] = std::discrete_distribution<char>(rates.begin(), rates.end());
        }
    }

    void set_mutation_rate(double mu) {
        (*this) *= (mu / mutation_rate);
        mutation_rate = mu;
    }


    bool is_reversible() const {
        bool reversible = true;
        for (int row{0}; row < 4; row++) {
            for (int col{0}; col < 4; col++) {
                double diff = nuc_frequencies(row) * (*this)(row, col) -
                              nuc_frequencies(col) * (*this)(col, row);
                if (std::abs(diff) > 1e-6) { reversible = false; }
            }
        }
        if (reversible) {
            std::cout << "The nucleotide rate matrix is time-reversible." << std::endl;
        } else {
            std::cerr << "The nucleotide rate matrix is not time-reversible." << std::endl;
        }
        return reversible;
    }

    // Compute the kernel of the mutation-rate matrix.
    // This kernel is a vector of the nucleotides frequencies (not normalized to 1) at equilibrium.
    void equilibrium_frequencies() {
        Eigen::Matrix<double, 4, Eigen::Dynamic> kernel = this->transpose().fullPivLu().kernel();

        assert(kernel.cols() == 1);
        nuc_frequencies = kernel.col(0) / kernel.col(0).sum();
    }

    void normalize_matrix() {
        double events = -(nuc_frequencies.transpose() * (*this).diagonal()).sum();
        (*this) /= events;
    }

    void add_to_trace(Trace &trace) const {
        for (char nuc_from{0}; nuc_from < 4; nuc_from++) {
            for (char nuc_to{0}; nuc_to < 4; nuc_to++) {
                if (nuc_from != nuc_to) {
                    std::string key = "q_";
                    key += Codon::nucleotides[nuc_from];
                    key += Codon::nucleotides[nuc_to];
                    trace.add(key, (*this)(nuc_from, nuc_to));
                }
            }
        }
        trace.add("NucleotideMatrixReversible", is_reversible());
    }
};

class LogMultivariate : public Vector3x1 {
    int dimensions = 3;

  public:
    explicit LogMultivariate(unsigned population_size, double generation_time, double mu)
        : Vector3x1(Vector3x1::Zero()) {
        set_population_size(population_size);
        set_generation_time(generation_time);
        set_mu(mu);
    }

    void set_population_size(unsigned population_size) { (*this)(0) = std::log(population_size); }
    void set_generation_time(double generation_time) { (*this)(1) = std::log(generation_time); }
    void set_mu(double mu) { (*this)(2) = std::log(mu); }

    unsigned population_size() { return static_cast<unsigned>(std::exp((*this)(0))); }
    double generation_time() { return std::exp((*this)(1)); }
    double mu() { return std::exp((*this)(2)); }

    void add_to_trace(Trace &trace) const {
        for (int i = 0; i < dimensions; i++) {
            trace.add("logmultivariate_" + std::to_string(i), (*this).coeffRef(i));
        }
    }
};

class CorrelationMatrix : public Matrix3x3 {
  public:
    int dimensions = 3;

    explicit CorrelationMatrix(std::string const &input_filename) : Matrix3x3(Matrix3x3::Zero()) {
        std::ifstream input_stream(input_filename);

        if (input_filename.empty()) {
            std::cerr << "No correlation matrix file provided, use the default correlation matrix."
                      << std::endl;
        }

        if (!input_stream) {
            std::cerr << "Correlation matrix file " << input_filename
                      << " doesn't exist, use the default correlation matrix instead." << std::endl;
        }

        std::string line;
        getline(input_stream, line);
        std::stringstream header_stream(line);
        getline(input_stream, line);
        std::stringstream values_stream(line);
        std::string header_word, value;
        while (getline(header_stream, header_word, '\t')) {
            assert(header_word.size() == 4);
            assert(header_word.substr(0, 2) == "c_");
            getline(values_stream, value, '\t');
            int i = to_int(header_word[2]), j = to_int(header_word[3]);
            double v = stod(value);
            (*this)(i, j) = v;
            if (i != j) {
                (*this)(j, i) = v;
            } else {
                assert(v >= 0.0);
            }
        }
        std::cout << "The correlation matrix is:\n" << *this << std::endl;
        assert(this->transpose() == (*this));
    }

    void add_to_trace(Trace &trace) const {
        for (int i = 0; i < dimensions; i++) {
            for (int j = 0; j <= i; j++) {
                trace.add(
                    "cov_" + std::to_string(i) + "_" + std::to_string(j), (*this).coeffRef(i, j));
            }
        }
    }
};

// Method computing the equilibrium frequencies for one site.
std::array<double, 64> codon_frequencies(
    std::array<double, 20> const &aa_fitness_profil, NucleotideRateMatrix const &nuc_matrix, double scale = 1.0) {
    std::array<double, 64> codon_frequencies{0};
    // For each site of the vector of the site frequencies.
    for (char codon{0}; codon < 64; codon++) {
        double codon_freq = 1.0;

        // For all nucleotides in the codon
        for (auto const &nuc : Codon::codon_to_triplet_array[codon]) {
            codon_freq *= nuc_matrix.nuc_frequencies[nuc];
        }

        if (Codon::codon_to_aa_array[codon] != 20) {
            codon_frequencies[codon] =
                codon_freq * exp(aa_fitness_profil[Codon::codon_to_aa_array[codon]] * scale);
        } else {
            codon_frequencies[codon] = 0.;
        }
    }

    double sum_freq = std::accumulate(codon_frequencies.begin(), codon_frequencies.end(), 0.0);
    for (char codon{0}; codon < 64; codon++) { codon_frequencies[codon] /= sum_freq; }

    return codon_frequencies;
};

double rate_fixation(std::array<double, 20> const &preferences, char codon_from, char codon_to, double scale = 1.0) {
    double rate_fix = 1.0;
    // Selective strength between the mutated and original amino-acids.
    double s{0.};
    s = preferences[Codon::codon_to_aa_array[codon_to]];
    s -= preferences[Codon::codon_to_aa_array[codon_from]];
    s *= scale;
    // If the selective strength is 0, the rate of fixation is neutral.
    // Else, the rate of fixation is computed using population genetic formulas
    // (Kimura).
    if (fabs(s) > Codon::epsilon) {
        // The substitution rate is the mutation rate multiplied by the rate of
        // fixation.
        rate_fix = s / (1 - exp(-s));
    }
    return rate_fix;
}

// Theoretical computation of the predicted omega
std::tuple<double, double> predicted_dn_d0(
    std::vector<std::array<double, 20>> const &aa_fitness_profiles,
    NucleotideRateMatrix const &mutation_rate_matrix, double scale = 1.0) {
    // For all site of the sequence.
    double dn{0.}, d0{0.};
    for (auto const &aa_fitness_profile : aa_fitness_profiles) {
        // Codon original before substitution.
        std::array<double, 64> codon_freqs =
            codon_frequencies(aa_fitness_profile, mutation_rate_matrix, scale);

        for (char codon_from{0}; codon_from < 64; codon_from++) {
            if (Codon::codon_to_aa_array[codon_from] != 20) {
                // For all possible neighbors.
                for (auto &neighbor : Codon::codon_to_neighbors_array[codon_from]) {
                    // Codon after mutation, Nucleotide original and Nucleotide after mutation.
                    char codon_to{0}, n_from{0}, n_to{0};
                    std::tie(codon_to, n_from, n_to) = neighbor;

                    // If the mutated amino-acid is a stop codon, the rate of fixation is 0.
                    // Else, if the mutated and original amino-acids are non-synonymous, we
                    // compute the rate of fixation. Note that, if the mutated and original
                    // amino-acids are synonymous, the rate of fixation is 1.
                    if (Codon::codon_to_aa_array[codon_to] != 20 and
                        Codon::codon_to_aa_array[codon_from] !=
                            Codon::codon_to_aa_array[codon_to]) {
                        dn += codon_freqs[codon_from] * mutation_rate_matrix(n_from, n_to) *
                              rate_fixation(aa_fitness_profile, codon_from, codon_to, scale);
                        d0 += codon_freqs[codon_from] * mutation_rate_matrix(n_from, n_to);
                    }
                }
            }
        }
    }
    return std::make_tuple(dn, d0);
}

// Theoretical computation of the predicted omega
double predicted_omega(std::vector<std::array<double, 20>> const &aa_fitness_profiles,
    NucleotideRateMatrix const &mutation_rate_matrix) {
    // For all site of the sequence.
    double dn{0.}, d0{0.};
    std::tie(dn, d0) = predicted_dn_d0(aa_fitness_profiles, mutation_rate_matrix);
    return dn / d0;
}