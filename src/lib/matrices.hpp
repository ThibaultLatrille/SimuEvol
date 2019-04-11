#pragma once

#include <cassert>
#include <random>
#include "Eigen/Dense"
#include "lib/tree.hpp"
#include "statistic.hpp"

typedef Eigen::Matrix<double, 3, 3> Matrix3x3;
typedef Eigen::Matrix<double, 4, 4> Matrix4x4;
typedef Eigen::Matrix<double, 3, 1> Vector3x1;
typedef Eigen::Matrix<double, 4, 1> Vector4x1;
typedef Eigen::MatrixXd EMatrix;
typedef Eigen::VectorXd EVector;

int to_int(char c) { return c - '0'; }

double GetEntropy(const Vector4x1 &profile) {
    double tot = 0;
    for (Eigen::Index i = 0; i < profile.size(); i++) {
        tot -= (profile[i] < 1e-6) ? 0 : profile[i] * log(profile[i]);
    }
    return tot;
}


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

        (*this) *= mutation_rate;

        std::cout << "The nucleotide rate matrix (" << Codon::nucleotides << ") is:\n"
                  << *this << std::endl;
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
        bool reversible = is_reversible();
        trace.add("NucleotideMatrixReversible", reversible);

        trace.add("NucStatEntropy", GetEntropy(nuc_frequencies));
        if (reversible) {
            std::cout << "The nucleotide rate matrix is time-reversible." << std::endl;
        } else {
            std::cerr << "The nucleotide rate matrix is not time-reversible." << std::endl;
        }
    }
};

static int dim_population_size{0};
static int dim_mutation_rate_per_generation{1};
static int dim_generation_time{2};
static int dimensions{3};

class LogMultivariate : public Vector3x1 {
  public:
    explicit LogMultivariate() : Vector3x1(Vector3x1::Zero()) {}

    explicit LogMultivariate(
        u_long population_size, double mutation_rate_per_generation, double generation_time)
        : Vector3x1(Vector3x1::Zero()) {
        set_population_size(population_size);
        set_mutation_rate_per_generation(mutation_rate_per_generation);
        set_generation_time(generation_time);
    }

    explicit LogMultivariate(
        double beta, double mutation_rate_per_generation, double generation_time)
        : Vector3x1(Vector3x1::Zero()) {
        set_beta(beta);
        set_mutation_rate_per_generation(mutation_rate_per_generation);
        set_generation_time(generation_time);
    }

    void set_population_size(u_long population_size) {
        (*this)(dim_population_size) = std::log(population_size);
    }
    void set_beta(double beta) { (*this)(dim_population_size) = std::log(beta); }
    void set_mutation_rate_per_generation(double mutation_rate_per_generation) {
        (*this)(dim_mutation_rate_per_generation) = std::log(mutation_rate_per_generation);
    }
    void set_generation_time(double generation_time) {
        (*this)(dim_generation_time) = std::log(generation_time);
    }

    u_long population_size() const {
        return static_cast<u_long>(std::exp((*this)(dim_population_size)));
    }
    double beta() const { return std::exp((*this)(dim_population_size)); }
    double mutation_rate_per_generation() const {
        return std::exp((*this)(dim_mutation_rate_per_generation));
    }
    double generation_time() const { return std::exp((*this)(dim_generation_time)); }
};

class PieceWiseMultivariate {
  private:
    std::vector<double> times;
    std::vector<LogMultivariate> logmultivariates;

  public:
    PieceWiseMultivariate() = default;
    ~PieceWiseMultivariate() = default;

    void AddMultivariate(double time_elapsed, LogMultivariate const &multivariate) {
        if (!logmultivariates.empty() and logmultivariates.back() == multivariate) {
            times.back() += time_elapsed;
        } else {
            times.push_back(time_elapsed);
            logmultivariates.push_back(multivariate);
        }
    }

    void clear() {
        times.clear();
        logmultivariates.clear();
    }

    double ArithmeticPopSize() const { return ArithmeticDim(dim_population_size); }
    double ArithmeticMutRate() const { return ArithmeticDim(dim_mutation_rate_per_generation); }
    double ArithmeticGenTime() const { return ArithmeticDim(dim_generation_time); }

    double GeometricPopSize() const { return GeometricDim(dim_population_size); }
    double GeometricMutRate() const { return GeometricDim(dim_mutation_rate_per_generation); }
    double GeometricGenTime() const { return GeometricDim(dim_generation_time); }

    double HarmonicPopSize() const { return HarmonicDim(dim_population_size); }

    double ArithmeticDim(int dim) const {
        double tot_log{0}, tot_time{0};
        for (size_t i = 0; i < times.size(); i++) {
            tot_time += times[i];
            tot_log += times[i] * exp(logmultivariates[i](dim));
        }
        return tot_log / tot_time;
    }

    double GeometricDim(int dim) const {
        double tot_log{0}, tot_time{0};
        for (size_t i = 0; i < times.size(); i++) {
            tot_time += times[i];
            tot_log += times[i] * logmultivariates[i](dim);
        }
        return exp(tot_log / tot_time);
    }

    double HarmonicDim(int dim) const {
        double tot_inv{0}, tot_time{0};
        for (size_t i = 0; i < times.size(); i++) {
            tot_time += times[i];
            tot_inv += times[i] * exp(-logmultivariates[i](dim));
        }
        return tot_time / tot_inv;
    }

    void add_to_tree(Tree &tree, Tree::NodeIndex node, double geom_pop_size) const {
        tree.set_tag(node, "Branch_population_size", d_to_string(geom_pop_size));
        tree.set_tag(node, "Branch_LogNe", d_to_string(log10(geom_pop_size)));
        tree.set_tag(node, "Branch_generation_time", d_to_string(GeometricGenTime()));
        tree.set_tag(node, "Branch_mutation_rate_per_generation", d_to_string(GeometricMutRate()));

        tree.set_tag(node, "Branch_arithmetic_population_size", d_to_string(ArithmeticPopSize()));
        tree.set_tag(node, "Branch_arithmetic_generation_time", d_to_string(ArithmeticGenTime()));
        tree.set_tag(node, "Branch_arithmetic_mutation_rate_per_generation",
            d_to_string(ArithmeticMutRate()));

        tree.set_tag(node, "Branch_harmonic_population_size", d_to_string(HarmonicPopSize()));
    }
};

class CorrelationMatrix : public Matrix3x3 {
  public:
    Matrix3x3 precision_matrix = Matrix3x3::Zero();

    explicit CorrelationMatrix() : Matrix3x3(Matrix3x3::Zero()) {}

    explicit CorrelationMatrix(std::string const &precision_filename, bool fix_pop_size,
        bool fix_mut_rate, bool fix_gen_time)
        : Matrix3x3(Matrix3x3::Zero()) {
        if (fix_pop_size and fix_mut_rate and fix_gen_time) {
            std::cout << "The correlation matrix is 0 (all fixed effects)" << std::endl;
            return;
        }

        std::ifstream input_stream(precision_filename);

        if (precision_filename.empty()) {
            std::cerr << "No precision matrix file provided, use the default precision matrix."
                      << std::endl;
            return;
        }

        if (!input_stream) {
            std::cerr << "Precision matrix file " << precision_filename
                      << " doesn't exist, use the default precision matrix instead." << std::endl;
            return;
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
            precision_matrix(i, j) = v;
            if (i != j) {
                precision_matrix(j, i) = v;
            } else {
                assert(v >= 0.0);
            }
        }
        std::cout << "The input precision matrix is:\n" << precision_matrix << std::endl;
        Matrix3x3 correlation_matrix = precision_matrix.inverse();
        assert(correlation_matrix.transpose() == correlation_matrix);
        for (int i = 0; i < dimensions; i++) {
            if (fix_pop_size and i == dim_population_size) { continue; }
            if (fix_mut_rate and i == dim_mutation_rate_per_generation) { continue; }
            if (fix_gen_time and i == dim_generation_time) { continue; }
            for (int j = 0; j < dimensions; j++) {
                if (fix_pop_size and j == dim_population_size) { continue; }
                if (fix_mut_rate and j == dim_mutation_rate_per_generation) { continue; }
                if (fix_gen_time and j == dim_generation_time) { continue; }
                (*this)(i, j) = correlation_matrix(i, j);
            }
        }
        std::cout << "The correlation (taking fixed effect into account) matrix is:\n"
                  << *this << std::endl;
    }

    void add_to_trace(Trace &trace) const {
        for (int i = 0; i < dimensions; i++) {
            for (int j = 0; j <= i; j++) {
                trace.add("Precision_" + std::to_string(i) + "_" + std::to_string(j),
                    precision_matrix.coeffRef(i, j));
            }
        }
    }
};

// Method computing the equilibrium frequencies for one site.
std::array<double, 64> codon_frequencies(std::array<double, 20> const &aa_fitness_profil,
    NucleotideRateMatrix const &nuc_matrix, double beta) {
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
                codon_freq * exp(aa_fitness_profil[Codon::codon_to_aa_array[codon]] * beta);
        } else {
            codon_frequencies[codon] = 0.;
        }
    }

    double sum_freq = std::accumulate(codon_frequencies.begin(), codon_frequencies.end(), 0.0);
    for (char codon{0}; codon < 64; codon++) { codon_frequencies[codon] /= sum_freq; }

    return codon_frequencies;
};

double rate_fixation(
    std::array<double, 20> const &preferences, char codon_from, char codon_to, double scale) {
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
std::tuple<double, double> predicted_dn_dn0(
    std::vector<std::array<double, 20>> const &aa_fitness_profiles,
    NucleotideRateMatrix const &mutation_rate_matrix, double scale) {
    // For all site of the sequence.
    double dn{0.}, dn0{0.};
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
                        double rate = codon_freqs[codon_from] * mutation_rate_matrix(n_from, n_to);
                        dn0 += rate;
                        rate *= rate_fixation(aa_fitness_profile, codon_from, codon_to, scale);
                        dn += rate;
                    }
                }
            }
        }
    }
    return std::make_tuple(dn, dn0);
}

// Theoretical computation of the predicted omega
std::tuple<double, double> flow_dn_dn0(
    std::vector<std::array<double, 20>> const &aa_fitness_profiles,
    std::vector<char> const &codon_seq, NucleotideRateMatrix const &mutation_rate_matrix,
    double scale = 1.0) {
    assert(aa_fitness_profiles.size() == codon_seq.size());
    double dn{0.}, dn0{0.};
    // For all site of the sequence.
    for (size_t site{0}; site < aa_fitness_profiles.size(); site++) {
        // For all possible neighbors.
        for (auto &neighbor : Codon::codon_to_neighbors_array[codon_seq[site]]) {
            // Codon after mutation, Nucleotide original and Nucleotide after mutation.
            char codon_to{0}, n_from{0}, n_to{0};
            std::tie(codon_to, n_from, n_to) = neighbor;

            // If the mutated amino-acid is a stop codon, the rate of fixation is 0.
            // Else, if the mutated and original amino-acids are non-synonymous, we
            // compute the rate of fixation. Note that, if the mutated and original
            // amino-acids are synonymous, the rate of fixation is 1.
            if (Codon::codon_to_aa_array[codon_to] != 20 and
                Codon::codon_to_aa_array[codon_seq[site]] != Codon::codon_to_aa_array[codon_to]) {
                double rate = mutation_rate_matrix(n_from, n_to);
                dn0 += rate;
                rate *= rate_fixation(aa_fitness_profiles[site], codon_seq[site], codon_to, scale);
                dn += rate;
            }
        }
    }
    return std::make_tuple(dn, dn0);
}