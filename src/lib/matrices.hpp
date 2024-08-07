#pragma once

#include <cassert>
#include <random>
#include "Eigen/Dense"
#include "codon.hpp"
#include "io.hpp"
#include "statistic.hpp"
#include "tree.hpp"

typedef Eigen::Matrix<double, 3, 3> Matrix3x3;
typedef Eigen::Matrix<double, 4, 4> Matrix4x4;
typedef Eigen::Matrix<double, 4, 1> Vector4x1;
typedef Eigen::MatrixXd EMatrix;
typedef Eigen::VectorXd EVector;

int to_int(char c) { return c - '0'; }

double GetEntropy(Vector4x1 const &profile) {
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
    mutable std::array<std::discrete_distribution<short>, 4> mutation_distr;

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
            mutation_distr[nuc_from] = std::discrete_distribution<short>(rates.begin(), rates.end());
        }

        (*this) *= mutation_rate;

        std::cout << "The nucleotide rate matrix (" << Codon::nucleotides << ") is:\n"
                  << *this << std::endl;
    }

    void set_mutation_rate(double mu) {
        (*this) *= (mu / mutation_rate);
        mutation_rate = mu;
    }

    double normalized_rate(char from, char to) const { return (*this)(from, to) / mutation_rate; }

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
                    key += Codon::nucleotide(nuc_from);
                    key += Codon::nucleotide(nuc_to);
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
static int dim_gBGC{3};
static int dimensions{4};

class LogMultivariate : public Vector4x1 {
  public:
    explicit LogMultivariate() : Vector4x1(Vector4x1::Zero()) {}

    explicit LogMultivariate(u_long population_size, double mutation_rate_per_generation,
        double generation_time, double gBGC)
        : Vector4x1(Vector4x1::Zero()) {
        set_population_size(population_size);
        set_mutation_rate_per_generation(mutation_rate_per_generation);
        set_generation_time(generation_time);
        set_gBGC(gBGC);
    }

    explicit LogMultivariate(
        double beta, double mutation_rate_per_generation, double generation_time, double gBGC)
        : Vector4x1(Vector4x1::Zero()) {
        set_beta(beta);
        set_mutation_rate_per_generation(mutation_rate_per_generation);
        set_generation_time(generation_time);
        set_gBGC(gBGC);
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
    void set_gBGC(double gBGC) { (*this)(dim_gBGC) = std::log(gBGC); }
    u_long population_size() const {
        return static_cast<u_long>(std::exp((*this)(dim_population_size)));
    }

    double beta() const { return std::exp((*this)(dim_population_size)); }
    double mutation_rate_per_generation() const {
        return std::exp((*this)(dim_mutation_rate_per_generation));
    }
    double generation_time() const { return std::exp((*this)(dim_generation_time)); }
    double gBGC() const { return std::exp((*this)(dim_gBGC)); }

    double log_population_size() const { return (*this)(dim_population_size); }
    double log_mutation_rate_per_generation() const {
        return (*this)(dim_mutation_rate_per_generation);
    }
    double log_generation_time() const { return (*this)(dim_generation_time); }
    double log_gBGC() const { return (*this)(dim_gBGC); }
};

class OrnsteinUhlenbeck {
  public:
    double x{0.0};
    double sigma{0.1};
    double theta{0.9};
    std::default_random_engine &generator;
    std::normal_distribution<double> normal_distrib;

    OrnsteinUhlenbeck(double sigma, double theta, std::default_random_engine &gen)
        : sigma{sigma}, theta{theta}, generator{gen} {
        normal_distrib = std::normal_distribution<double>(0.0, 1.0);
    }

    void Next() { x += sigma * normal_distrib(generator) - theta * x; }
    double GetExpVal() const { return exp(x); }
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
    double ArithmeticgBGC() const { return ArithmeticDim(dim_gBGC); }

    double GeometricPopSize() const { return GeometricDim(dim_population_size); }
    double GeometricMutRate() const { return GeometricDim(dim_mutation_rate_per_generation); }
    double GeometricGenTime() const { return GeometricDim(dim_generation_time); }
    double GeometricgBGC() const { return GeometricDim(dim_gBGC); }

    double HarmonicPopSize() const { return HarmonicDim(dim_population_size); }

    double ArithmeticDim(int dim) const {
        double tot_log{0}, tot_time{0};
        for (std::size_t i = 0; i < times.size(); i++) {
            tot_time += times[i];
            tot_log += times[i] * exp(logmultivariates[i](dim));
        }
        return tot_log / tot_time;
    }

    double GeometricDim(int dim) const {
        double tot_log{0}, tot_time{0};
        for (std::size_t i = 0; i < times.size(); i++) {
            tot_time += times[i];
            tot_log += times[i] * logmultivariates[i](dim);
        }
        return exp(tot_log / tot_time);
    }

    double HarmonicDim(int dim) const {
        double tot_inv{0}, tot_time{0};
        for (std::size_t i = 0; i < times.size(); i++) {
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
        tree.set_tag(node, "Branch_gBGC", d_to_string(GeometricgBGC()));

        tree.set_tag(node, "Branch_arithmetic_population_size", d_to_string(ArithmeticPopSize()));
        tree.set_tag(node, "Branch_arithmetic_generation_time", d_to_string(ArithmeticGenTime()));
        tree.set_tag(node, "Branch_arithmetic_mutation_rate_per_generation",
            d_to_string(ArithmeticMutRate()));
        tree.set_tag(node, "Branch_arithmetic_gBGC", d_to_string(ArithmeticgBGC()));

        tree.set_tag(node, "Branch_harmonic_population_size", d_to_string(HarmonicPopSize()));
    }
};

class CorrelationMatrix : public Matrix4x4 {
  public:
    Matrix4x4 precision_matrix = Matrix4x4::Zero();

    explicit CorrelationMatrix() : Matrix4x4(Matrix4x4::Zero()) {}

    explicit CorrelationMatrix(std::string const &precision_filename, bool fix_pop_size,
        bool fix_mut_rate, bool fix_gen_time, bool fix_gBGC)
        : Matrix4x4(Matrix4x4::Zero()) {
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
        Matrix4x4 covariance_matrix = precision_matrix.inverse();
        assert(covariance_matrix.transpose() == covariance_matrix);
        for (int i = 0; i < dimensions; i++) {
            if (fix_pop_size and i == dim_population_size) { continue; }
            if (fix_mut_rate and i == dim_mutation_rate_per_generation) { continue; }
            if (fix_gen_time and i == dim_generation_time) { continue; }
            if (fix_gBGC and i == dim_gBGC) { continue; }
            for (int j = 0; j < dimensions; j++) {
                if (fix_pop_size and j == dim_population_size) { continue; }
                if (fix_mut_rate and j == dim_mutation_rate_per_generation) { continue; }
                if (fix_gen_time and j == dim_generation_time) { continue; }
                if (fix_gBGC and j == dim_gBGC) { continue; }
                (*this)(i, j) = covariance_matrix(i, j);
            }
        }
        auto part_corr_matrix = precision_matrix;
        auto corr_matrix = covariance_matrix;
        for (int i = 0; i < dimensions; i++) {
            corr_matrix.row(i) /= std::sqrt(covariance_matrix(i, i));
            corr_matrix.col(i) /= std::sqrt(covariance_matrix(i, i));
            part_corr_matrix.row(i) /= -std::sqrt(precision_matrix(i, i));
            part_corr_matrix.col(i) /= std::sqrt(precision_matrix(i, i));
        }
        std::cout << "Covariances (taking fixed effect into account):\n" << *this << std::endl;
        std::cout << "Correlation coefficients (taking fixed effect into account):\n"
                  << corr_matrix << std::endl;
        if (std::abs(corr_matrix.maxCoeff()) > 1.01) {
            std::cerr
                << "The precision matrix is not properly defined, specifically the "
                   "correlation matrix has at least 1 entry greater than one (in absolute value)."
                << std::endl;
            exit(1);
        }
        std::cout << "Precisions:\n" << precision_matrix << std::endl;
        std::cout << "Partial correlation coefficients:\n" << part_corr_matrix << std::endl;
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

class BiasMultivariate : public Vector4x1 {
  public:
    bool step_wise_pop_size{false};
    explicit BiasMultivariate() : Vector4x1(Vector4x1::Zero()) {}

    explicit BiasMultivariate(double bias_population_size, double bias_mutation_rate_per_generation,
        double bias_generation_time, double bias_gBGC, bool step_wise_pop_size)
        : Vector4x1(Vector4x1::Zero()), step_wise_pop_size{step_wise_pop_size} {
        if (step_wise_pop_size and bias_population_size == 0.) {
            std::cerr << "Argument '--step_wise_pop_size' has been set to true, but the argument "
                         "'--bias_pop_size' is set to 0, thus the step_wise has no effect"
                      << std::endl;
        }
        (*this)(dim_population_size) = bias_population_size;
        (*this)(dim_mutation_rate_per_generation) = bias_mutation_rate_per_generation;
        (*this)(dim_generation_time) = bias_generation_time;
        (*this)(dim_gBGC) = bias_gBGC;
    };
};