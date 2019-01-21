#pragma once

#include <cassert>
#include <random>
#include "Eigen/Dense"
typedef Eigen::Matrix<double, 4, 4> Matrix4x4;
typedef Eigen::Matrix<double, 4, 1> Vector4x1;


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
                if (nuc_from != nuc_to) {
                    rates[nuc_to] = (*this)(nuc_from, nuc_to);
                }
            }
            mutation_distr[nuc_from] = std::discrete_distribution<char>(rates.begin(), rates.end());
        }
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
};