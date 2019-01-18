#pragma once

#include <cassert>
#include "Eigen/Dense"

typedef Eigen::Matrix<double, 4, 4> Matrix4x4;
typedef Eigen::Matrix<double, 4, 1> Vector4x1;

// Compute the kernel of the mutation-rate matrix.
// This kernel is a vector of the nucleotides frequencies (not normalized to 1) at equilibrium.
Vector4x1 equilibrium_frequencies(Matrix4x4 const &mutation_matrix) {
    Eigen::Matrix<double, 4, Eigen::Dynamic> kernel =
        mutation_matrix.transpose().fullPivLu().kernel();
    Vector4x1 nuc_frequencies;

    assert(kernel.cols() == 1);
    nuc_frequencies = kernel.col(0);
    nuc_frequencies /= nuc_frequencies.sum();

    return nuc_frequencies;
}


void normalize_nucmatrix(Matrix4x4 &mutation_matrix) {
    for (int diag{0}; diag < 4; diag++) { mutation_matrix(diag, diag) = 0.0; }
    mutation_matrix -= mutation_matrix.rowwise().sum().asDiagonal();

    double events{0.};
    Vector4x1 nuc_frequencies = equilibrium_frequencies(mutation_matrix);
    std::cout << "The equilibrium nucleotides frequencies (" << Codon::nucleotides << ") are:\n"
              << nuc_frequencies << std::endl;


    bool reversible = true;
    for (int row{0}; row < 4; row++) {
        for (int col{0}; col < 4; col++) {
            double diff = nuc_frequencies(row) * mutation_matrix(row, col) -
                          nuc_frequencies(col) * mutation_matrix(col, row);
            if (std::abs(diff) > 1e-6) { reversible = false; }
        }
    }
    if (reversible) {
        std::cout << "The nucleotide rate matrix is time-reversible." << std::endl;
    } else {
        std::cerr << "The nucleotide rate matrix is not time-reversible." << std::endl;
    }

    for (int diag{0}; diag < 4; diag++) {
        events -= nuc_frequencies(diag) * mutation_matrix(diag, diag);
    }
    mutation_matrix /= events;
}

Matrix4x4 input_matrix(std::string const &input_filename) {
    Matrix4x4 mutation_matrix = Matrix4x4::Ones();
    std::ifstream input_stream(input_filename);

    if (input_filename.empty()) {
        std::cerr << "No rate matrix file provided, use the default nucleotide rate matrix." << std::endl;
        return mutation_matrix;
    }

    if (!input_stream) {
        std::cerr << "Rate matrix file " << input_filename
                  << "doesn't exist, use the default rate matrix nucleotide instead." << std::endl;
        return mutation_matrix;
    }

    std::string line;
    getline(input_stream, line);
    line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
    assert(line == Codon::nucleotides);

    int row = 0;
    while (getline(input_stream, line)) {
        std::stringstream line_stream(line);
        std::string word;
        int col = 0;

        while (getline(line_stream, word, '\t')) {
            if (word != "-") { mutation_matrix(row, col) = stod(word); }
            col++;
        }
        row++;
    }

    normalize_nucmatrix(mutation_matrix);

    std::cout << "The mutation transition matrix (" << Codon::nucleotides << ") is:\n"
              << mutation_matrix << std::endl;
    return mutation_matrix;
}