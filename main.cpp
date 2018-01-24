#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <random>
#include <fstream>

using namespace std;
double epsilon = numeric_limits<double>::epsilon();

#include "Eigen/Dense"

using namespace Eigen;

typedef Matrix<double, 4, 4> Matrix4x4;

#define DOCOPT_HEADER_ONLY
#include "docopt.cpp/docopt.h"

// Definitions:
// Nucleotide: a char from 0 to 3 (included) encoding one of the nucleotide (ATGC).
// triplet: an array of 3 nucleotides (first, second and third position).
// codon: a char from 0 to 63 (included) encoding for a triplet.
// neighboring codons: two codons differing by only 1 mutation.
// neighbor: a 3-tuple containing - in 1st position the codon after mutation,
//                                - in 2nd position the nucleotide before mutation,
//                                - in 3rd position the nucleotide after mutation.

// String of all nucleotides.
string const nucleotides{"ACGT"};

// Function to map a triplet of 3 nucleotides (1st, 2nd and 3rd positions) to the corresponding codon.
// n_1 : the nucleotide in 1st position
// n_2 : the nucleotide in 2nd position
// n_3 : the nucleotide in 3rd position
char triplet_to_codon(char n_1, char n_2, char n_3) {
    assert(0 <= n_1 <= 3);
    assert(0 <= n_2 <= 3);
    assert(0 <= n_3 <= 3);
    int codon{n_1 * 16 + n_2 * 4 + n_3};
    assert(0 <= codon <= 63);
    return char(codon);
}

// Function to build an array mapping each of the 64 codons to their respective triplet.
array<array<char, 3>, 64> build_codon_to_triplet_array() {
    array<array<char, 3>, 64> codon_to_triplet_vec = {0};

    for (char n_1{0}; n_1 < 4; n_1++) {
        // n_1 is the first position of the triplet.
        for (char n_2{0}; n_2 < 4; n_2++) {
            // n_2 is the second position of the triplet.
            for (char n_3{0}; n_3 < 4; n_3++) {
                // n_3 is the third position of the triplet.
                array<char, 3> triplet = {n_1, n_2, n_3};
                codon_to_triplet_vec[triplet_to_codon(n_1, n_2, n_3)] = triplet;
            }
        }
    }
    return codon_to_triplet_vec;
}

// Function to build an array which maps each of the 64 codons to the set of their 9 respective neighbors.
// The set of neighbors contains the 9 codons which differ by only 1 nucleotide (3 positions * 3 possible mutations).
// codon_to_triplet: array mapping each of the 64 codons to their respective triplet (see build_codon_to_triplet_array)
array<array<tuple<char, char, char>, 9>, 64>
build_codon_to_neighbors_array(array<array<char, 3>, 64> &codon_to_triplet) {
    // The array which maps the 64 codons to the set of their 9 respective neighbors.
    array<array<tuple<char, char, char>, 9>, 64> codon_to_neighbors;

    // For each possible codon.
    for (char codon{0}; codon < 64; codon++) {

        // The triplet corresponding to this codon.
        array<char, 3> triplet_from = codon_to_triplet[codon];

        // For each position in the triplet (3 positions).
        for (char position{0}; position < 3; position++) {
            // The original nucleotide which is going to be mutated.
            char n_from = triplet_from[position];

            // For each possible mutation (3 possible mutations).
            for (char mutation{0}; mutation < 3; mutation++) {
                // The mutated nucleotide.
                char n_to = n_from;
                n_to++;
                n_to += mutation;
                n_to %= 4;

                // The mutated triplet.
                array<char, 3> triplet_to = triplet_from;
                triplet_to[position] = n_to;

                // The mutated codon.
                char codon_to = triplet_to_codon(triplet_to[0], triplet_to[1], triplet_to[2]);

                // Assign the neighbor to the array of neighbors.
                codon_to_neighbors[codon][3 * position + mutation] = make_tuple(codon_to, n_from, n_to);
            }
        }
    }
    return codon_to_neighbors;
}

// Array mapping each of the 64 codons to their respective triplet.
auto codon_to_triplet_array = build_codon_to_triplet_array();

// The set of neighbors contains the 9 codons which differ by only 1 nucleotide (3 positions * 3 possible mutations).
// Array which maps each of the 64 codons to the set of their 9 respective neighbors.
auto const codon_to_neighbors_array = build_codon_to_neighbors_array(codon_to_triplet_array);

// String containing the 21 amino-acids (20 + stop) for translations from index to amino-acid character
string const amino_acids{"ACDEFGHIKLMNPQRSTVWYX"};

// Map from codons to amino-acids
map<string, char> triplet_str_to_aa_char = {{"GCA", 'A'},
                                            {"GAA", 'E'},
                                            {"ACT", 'T'},
                                            {"CAT", 'H'},
                                            {"ACG", 'T'},
                                            {"GGT", 'G'},
                                            {"GCG", 'A'},
                                            {"GAG", 'E'},
                                            {"CGC", 'R'},
                                            {"TGA", 'X'},
                                            {"CGG", 'R'},
                                            {"GCC", 'A'},
                                            {"TGC", 'C'},
                                            {"GAC", 'D'},
                                            {"CAA", 'Q'},
                                            {"CGT", 'R'},
                                            {"GAT", 'D'},
                                            {"TCA", 'S'},
                                            {"CAC", 'H'},
                                            {"ATC", 'I'},
                                            {"CGA", 'R'},
                                            {"ATA", 'I'},
                                            {"GCT", 'A'},
                                            {"CAG", 'Q'},
                                            {"TGG", 'W'},
                                            {"GGC", 'G'},
                                            {"TTC", 'F'},
                                            {"CCA", 'P'},
                                            {"ACC", 'T'},
                                            {"TAC", 'Y'},
                                            {"GTG", 'V'},
                                            {"AAC", 'N'},
                                            {"AAG", 'K'},
                                            {"CCT", 'P'},
                                            {"TCC", 'S'},
                                            {"CCC", 'P'},
                                            {"CTC", 'L'},
                                            {"GTT", 'V'},
                                            {"AGC", 'S'},
                                            {"ATT", 'I'},
                                            {"ACA", 'T'},
                                            {"TTG", 'L'},
                                            {"GTC", 'V'},
                                            {"AGT", 'S'},
                                            {"CTG", 'L'},
                                            {"TCG", 'S'},
                                            {"TAT", 'Y'},
                                            {"TTT", 'F'},
                                            {"AAT", 'N'},
                                            {"CCG", 'P'},
                                            {"TTA", 'L'},
                                            {"TGT", 'C'},
                                            {"GGA", 'G'},
                                            {"CTA", 'L'},
                                            {"AAA", 'K'},
                                            {"GGG", 'G'},
                                            {"ATG", 'M'},
                                            {"GTA", 'V'},
                                            {"TCT", 'S'},
                                            {"AGA", 'R'},
                                            {"TAA", 'X'},
                                            {"TAG", 'X'},
                                            {"AGG", 'R'},
                                            {"CTT", 'L'}};

// Function to map a particular amino-acid to its index in the string containing all the amino-acids
// Equivalent to the inverse function of translating index to amino-acid character
// aa_char: character of an amino-acid
// amino_acids: string containing all the amino-acids
char aa_char_to_aa(const char aa_char, const string &amino_acids) {
    return char(amino_acids.find(aa_char));
}

// Function to build an array mapping each of the 64 codons to their respective amino-acid.
// triplet_str_to_aa_char: map from codons to amino-acids
// amino_acids: string containing all the amino-acids
// codon_to_triplet: array mapping each of the 64 codons to their respective triplet (see build_codon_to_triplet_array)
array<char, 64> build_codon_to_aa_array(map<string, char> &triplet_str_to_aa_char,
                                        const string &amino_acids,
                                        array<array<char, 3>, 64> &codon_to_triplet) {
    // Array mapping each of the 64 codons to their respective amino-acid.
    array<char, 64> codon_to_aa = {0};

    // For each codon
    for (char codon{0}; codon < 64; codon++) {

        // Triplet corresponding to the codon
        array<char, 3> triplet = codon_to_triplet_array[codon];

        // String containing the translated triplet (e.g "ATG")
        string triplet_str(3, ' ');
        for (char position{0}; position < 3; position++) {
            triplet_str[position] = nucleotides[triplet[position]];
        }

        // Amino-acid corresponding to the translated string
        codon_to_aa[codon] = aa_char_to_aa(triplet_str_to_aa_char[triplet_str], amino_acids);
    }
    return codon_to_aa;
}

// Array mapping each of the 64 codons to their respective amino-acid.
array<char, 64> const codon_to_aa_array = build_codon_to_aa_array(triplet_str_to_aa_char, amino_acids,
                                                                  codon_to_triplet_array);

// Random generator engine with seed 0.
default_random_engine re(0);

struct Substitution {
    unsigned site{0};
    char codon_from{0};
    char codon_to{0};
};

// Class representing DNA sequence.
class Sequence_dna {
private:
    // The number of sites in the sequence (each position is a codon, thus the DNA sequence is 3 times greater).
    const unsigned nbr_sites;

    // The sequence of codons.
    vector<char> codon_seq;

    // The matrix of mutation rates between nucleotides.
    Matrix4x4 mutation_rate_matrix;


    // The fitness profiles of amino-acids.
    vector<array<double, 20>> aa_fitness_profiles;
public:
    // Constructor of Sequence_dna.
    // size: the size of the DNA sequence.
    explicit Sequence_dna(const unsigned size) : nbr_sites{size}, codon_seq(size, 0), aa_fitness_profiles{0} {
        mutation_rate_matrix.Zero();
    }

    // Method translating the codon sequence to amino-acid sequence.
    string translate() const {
        string protein(nbr_sites, ' ');

        // For each site
        for (unsigned site{0}; site < nbr_sites; site++) {
            // Use the codon_to_aa_array to translate site to amino-acid.
            protein[site] = amino_acids[codon_to_aa_array[codon_seq[site]]];
        }

        return protein; // return the amino-acid sequence as a string.
    }

    // Method computing the next substitution event to occur, and the time for it to happen.
    // This method takes a time as input. If there is enough time given, the substitution event is computed.
    // This method also returns the time (given as input) decremented by the time needed for the substitution event to occur.
    // If there is not enough time given, no substitution event is computed and the method returns 0.
    // time_left: The time available for a substitution event to occur.
    double next_substitution(double &time_left, vector<Substitution> &substitutions) {

        // Number of possible substitutions is 9 times the number of sites (3 substitutions for each 3 possible positions).
        unsigned nbr_substitutions{9 * nbr_sites};

        // Vector of substitution rates.
        vector<double> substitution_rates(nbr_substitutions, 0);

        // Sum of substitution rates.
        double total_substitution_rates{0.};

        // For all site of the sequence.
        for (unsigned site{0}; site < nbr_sites; site++) {
            // Codon original before substitution.

            char codon_from = codon_seq[site];

            // Array of neighbors of the original codon (codons differing by only 1 mutation).
            array<tuple<char, char, char>, 9> neighbors = codon_to_neighbors_array[codon_from];

            // For all possible neighbors.
            for (char neighbor{0}; neighbor < 9; neighbor++) {

                // Codon after mutation, Nucleotide original and Nucleotide after mutation.
                char codon_to{0}, n_from{0}, n_to{0};
                tie(codon_to, n_from, n_to) = neighbors[neighbor];

                // Assign the substitution rate given by the method substitution rate.
                // Rate of fixation initialized to 1 (neutral mutation)
                double rate_fixation{1.};

                // If the mutated amino-acid is a stop codon, the rate of fixation is 0.
                // Else, if the mutated and original amino-acids are non-synonymous, we compute the rate of fixation.
                // Note that, if the mutated and original amino-acids are synonymous, the rate of fixation is 1.
                if (codon_to_aa_array[codon_to] == 20) {
                    rate_fixation = 0.;
                } else if (codon_to_aa_array[codon_from] != codon_to_aa_array[codon_to]) {
                    // Selective strength between the mutated and original amino-acids.
                    double s{0.};
                    s = aa_fitness_profiles[site][codon_to_aa_array[codon_to]];
                    s -= aa_fitness_profiles[site][codon_to_aa_array[codon_from]];
                    // If the selective strength is 0, the rate of fixation is neutral.
                    // Else, the rate of fixation is computed using population genetic formulas (Kimura).
                    if (fabs(s) <= epsilon) {
                        rate_fixation = 1;
                    } else {
                        rate_fixation = s / (1 - exp(-s));
                    }
                }

                // The substitution rate is the mutation rate multiplied by the rate of fixation.
                substitution_rates[9 * site + neighbor] = rate_fixation * mutation_rate_matrix(n_from, n_to);
                // Increment the sum of substitution rates
                total_substitution_rates += substitution_rates[9 * site + neighbor];

            }
        }

        // Decrement the time by inverse of the sum of substitution rates.
        time_left -= 1. / total_substitution_rates;

        // Substitute the sequence if the time is positive, else there is no substitution but the time left is set to 0.
        if (time_left > 0. and total_substitution_rates != 0.) {

            // Uniform random generator between 0 and the total sum of substitution rates.
            uniform_real_distribution<double> unif(0, total_substitution_rates);
            double random_cumulative_substitution_rates = unif(re);
            double cumulative_substitution_rates{0.};

            unsigned index = 0;
            for (unsigned t{0}; t < nbr_substitutions; t++) {
                // Iterate through the cumulative substitution rates and break the loop when it is greater than the random cumulative substitution rates
                cumulative_substitution_rates += substitution_rates[t];
                if (random_cumulative_substitution_rates < cumulative_substitution_rates) {
                    index = t;
                    break;
                }
            }

            // Array of neighbors of the original codon (codons differing by only 1 mutation).
            array<tuple<char, char, char>, 9> neighbors = codon_to_neighbors_array[codon_seq[index / 9]];
            // Codon after mutation, Nucleotide original and Nucleotide after mutation.
            char codon_to{0}, n_from{0}, n_to{0};

            tie(codon_to, n_from, n_to) = neighbors[index % 9];
            unsigned site = index / 9;

            Substitution substitution = {site, codon_seq[site], codon_to};
            substitutions.push_back(substitution);
            codon_seq[site] = codon_to;

        } else if (time_left < 0.) {
            time_left = 0.;
        }
        return time_left;
    }


    // Method computing all substitution event occuring during a given time-frame.
    // t: time during which substitution events occur (typically branch length).
    void run_substitutions(double t, vector<Substitution> &substitutions) {
        while (t > 0) {
            t = next_substitution(t, substitutions);
        }
    }


    // Set the the DNA sequence to the mutation-selection equilibrium.
    void at_equilibrium() {

        // Compute the kernel of the mutation-rate matrix.
        // This kernel is a vector of the nucleotides frequencies (not normalized to 1) at equilibrium.
        Matrix<double, 4, Dynamic> kernel = mutation_rate_matrix.transpose().fullPivLu().kernel();
        Matrix<double, 4, 1> nuc_frequencies;

        if (kernel.cols() > 1) {
            cerr << "The kernel has " << kernel.cols() << " dimensions, this is weird ! " << endl;
            uniform_int_distribution<unsigned> unif_int(0, unsigned(kernel.cols()) - 1);
            unsigned chosen_row = unif_int(re);
            nuc_frequencies = kernel.col(chosen_row);

        } else {
            nuc_frequencies = kernel.col(0);
        }

        array<double, 64> codon_frequencies{0};

        // For each site of the sequence.
        for (unsigned site{0}; site < nbr_sites; site++) {

            // Initialize the total sum of equilibrium frequency at 0.
            double total_frequencies{0.};

            // For each site of the vector of the site frequencies.
            for (char codon{0}; codon < 64; codon++) {
                double codon_freq{1.};

                // Translate the site to a triplet of DNA nucleotides
                array<char, 3> triplet = codon_to_triplet_array[codon];
                for (char position{0}; position < 3; position++) {
                    codon_freq *= nuc_frequencies(triplet[position]);
                }

                if (codon_to_aa_array[codon] != 20) {
                    codon_frequencies[codon] = exp(aa_fitness_profiles[site][codon_to_aa_array[codon]]);
                } else {
                    codon_frequencies[codon] = 0.;
                }


                // Increment the total sum of equilibrium frequencies.
                total_frequencies += codon_frequencies[codon];
            }

            // Uniform random generator between 0 and the total sum of equilibrium frequencies.
            uniform_real_distribution<double> unif(0, total_frequencies);
            double random_cumulative_frequencies = unif(re);
            double cumulative_frequencies{0.};

            char index{0};
            for (char m{0}; m < 64; m++) {
                // Iterate through the cumulative frequencies and break the loop when it is greater than the random cumulative frequencies.
                cumulative_frequencies += codon_frequencies[m];
                if (random_cumulative_frequencies < cumulative_frequencies) {
                    index = m;
                    break;
                }
            }

            // Substitute the site with the substitution given by the loop break.
            codon_seq[site] = index;
        }
    }

    // Set attribute method for the codon sequence.
    void set_parameters(Sequence_dna const& sequence_dna) {
        codon_seq = sequence_dna.codon_seq;
        mutation_rate_matrix = sequence_dna.mutation_rate_matrix;
        aa_fitness_profiles = sequence_dna.aa_fitness_profiles;
    }

    // Set attribute method for the mutation rate matrix.
    void set_mutation_rate(Matrix4x4 const &mutation_rate) {
        mutation_rate_matrix = mutation_rate;
    }

    // Set attribute method for the amino-acid fitness profil.
    void set_fitness_profiles(vector<array<double, 20>> const &fitness_profiles) {
        assert(nbr_sites == fitness_profiles.size());
        aa_fitness_profiles = fitness_profiles;
    }

    // Method returning the DNA string corresponding to the codon sequence.
    string get_dna_str() const {

        // The DNA string is 3 times larger than the codon sequence.
        string dna_str(nbr_sites * 3, ' ');

        // For each site of the sequence.
        for (unsigned site{0}; site < nbr_sites; site++) {

            // Assert there is no stop in the sequence.
            assert(codon_to_aa_array[codon_seq[site]] != 20);

            // Translate the site to a triplet of DNA nucleotides
            array<char, 3> triplet = codon_to_triplet_array[codon_seq[site]];
            for (char position{0}; position < 3; position++) {
                dna_str[3 * site + position] = nucleotides[triplet[position]];
            }
        }
        return dna_str; // return the DNA sequence as a string.
    }
};


// Class representing nodes of a tree.
class Node {
private:
    string name;  // The species name of the node.
    double length;  // The length of the branch attached to the node (ascending).
    string newick;  // The newick tree descending from the node.
    vector<Node> children;  // Vector of direct children (first order, no grand-children).
    Sequence_dna sequence_dna;  // The DNA sequence attached to the node.
    vector<Substitution> substitutions;  // The number of substitutions in the branch attached to the node.

public:

    // Constructor
    Node(string name, const string &len, string newick, Sequence_dna seq) :
            name{move(name)}, length{stod(len)}, newick{move(newick)}, sequence_dna{move(seq)}, substitutions{0} {

        // Parse the newick tree descending of the node.
        parse_newick();
    }

    Node(string newick, const unsigned nbr_sites) :
            name{"Root"}, length{0.}, newick{move(newick)}, sequence_dna(nbr_sites), substitutions{0} {

        // Parse the newick tree descending of the node.
        parse_newick();
    }

    // Is true if the node don't have children.
    bool is_leaf() const {
        return newick.length() == 0;
    }

    // Add a node as the vector of children.
    void add_child(const Node &node) {
        children.push_back(node);
    }

    // Recursively iterate through the subtree and count the number of nodes.
    double tot_length() {
        double tot_length = length;

        if (!is_leaf()) {
            // Else, if the node is internal, iterate through the direct children.
            for (auto &child : children) {
                tot_length += child.tot_length();
            }
        }
        return tot_length;
    }

    // Recursively iterate through the subtree and count the number of nodes.
    unsigned nbr_nodes() {
        unsigned nbr_nodes = 1;

        if (!is_leaf()) {
            // Else, if the node is internal, iterate through the direct children.
            for (auto &child : children) {
                nbr_nodes += child.nbr_nodes();
            }
        }
        return nbr_nodes;
    }

    // Recursively iterate through the subtree and count the number of leaves.
    unsigned nbr_leaves() {
        unsigned nbr_leaves = 0;

        if (is_leaf()) {
            // If the node is a leaf, return 1.
            nbr_leaves = 1;
        } else {
            // Else, if the node is internal, iterate through the direct children.
            for (auto &child : children) {
                nbr_leaves += child.nbr_leaves();
            }
        }
        return nbr_leaves;
    }

    // Recursively iterate through the subtree and count the number of substitutions.
    unsigned long nbr_substitutions() {
        unsigned long nbr_substitutions = substitutions.size();

        if (!is_leaf()) {
            // Else, if the node is internal, iterate through the direct children.
            for (auto &child : children) {
                nbr_substitutions += child.nbr_substitutions();
            }
        }
        return nbr_substitutions;
    }

    // Recursively iterate through the subtree.
    void traverse(string & output_filename) {
        // Substitutions of the DNA sequence is generated.
        sequence_dna.run_substitutions(length, substitutions);

        if (is_leaf()) {
            // If the node is a leaf, output the DNA sequence and name.
            ofstream output_file;
            output_file.open(output_filename, ios_base::app);
            output_file << name << " " << sequence_dna.get_dna_str() << endl;
            output_file.close();
        } else {
            // If the node is internal, iterate through the direct children.
            for (auto &child : children) {
                child.sequence_dna.set_parameters(sequence_dna);
                child.traverse(output_filename);
            }
        }
    }

    // Set method for the parameters of evolution of the sequence
    void set_evolution_parameters(Matrix4x4 const &mutation_rate, vector<array<double, 20>> const &fitness_profiles) {
        sequence_dna.set_mutation_rate(mutation_rate);
        sequence_dna.set_fitness_profiles(fitness_profiles);
    }

    // Set the the DNA sequence to the mutation-selection equilibrium.
    void at_equilibrium() {
        sequence_dna.at_equilibrium();
    }

    void parse_newick() {

        if (!is_leaf()) {
            // The size of the string of newick tree.
            size_t max_position{newick.size()};
            // The current position in the string of the newick tree.
            size_t position{0};

            // While the current position is lower than the size of the string, their is at least one node to parse.
            while (position < max_position) {
                // 'subtree' is the left hand side of the node name, it can be a subtree or nothing if the node is a leaf.
                string subtree{};
                if (newick[position] == '(') {

                    size_t postpoint{position};
                    unsigned nbr_open{1};

                    for (size_t i{position + 1}; i < max_position; i++) {
                        if (nbr_open == 0) {
                            postpoint = i;
                            break;
                        } else if (newick[i] == '(') {
                            nbr_open++;
                        } else if (newick[i] == ')') {
                            nbr_open--;
                        };
                    }
                    subtree = newick.substr(position + 1, postpoint - position - 2);
                    position = postpoint;
                }

                // 'name_suffix' contains the name of the node and the branch length.
                string name_suffix{};

                size_t next_sep = newick.substr(position).find(',');
                if (next_sep == string::npos) {
                    name_suffix = newick.substr(position);
                    position = max_position;
                } else {
                    name_suffix = newick.substr(position, next_sep);
                    position = position + next_sep + 1;
                }

                // 'length' contains the name of the node.
                string length{};
                // 'name' contains the branch length of the node.
                string name{};
                size_t ddot = name_suffix.rfind(':');
                if (ddot != string::npos) {
                    length = name_suffix.substr(ddot + 1);
                    name = name_suffix.substr(0, ddot);
                } else {
                    name = name_suffix;
                    length = "0";
                }

                // New node from 'subtree', 'name' and 'length' using the DNA sequence of this node.
                add_child(Node(name, length, subtree, sequence_dna));
            }
        }
    }
};


string open_newick(string const &file_name) {
    ifstream input_stream(file_name);
    if (!input_stream) cerr << "Can't open newick file!";

    string line;
    getline(input_stream, line);

    return line;
}

vector<array<double, 20>> open_preferences(string const &file_name) {
    vector<array<double, 20>> fitness_profiles{0};

    ifstream input_stream(file_name);
    if (!input_stream) cerr << "Can't open preferences file!";

    string line;

    // skip the header of the file
    getline(input_stream, line);

    while (getline(input_stream, line)) {

        array<double, 20> fitness_profil{0};
        string word;
        istringstream line_stream(line);
        unsigned counter{0};

        while (getline(line_stream, word, ' ')) {
            if (counter > 2) {
                fitness_profil[counter - 3] = log(stod(word));
            }
            counter++;
        }

        fitness_profiles.push_back(fitness_profil);
    }
    return fitness_profiles;
}

static const char USAGE[] =
R"(
Usage:
      SimuEvol [--protein=<name>] [--output=<filename>]
      SimuEvol --help
      SimuEvol --version

Options:
-h --help              show this help message and exit
--version              show version and exit
--protein=<name>       specify input protein name [default: gal4]
--output=<filename>    specify input protein name [default: protein.ali]
)";

int main(int argc, char* argv[]) {

    auto args = docopt::docopt(USAGE,
                             { argv + 1, argv + argc },
                             true,               // show help if requested
                             "SimuEvol 0.1a");  // version string

    string protein{"gal4"};
    if (args["--protein"]){
        protein = args["--protein"].asString();
    }

    string output{"../data/"};
    if (args["--output"]){
        output += args["--output"].asString();
    } else {
        output += protein + ".ali";
    }

    vector<array<double, 20>> fitness_profiles = open_preferences("../data/" + protein + ".txt");
    auto nbr_sites = static_cast<unsigned>(fitness_profiles.size());

    string newick_tree = open_newick("../data/" + protein + ".newick");
    Node root(newick_tree, nbr_sites);

    Matrix4x4 mutation_rate;
    double mu = 5 * pow(10, -1);
    double lambda = 3;
    mutation_rate << 0, 1, 1, lambda,
            lambda, 0, 1, lambda,
            lambda, 1, 0, lambda,
            lambda, 1, 1, 0;
    mutation_rate *= mu;
    mutation_rate -= mutation_rate.rowwise().sum().asDiagonal();

    // If the node is a leaf, output the DNA sequence and name.
    ofstream output_file;
    output_file.open(output);
    output_file << root.nbr_leaves() << " " << nbr_sites * 3 << endl;
    output_file.close();

    cout << "The tree contains " << root.nbr_nodes() << " nodes for " << root.nbr_leaves() << " species at the tips." << endl;
    cout << "The tree has a total branch length of " << root.tot_length() << "." << endl;
    cout << "The DNA sequence is " << nbr_sites * 3 << " base pairs long." << endl;

    root.set_evolution_parameters(mutation_rate, fitness_profiles);
    root.at_equilibrium();
    root.traverse(output);
    cout << "The simulation mapped " << root.nbr_substitutions() << " substitutions along the tree." << endl;
    return 0;
}
