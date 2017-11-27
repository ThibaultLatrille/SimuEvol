#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <random>
#include "Eigen/Dense"

using namespace Eigen;
typedef Matrix<double, 64, 64> Matrix64x64;
typedef Matrix<double, 64, 1> Vector64;

using namespace std;

// Definitions:
// Nucleotide: a char from 0 to 3 (included) encoding one of the nucleotide (ATGC).
// triplet: an array of 3 nucleotides (first, second and third position).
// codon: a char from 0 to 63 (included) encoding for a triplet.
// neighboring codons: two codons differing by only 1 mutation.
// neighbor: a 3-tuple containing - in 1st position the codon after mutation,
//                                - in 2nd position the nucleotide before mutation,
//                                - in 3rd position the nucleotide after mutation.

// String of all nucleotides.
string const nucleotides{"ATGC"};

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
    array<array<char, 3>, 64> codon_to_triplet_vec;

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

            // For each possible mutation (3 possible).
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
// The set of neighbors contains the 9 codons which differ by only 1 nucleotide (3 positions * 3 possible mutations).
auto codon_to_triplet_array = build_codon_to_triplet_array();

// Array which maps each of the 64 codons to the set of their 9 respective neighbors.
auto const codon_to_neighbors_array = build_codon_to_neighbors_array(codon_to_triplet_array);

// String containing the 21 amino-acids (20 + stop) for translations from index to amino-acid character
string const amino_acids{"ARNDCQEGHILKMFPSTWYVX"};

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
    array<char, 64> codon_to_aa;

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

// Class representing DNA sequence.
class Sequence_dna {
private:
    // The sequence length (each position is a codon, thus the DNA sequence is 3 times greater).
    const unsigned length;

    // The sequence of codons.
    vector<char> codon_seq;

    // The matrix of mutation rates between nucleotides.
    Matrix<double, 4, 4> mutation_rate_matrix;

    // The fitness profil of amino-acids.
    vector<double> aa_fitness_profil;

public:
    // Constructor of Sequence_dna.
    // len: the size of the DNA sequence.
    Sequence_dna(const unsigned len) : length{len}, codon_seq(len, 0), aa_fitness_profil(21, 0.) {

        // Uniform random generator between 0 (included) and 63 (included) for assigning codons.
        uniform_int_distribution<char> unif_int(0, 63);

        // Assign each codon randomly.
        for (char codon{0}; codon < length; codon++) {
            codon_seq[codon] = unif_int(re);
        }

        // Uniform random generator between 0. and 1. for assigning mutation rates.
        uniform_real_distribution<double> unif_real(0., 1);

        // For each original nucleotide.
        for (char from{0}; from < 4; from++) {
            double total{0};
            // For each mutated nucleotide.
            for (char to{0}; to < 4; to++) { 
                // Assign randomly mutation rates if their is really a mutation (not from G to G for example).
                if (from != to) {
                    mutation_rate_matrix(from, to) = unif_real(re);
                    total -= mutation_rate_matrix(from, to);
                }
            }
            mutation_rate_matrix(from, from) = total;
        }

        // Assign randomly the fitness of amino-acids.
        for (char aa_index{0}; aa_index < 20; aa_index++) {
            aa_fitness_profil[aa_index] = unif_real(re);
        }

    }

    // Method translating the codon sequence to amino-acid sequence.
    string translate() const {
        string protein(length, ' ');

        // For each codon
        for (unsigned codon{0}; codon < length; codon++) {
            // Use the codon_to_aa_array to translate codon to amino-acid.
            protein[codon] = amino_acids[codon_to_aa_array[codon_seq[codon]]];
        }

        return protein; // return the amino-acid sequence as a string.
    }

    // Method returning the substitution rate (mutation rate * rate of fixation) between an original an a mutated codon.
    // The original and mutated nucleotides are also necessary !
    // codon_from: codon original.
    // codon_to: codon mutated.
    // n_from: nucleotide original.
    // n_to: nucleotide mutated.
    double substitution_rate(char codon_from, char codon_to, char n_from, char n_to) {
        // Rate of fixation initialized to 1 (neutral mutation)
        double rate_fixation{1.};

        // If the mutated amino-acid is a stop codon, the rate of fixation is 0.
        // Else, if the mutated and original amino-acids are non-synonymous, we compute the rate of fixation.
        // Note that, if the mutated and original amino-acids are synonymous, the rate of fixation is 1.
        if (codon_to_aa_array[codon_to] == 20) {
            rate_fixation = 0.;
        } else if (codon_from != codon_to) {
            // Selective strength between the mutated and original amino-acids.
            double s{0};
            s = aa_fitness_profil[codon_to_aa_array[codon_to]] - aa_fitness_profil[codon_to_aa_array[codon_from]];

            // If the selective strength is 0, the rate of fixation is neutral.
            // Else, the rate of fixation is computed using population genetic formulas (Kimura).
            if (s == .0) {
                rate_fixation = 1;
            } else {
                rate_fixation = s / (1 - exp(-s));
            }
        }

        // The substitution rate is the mutation rate multiplied by the rate of fixation.
        return rate_fixation * mutation_rate_matrix(n_from, n_to);
    }

    // Method computing the next substitution event to occur, and the time for it to happen.
    // This method takes a time as input. If there is enough time given, the substitution event is computed.
    // This method also returns the time (given as input) decremented by the time needed for the substitution event to occur.
    // If there is not enough time given, no substitution event is computed and the method returns 0.
    // time_left: The time available for a substitution event to occur.
    double next_substitution(double &time_left) {

        // Number of possible substitutions is 9 times the sequence length (3 substitutions for each 3 possible positions).
        unsigned nbr_substitutions{9 * length};

        // Vector of substitution rates.
        vector<double> substitution_rates(nbr_substitutions, 0);

        // Sum of substitution rates.
        double total_substitution_rates{0.};

        // For all codon of the sequence.
        for (unsigned codon_index{0}; codon_index < length; codon_index++) {
            // Codon original before substitution.
            char codon_from = codon_seq[codon_index];

            // Array of neighbors of the original codon (codons differing by only 1 mutation).
            array<tuple<char, char, char>, 9> neighbors = codon_to_neighbors_array[codon_from];

            // For all possible neighbors.
            for (char neighbor{0}; neighbor < 9; neighbor++) {

                // Codon after mutation, Nucleotide original and Nucleotide after mutation.
                char codon_to{0}, n_from{0}, n_to{0};
                tie(codon_to, n_from, n_to) = neighbors[neighbor];

                // Assign the substitution rate given by the method substitution rate.
                substitution_rates[9 * codon_index + neighbor] = substitution_rate(codon_from, codon_to, n_from, n_to);
                // Increment the sum of substitution rates
                total_substitution_rates += substitution_rates[9 * codon_index + neighbor];

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
            codon_seq[index / 9] = codon_to;

        } else if (time_left < 0.) {
            time_left = 0.;
        }
        return time_left;
    }


    // Method computing all substitution event occuring during a given time-frame.
    // t: time during which substitution events occur (typically branch length).
    void run_substitutions(double t) {
        while (t > 0) {
            t = next_substitution(t);
        }
    }


    // Set the the DNA sequence to the mutation-selection equilibrium.
    void at_equilibrium() {

        // Matrix of substitution rates between codons.
        Matrix64x64 codon_matrix(Matrix64x64::Zero());

        // For each possible codon before mutation.
        for (char codon_from{0}; codon_from < 64; codon_from++) {

            // Array of neighbors of the codon (codons differing by only 1 mutation).
            array<tuple<char, char, char>, 9> neighbors = codon_to_neighbors_array[codon_from];

            // Total rate of substitution from codon to all neighbors.
            double total{0};

            // For all neighbors.
            for (unsigned neighbor{0}; neighbor < 9; neighbor++) {

                // Codon after mutation, Nucleotide original and Nucleotide after mutation.
                char codon_to{0}, n_from{0}, n_to{0};
                tie(codon_to, n_from, n_to) = neighbors[neighbor];

                // Assign the substitution rate given by the method substitution rate.
                // Note that the matrix is transposed !
                codon_matrix(codon_to, codon_from) = substitution_rate(codon_from, codon_to, n_from, n_to);

                // Increment the total rate of substitutions.
                total -= codon_matrix(codon_to, codon_from);
            }
            codon_matrix(codon_from, codon_from) = total;
        }

        // Compute the kernel of the substitution rate matrix.
        // This kernel is a vector of the codon frequencies (not normalized to 1) at equilibrium.
        FullPivLU<Matrix64x64> lu_decomp(codon_matrix);
        Vector64 codon_frequencies(lu_decomp.kernel());

        // For each codon of the sequence.
        for (unsigned codon{0}; codon < length; codon++) {

            // Initialize the total sum of equilibrium frequency at 0.
            double total_frequencies{0.};

            // For each codon frequency of the vector of the codon frequencies.
            for (char codon_frequency{0}; codon_frequency < 64; codon_frequency++) {

                // Increment the total sum of equilibrium frequencies.
                total_frequencies += codon_frequencies(codon_frequency);
            }

            // Uniform random generator between 0 and the total sum of equilibrium frequencies.
            uniform_real_distribution<double> unif(0, total_frequencies);
            double random_cumulative_frequencies = unif(re);
            double cumulative_frequencies{0.};

            char index{0};
            for (char m{0}; m < 64; m++) {
                // Iterate through the cumulative frequencies and break the loop when it is greater than the random cumulative frequencies.
                cumulative_frequencies += codon_frequencies(m);
                if (random_cumulative_frequencies < cumulative_frequencies) {
                    index = m;
                    break;
                }
            }

            // Substitute the codon with the substitution given by the loop break.
            codon_seq[codon] = index;
        }
    }

    // Get attribute method for the codon sequence.
    vector<char> get_codon_seq() const { return codon_seq; }

    // Set attribute method for the codon sequence.
    void set_codon_seq(vector<char> const codon) {
        codon_seq = codon;
    }

    // Method returning the DNA string corresponding to the codon sequence.
    string get_dna_str() const {

        // The DNA string is 3 times larger than the codon sequence.
        string dna_str(length * 3, ' ');

        // For each codon of the sequence.
        for (unsigned codon{0}; codon < length; codon++) {

            // Assert there is no stop in the sequence.
            assert(codon_to_aa_array[codon_seq[codon]] != 20);

            // Translate the codon to a triplet of DNA nucleotides
            array<char, 3> triplet = codon_to_triplet_array[codon_seq[codon]];
            for (char position{0}; position < 3; position++) {
                dna_str[3 * codon + position] = nucleotides[triplet[position]];
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

public:

    // Constructor
    Node(const string &name, const string &len, const string &newick, const Sequence_dna &seq) :
            name{name}, length{stod(len)}, newick{newick}, sequence_dna{seq} {

        // Parse the newick tree descending of the node.
        parse_newick();
    }

    Node(const string &newick, const unsigned length) :
            name{"Root"}, length{0.}, newick{newick}, sequence_dna(length) {

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

    // Recursively iterate through the subtree.
    void traverse() {
        // Substitutions of the DNA sequence is generated.
        sequence_dna.run_substitutions(length);

        if (is_leaf()) {
            // If the node is a leaf, output the DNA sequence and name.
            cout << ">" << name << endl;
            cout << sequence_dna.get_dna_str() << endl;
        } else {
            // If the node is internal, iterate through the direct children.
            for (auto child : children) {
                child.sequence_dna.set_codon_seq(sequence_dna.get_codon_seq());
                child.traverse();
            }
        }
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
                string subtree{""};
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
                string name_suffix{""};

                size_t next_sep = newick.substr(position).find(",");
                if (next_sep == string::npos) {
                    name_suffix = newick.substr(position);
                    position = max_position;
                } else {
                    name_suffix = newick.substr(position, next_sep);
                    position = position + next_sep + 1;
                }

                // 'length' contains the name of the node.
                string length{""};
                // 'name' contains the branch length of the node.
                string name{""};
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


int main() {
    string mammals{
            "(Homo_sapiens:1.014243,(Canis_lupus_familiaris:0.844021,(Equus_caballus:0.805325,((Bos_taurus:0.182884"
                    ",((Capra_hircus:0.011884,Capra_aegagrus:0.011884)'0.7-1.9':0.044887,Ovis_aries:0.056771)'4.6-7"
                    ".0':0.126113)'13.6-25.4':0.455399,Sus_scrofa:0.638283)'63.1-64.8':0.167041)'75.6-88.1':0.038697"
                    ")'79.1-92.4':0.170222)'92.5-115.2'"};

    Node root(mammals, 300);
    root.at_equilibrium();
    root.traverse();

    return 0;
}
