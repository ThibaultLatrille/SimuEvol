// Program to print BFS traversal from a given source vertex. BFS(int s) 
// traverses vertices reachable from s.
#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <random>
#include <tuple>
#include "Eigen/Dense"

using namespace Eigen;
typedef Matrix<double, 64, 64> Matrix64x64;
typedef Matrix<double, 64, 1> Vector64;

using namespace std;

// The nucleotides
string const nucleotides{"ATGC"};
// Mapping nucleotides to their positions in the string 'nucleotides'


int triplet_to_codon(int n_1, int n_2, int n_3){
    return n_1*16+n_2*4+n_3;
}

array<tuple<int, int, int>, 64> build_codon_to_triplet_vector() {
    array<tuple<int, int, int>, 64> codon_to_triplet_vec;
    for (unsigned n_1{0}; n_1 < 4; n_1++) {
        for (unsigned n_2{0}; n_2 < 4; n_2++) {
            for (unsigned n_3{0}; n_3 < 4; n_3++) {
                codon_to_triplet_vec[triplet_to_codon(n_1,n_2,n_3)] = make_tuple(n_1,n_2,n_3);
            }
        }
    }
    return codon_to_triplet_vec;
}

array<array<tuple<int, int, int>, 9>, 64> build_codon_to_neighbors_vector(array<tuple<int, int, int>, 64>& codon_to_triplet) {
    array<array<tuple<int, int, int>, 9>, 64> codon_to_neighbors;
    int n_1_from{0}, n_2_from{0}, n_3_from{0};
    for (unsigned codon{0}; codon < 64; codon++) {
        tie(n_1_from, n_2_from, n_3_from) = codon_to_triplet[codon];
        int flag{0};
        for (unsigned n_1_to{0}; n_1_to < 4; n_1_to++) {
            if (n_1_from != n_1_to) {
                codon_to_neighbors[codon][flag] = make_tuple(triplet_to_codon(n_1_to,n_2_from,n_3_from), n_1_from, n_1_to);
                flag++;
            }
        }
        for (unsigned n_2_to{0}; n_2_to < 4; n_2_to++) {
            if (n_2_from != n_2_to) {
                codon_to_neighbors[codon][flag] = make_tuple(triplet_to_codon(n_1_from,n_2_to,n_3_from), n_2_from, n_2_to);
                flag++;
            }
        }
        for (unsigned n_3_to{0}; n_3_to < 4; n_3_to++) {
            if (n_3_from != n_3_to) {
                codon_to_neighbors[codon][flag] = make_tuple(triplet_to_codon(n_1_from,n_2_from,n_3_to), n_3_from, n_3_to);
                flag++;
            }
        }
    }
    return codon_to_neighbors;
}

auto codon_to_triplet_vector(build_codon_to_triplet_vector());
auto const codon_to_neighbors_vector(build_codon_to_neighbors_vector(codon_to_triplet_vector));

string const amino_acids{"ARNDCQEGHILKMFPSTWYVX"};
// Mapping codons to amino-acids
map<string, char> triplet_str_to_aa_char{{"GCA", 'A'},
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

int aa_char_to_aa(const char aa_char, const string& amino_acids){
    return int(amino_acids.find(aa_char));
}

array<int, 64> build_codon_to_aa_vector(map<string, char>& triplet_str_to_aa_char, const string& amino_acids) {
    array<int, 64> codon_to_aa;
    string triplet_str(3, ' ');
    for (unsigned n_1{0}; n_1 < 4; n_1++) {
        triplet_str[0]=nucleotides[n_1];
        for (unsigned n_2{0}; n_2 < 4; n_2++) {
            triplet_str[1]=nucleotides[n_2];
            for (unsigned n_3{0}; n_3 < 4; n_3++) {
                triplet_str[2]=nucleotides[n_3];
                codon_to_aa[triplet_to_codon(n_1,n_2,n_3)] = aa_char_to_aa(triplet_str_to_aa_char[triplet_str], amino_acids);
            }
        }
    }
    return codon_to_aa;
}

array<int, 64> const codon_to_aa_vector(build_codon_to_aa_vector(triplet_str_to_aa_char, amino_acids));

// Random generator engine with seed 0
default_random_engine re(0);

// This class represents a DNA sequence
class Sequence_dna {
private:
    const unsigned length; // The DNA sequence length
    vector<int> codon_seq; // The DNA sequence string (contains only A,T,G and C)
    Matrix<double, 4, 4> mutation_rate_matrix; // The matrix of mutation rates between nucleotides
    vector<double> aa_fitness_profil; // The fitness profil of amino-acids

public:
    // Constructor of Sequence_dna
    Sequence_dna(const unsigned len) : length{len}, codon_seq(len, 0), aa_fitness_profil(21, 0.) {
        // Length is the size of the DNA sequence

        // Uniform random generator between 0 (included) and 3 (included) for assigning nucleotides
        uniform_int_distribution<unsigned char> unif_int(0, 63);

        for (unsigned i{0}; i < length; i++) {
            codon_seq[i] = unif_int(re);
        }

        // Uniform random generator between 0. and 1. for assigning mutation rates
        uniform_real_distribution<double> unif_real(0., 1);

        for (unsigned from{0}; from < 4; from++) { // For each nucleotide (from)
            double total{0};
            for (unsigned to{0}; to < 4; to++) { // For each nucleotide (to)
                // Assign random mutation rates if their is really a mutation (not from G to G for example)
                if (from != to) {
                    mutation_rate_matrix(from, to) = unif_real(re);
                    total -= mutation_rate_matrix(from, to);
                }
            }
            mutation_rate_matrix(from, from) = total;
        }

        for (unsigned aa_index{0}; aa_index < 20; aa_index++) {
            aa_fitness_profil[aa_index] = unif_real(re);
        }

    }

    // Translate the DNA sequence to amino-acid sequence
    string translate() const {
        string protein(length, ' ');

        for (unsigned i{0}; i < length; i++) {
            // Use the codon_to_aa_vector to translate codon to amino-acid
            protein[i] = amino_acids[codon_to_aa_vector[codon_seq[i]]];
        }

        return protein; // return the amino-acid sequence as a string
    }

    double substitution_rate(int codon_from, int codon_to, int n_from, int n_to) {
        double rate_fixation{1.};
        double s{0};

        if (codon_to_aa_vector[codon_to] == 20){
            rate_fixation=0.;
        } else if (codon_from != codon_to){
            s = aa_fitness_profil[codon_to_aa_vector[codon_to]] - aa_fitness_profil[codon_to_aa_vector[codon_from]];
            if (s==.0){
                rate_fixation=1;
            } else {
                rate_fixation=s/(1-exp(-s));
            }
        } else {
            rate_fixation = 1.;
        }

        return rate_fixation * mutation_rate_matrix(n_from, n_to);
    }

    double next_substitution(double &time_left) {

        // The number of transitions possible is 4 times the sequence length (4 nucleotides possible per position)
        unsigned nbr_transitions{9 * length};

        // Initialize the vector of transition rates at 0
        vector<double> transition_rates(nbr_transitions, 0);

        // Initialize the total sum of transition rates at 0
        double total_transition_rates{0.};

        array<tuple<int, int, int>, 9> neighbors;
        int codon_from{0}, codon_to{0}, n_from{0}, n_to{0};

        for (unsigned n{0}; n < length; n++) { // For all indices of the sequence length
            codon_from = codon_seq[n];
            neighbors = codon_to_neighbors_vector[codon_from];
            for (unsigned tr{0}; tr<9; tr++) { // For all possible transitions
                // Don't count if transition to the same nucleotide as already present
                tie(codon_to, n_from, n_to) = neighbors[tr];

                // Assign to the transition rates given by the method substitution rate
                transition_rates[9 * n + tr] = substitution_rate(codon_from, codon_to, n_from, n_to);
                // Increment the total sum of transition rates using the mutation rate 2D array
                total_transition_rates += transition_rates[9 * n + tr];

            }
        }

        // Decrement the time by the total sum of transition rates
        time_left -= 1. / total_transition_rates;

        // Mutate the sequence if the time is positive, else there is no mutation but the time left is set to 0.
        if (time_left > 0. and total_transition_rates != 0.) {
            // Uniform random generator between 0 and the total sum of transition rates
            uniform_real_distribution<double> unif(0, total_transition_rates);
            double random_cumulative_transition_rates = unif(re);
            double cumulative_transition_rates{0.};

            unsigned index = 0;
            for (unsigned t{0}; t < nbr_transitions; t++) {
                // Iterate through the cumulative transition rates and break the loop when it is greater than...
                // ...the random cumulative transition rates
                cumulative_transition_rates += transition_rates[t];
                if (random_cumulative_transition_rates < cumulative_transition_rates) {
                    index = t;
                    break;
                }
            }

            // Mutate the chosen nucleotide with the transition given by the loop break
            neighbors = codon_to_neighbors_vector[codon_seq[index / 9]];
            tie(codon_to, n_from, n_to) = neighbors[index % 9];
            codon_seq[index / 9] = codon_to;

        } else if (time_left < 0.) {
            time_left = 0.;
        }
        return time_left;
    }

    void run_substitutions(double t) {
        while (t > 0) {
            t = next_substitution(t);
        }
    }



    // Generate substitutions until equilibrium by calculating explicitly the equilibrium (fast)
    void burn_in_fast() {
        Matrix64x64 codon_matrix(Matrix64x64::Zero());
        array<tuple<int, int, int>, 9> neighbors;
        int codon_to{0}, n_from{0}, n_to{0};

        for (unsigned codon_from{0}; codon_from < 64; codon_from++) {
            neighbors = codon_to_neighbors_vector[codon_from];
            double total{0};
            for (unsigned tr{0}; tr<9; tr++) { // For all possible transitions
                // Don't count if transition to the same nucleotide as already present
                tie(codon_to, n_from, n_to) = neighbors[tr];
                codon_matrix(codon_to , codon_from) = substitution_rate(codon_from, codon_to, n_from, n_to);
                total -= codon_matrix(codon_to, codon_from);
                // Increment the total sum of transition rates using the mutation rate 2D array
            }
            codon_matrix(codon_from, codon_from) = total;
        }

        FullPivLU<Matrix64x64> lu_decomp(codon_matrix);
        Vector64 kernel(lu_decomp.kernel());

        for (unsigned n{0}; n < length; n++) { // For all indices of the sequence length

            // Initialize the total sum of transition rates at 0
            double total_transition_rates{0.};
            for (unsigned to{0}; to < 64; to++) { // For each nucleotide (to)
                // Increment the total sum of transition rates using the marginal transition rates
                total_transition_rates += kernel(to);
            }

            // Uniform random generator between 0 and the total sum of transition rates
            uniform_real_distribution<double> unif(0, total_transition_rates);
            double random_cumulative_transition_rates = unif(re);
            double cumulative_transition_rates{0.};

            int index{0};
            for (unsigned m{0}; m < 64; m++) {
                // Iterate through the cumulative transition rates and break the loop when it is greater than...
                // ...the random cumulative transition rates
                cumulative_transition_rates += kernel(m);
                if (random_cumulative_transition_rates < cumulative_transition_rates) {
                    index = m;
                    break;
                }
            }

            // Mutate the nucleotide with the transition given by the loop break
            codon_seq[n] = index;
        }
    }

    // Get attribute method for the DNA string
    vector<int> get_codon_seq() const { return codon_seq; }

    // Get attribute method for the DNA string
    string get_dna_str() const {
        string dna_str(length * 3, ' ');
        int n_1{0}, n_2{0}, n_3{0};
        for (unsigned i{0}; i < length; i++) {
            tie(n_1, n_2, n_3) = codon_to_triplet_vector[codon_seq[i]];
            dna_str[3*i] = nucleotides[n_1];
            dna_str[3*i+1] = nucleotides[n_2];
            dna_str[3*i+2] = nucleotides[n_3];
        }
        return dna_str; // return the DNA sequence as a string
    }

    // Set attribute method for the DNA string
    void set_codon_seq(vector<int> const codon) {
        codon_seq = codon;
    }
};


// This class represents nodes of a tree (or a graph)
class Node {
private:
    string name;  // The species name of the node
    double length;  // The length of the branch attached to the node (ascending)
    string newick;  // The newick tree descending from the node
    vector<Node> children;  // Vector of direct children (first order, no grand-children)
    Sequence_dna sequence_dna;  // The DNA sequence attached to the node

public:

    // Constructor

    Node(const string &name, const string &len, const string &newick, const Sequence_dna &seq) :
            name{name}, length{stod(len)}, newick{newick}, sequence_dna{seq} {

        // Parse the newick tree descending of the node
        parse_newick();
    }

    Node(const string &newick, const unsigned length) :
            name{"Root"}, length{0.}, newick{newick}, sequence_dna(length) {

        // Parse the newick tree descending of the node
        parse_newick();
    }

    // Is true if the node don't have children
    bool is_leaf() const {
        return newick.length() == 0;
    }

    // Add a node as the vector of children
    void add_child(const Node &node) {
        children.push_back(node);
    }

    // Recursively iterate through the subtree
    void traverse() {
        // Substitutions of the DNA sequence is generated
        sequence_dna.run_substitutions(length);

        if (is_leaf()) {
            // If the node is a leaf, output the DNA sequence and name
            cout << ">" << name << endl;
            cout << sequence_dna.get_dna_str() << endl;
        } else {
            // If the node is internal, iterate through the direct children
            for (auto child : children) {
                child.sequence_dna.set_codon_seq(sequence_dna.get_codon_seq());
                child.traverse();
            }
        }
    }

    // Set the the DNA sequence to the mutation-selection equilibrium
    void at_equilibrium() {
        sequence_dna.burn_in_fast();
    }

    void parse_newick() {

        if (!is_leaf()) {
            // The size of the string of newick tree
            size_t max_position{newick.size()};
            // The current position in the string of the newick tree
            size_t position{0};

            // While the current position is lower than the size of the string, their is at least one node to parse
            while (position < max_position) {
                // 'subtree' is the left hand side of the node name, it can be a subtree or nothing if the node is a leaf
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

                // 'name_suffix' contains the name of the node and the branch length
                string name_suffix{""};

                size_t next_sep = newick.substr(position).find(",");
                if (next_sep == string::npos) {
                    name_suffix = newick.substr(position);
                    position = max_position;
                } else {
                    name_suffix = newick.substr(position, next_sep);
                    position = position + next_sep + 1;
                }

                // 'length' contains the name of the node
                string length{""};
                // 'name' contains the branch length of the node
                string name{""};
                size_t ddot = name_suffix.rfind(':');
                if (ddot != string::npos) {
                    length = name_suffix.substr(ddot + 1);
                    name = name_suffix.substr(0, ddot);
                } else {
                    name = name_suffix;
                    length = "0";
                }

                // New node from 'subtree', 'name' and 'length' using the DNA sequence of this node
                add_child(Node(name, length, subtree, sequence_dna));
            }
        }
    }
};


// Driver program to test methods of graph class
int main() {
    string mammals{
            "(Homo_sapiens:1.014243,(Canis_lupus_familiaris:0.844021,(Equus_caballus:0.805325,((Bos_taurus:0.182884"
                    ",((Capra_hircus:0.011884,Capra_aegagrus:0.011884)'0.7-1.9':0.044887,Ovis_aries:0.056771)'4.6-7"
                    ".0':0.126113)'13.6-25.4':0.455399,Sus_scrofa:0.638283)'63.1-64.8':0.167041)'75.6-88.1':0.038697"
                    ")'79.1-92.4':0.170222)'92.5-115.2'"};

    Node root(mammals, 999);
    root.at_equilibrium();
    root.traverse();

    return 0;
}
