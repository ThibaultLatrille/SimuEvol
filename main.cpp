// Program to print BFS traversal from a given source vertex. BFS(int s) 
// traverses vertices reachable from s.
#include <iostream>
#include <vector>
#include <map>
#include <random>

using namespace std;

// The nucleotides
string nucleotides{"ATGC"};
// Mapping nucleotides to their positions in the string 'nucleotides'
map<char, unsigned> nucleotide_table{{'A', 0},
                                     {'T', 1},
                                     {'G', 2},
                                     {'C', 3}};

// Mapping codons to amino-acids
map<string, char> codon_table{{"GCA", 'A'},
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

// Random generator engine with seed 0
default_random_engine re(0);

// This class represents a DNA sequence
class Sequence_dna {
private:
    string dna_str; // The DNA sequence string (contains only A,T,G and C)
    const unsigned length; // The DNA sequence length
    double mutation_rate[4][4]; // The matrix of mutation rates between nucleotides

public:
    // Constructor of Sequence_dna
    Sequence_dna(const unsigned len) : length(len), dna_str(len, ' ') {
        // Length is the size of the DNA sequence

        // Uniform random generator between 0 (included) and 3 (included) for assigning nucleotides
        uniform_int_distribution<int> unif_int(0, 3);

        for (unsigned i{0}; i < length; i++) {
            dna_str[i] = nucleotides[unif_int(re)]; // Assign random nucleotides (A,T,G or C)
        }

        // Uniform random generator between 0. and 1. for assigning mutation rates
        uniform_real_distribution<double> unif_real(0., 1);

        for (unsigned from{0}; from < 4; from++) { // For each nucleotide (from)
            for (unsigned to{0}; to < 4; to++) { // For each nucleotide (to)
                // Assign random mutation rates if their is really a mutation (not from G to G for example)
                if (from != to) {
                    mutation_rate[from][to] = unif_real(re);
                } else {
                    mutation_rate[from][to] = 0.;
                }

            }
        }
    }

    // Is Sequence_dna a coding sequence ?
    bool is_coding() { return length % 3 == 0; }

    // Translate the DNA sequence to amino-acid sequence
    string translate() {
        // Assert that the sequence is a multiple of 3.
        if (!is_coding()) { throw logic_error("The DNA sequence is not coding (not a multiple of 3)"); }

        // The amino-acid sequence is 3 times shorter than the DNA sequence
        unsigned protein_size{length / 3};
        string protein(length / 3, ' ');

        for (unsigned i{0}; i < protein_size; i++) {
            // Use the codon_table to translate codon to amino-acid
            protein[i] = codon_table[dna_str.substr(3 * i, 3)];
        }

        return protein; // return the amino-acid sequence as a string
    }

    double next_substitution(double &time_left) {

        // The number of transitions possible is 4 times the sequence length (4 nucleotides possible per position)
        unsigned nbr_transitions{4 * length};

        // Initialize the vector of cumulative transition rates at 0
        vector<double> cumulative_transition_rates(nbr_transitions, 0.);

        // Initialize the total sum of transition rates at 0
        double total_transition_rates{0.};

        for (unsigned n{0}; n < length; n++) { // For all indices of the sequence length
            unsigned current_n = nucleotide_table[dna_str[n]]; // The nucleotide at position n

            for (unsigned t{0}; t < 4; t++) { // For all possible transitions
                // Don't count if transition to the same nucleotide as already present
                if (t != current_n) {
                    // Increment the total sum of transition rates using the mutation rate 2D array
                    total_transition_rates += mutation_rate[current_n][t];
                    // Assign to the cumulative transition rates the current total sum
                    cumulative_transition_rates[4 * n + t] = total_transition_rates;
                }
            }
        }

        // Decrement the time by the total sum of transition rates
        time_left -= 1. / total_transition_rates;

        // Mutate the sequence if the time is positive, else there is no mutation but the time left is set to 0.
        if (time_left > 0. and total_transition_rates != 0.) {
            // Uniform random generator between 0 and the total sum of transition rates
            uniform_real_distribution<double> unif(0, total_transition_rates);
            double random_cumulative_transition_rates = unif(re);

            unsigned index = 0;
            for (unsigned t{0}; t < nbr_transitions; t++) {
                // Iterate through the cumulative transition rates and break the loop when it is greater than...
                // ...the random cumulative transition rates
                if (random_cumulative_transition_rates < cumulative_transition_rates[t]) {
                    index = t;
                    break;
                }
            }

            // Mutate the chosen nucleotide with the transition given by the loop break
            dna_str[index / 4] = nucleotides[index % 4];

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
        for (unsigned n{0}; n < length; n++) { // For all indices of the sequence length

            // Initialize the total sum of transition rates at 0
            double total_transition_rates{0.};
            // Initialize the vector of cumulative transition rates at 0
            vector<double> cumulative_transition_rates(4, 0.);

            for (unsigned to{0}; to < 4; to++) { // For each nucleotide (to)
                // Initialize the marginal sum of transition rates at 0
                double marginal_transition_rates{0.};
                for (unsigned from{0}; from < 4; from++) { // For each nucleotide (from)
                    // Increment the marginal transition rates using the mutation rate 2D array
                    if (from != to) {
                        marginal_transition_rates += mutation_rate[from][to];
                    }
                }
                // Increment the total sum of transition rates using the marginal transition rates
                total_transition_rates += marginal_transition_rates;
                // Assign to the cumulative transition rates the current total sum
                cumulative_transition_rates[to] = total_transition_rates;
            }

            // Uniform random generator between 0 and the total sum of transition rates
            uniform_real_distribution<double> unif(0, total_transition_rates);
            double random_cumulative_transition_rates = unif(re);

            unsigned index = 0;
            for (unsigned m{0}; m < 4; m++) {
                // Iterate through the cumulative transition rates and break the loop when it is greater than...
                // ...the random cumulative transition rates
                if (random_cumulative_transition_rates < cumulative_transition_rates[m]) {
                    index = m;
                    break;
                }
            }

            // Mutate the nucleotide with the transition given by the loop break
            dna_str[n] = nucleotides[index];
        }
    }

    // Get attribute method for the DNA string
    string get_dna() { return dna_str; }

    // Set attribute method for the DNA string
    void set_dna(string dna) {dna_str=dna;}
};


// This class represents nodes of a tree (or a graph)
class Node {
private:
    string name;  // The species name of the node
    double length;  // The length of the branch attached to the node (ascending)
    string newick;  // The newick tree descending from the node
    Sequence_dna sequence_dna;  // The DNA sequence attached to the node
    vector<Node> children;  // Vector of direct children (first order, no grand-children)

public:
    // Constructor
    Node(const string &name, const string &len, const string &newick, const Sequence_dna &seq) :
            name(name), length(stod(len)), newick(newick), sequence_dna(seq) {

        // If the node is internal, parse the newick tree descending of the node
        if (!is_leaf()) {
            parse_newick();
        }
    }

    // Is true if the node don't have children
    bool is_leaf() {
        return newick.length() == 0;
    }

    // Add a node as the vector of children
    void add_child(Node node) {
        children.push_back(node);
    }

    // Recursively iterate through the subtree
    void traverse() {
        // Substitutions of the DNA sequence is generated
        sequence_dna.run_substitutions(length);

        if (is_leaf()) {
            // If the node is a leaf, output the DNA sequence and name
            cout << ">" << name << endl;
            cout << sequence_dna.get_dna() << endl;
        } else {
            // If the node is internal, iterate through the direct children
            for (auto child : children) {
                child.sequence_dna.set_dna(sequence_dna.get_dna());
                child.traverse();
            }
        }
    }

    // Set the the DNA sequence to the mutation-selection equilibrium
    void at_equilibrium() {
        sequence_dna.burn_in_fast();
    }

    void parse_newick() {
        // The size of the string of newick tree
        unsigned long max_position{newick.size()};
        // The current position in the string of the newick tree
        unsigned long position{0};

        // While the current position is lower than the size of the string, their is at least one node to parse
        while (position < max_position) {
            // 'subtree' is the left hand side of the node name, it can be a subtree or nothing if the node is a leaf
            string subtree{""};
            if (newick[position] == '(') {

                unsigned long postpoint{position};
                unsigned long nbr_open{1};

                for (unsigned long i{position + 1}; i < max_position; i++) {
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

            long next_sep = newick.substr(position).find(",");
            if (next_sep == -1) {
                name_suffix = newick.substr(position);
                position = max_position;
            } else {
                name_suffix = newick.substr(position, (unsigned long) next_sep);
                position = position + (unsigned long) next_sep + 1;
            }

            // 'length' contains the name of the node
            string length{""};
            // 'name' contains the branch length of the node
            string name{""};
            long ddot = name_suffix.rfind(':');
            if (ddot != -1) {
                length = name_suffix.substr((unsigned long) ddot + 1);
                name = name_suffix.substr(0, (unsigned long) ddot);
            } else {
                name = name_suffix;
                length = "0";
            }

            // New node from 'subtree', 'name' and 'length' using the DNA sequence of this node
            add_child(Node(name, length, subtree, sequence_dna));
        }
    }
};


// Driver program to test methods of graph class
int main() {
    string mammals{
            "(Homo_sapiens:1.014243,(Canis_lupus_familiaris:0.844021,(Equus_caballus:0.805325,((Bos_taurus:0.182884,((Capra_hircus:0.011884,Capra_aegagrus:0.011884)'0.7-1.9':0.044887,Ovis_aries:0.056771)'4.6-7.0':0.126113)'13.6-25.4':0.455399,Sus_scrofa:0.638283)'63.1-64.8':0.167041)'75.6-88.1':0.038697)'79.1-92.4':0.170222)'92.5-115.2'"};
    string newick{"(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F"};
    string wnewick{"((F:0.3,G:0.4)A:0.5,B:0.2,(C:0.3,D:0.4)E:0.5)F"};


    Node root("ROOT", "0", mammals, Sequence_dna(999));
    root.at_equilibrium();
    root.traverse();

    return 0;
}
