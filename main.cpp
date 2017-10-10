// Program to print BFS traversal from a given source vertex. BFS(int s) 
// traverses vertices reachable from s.
#include <iostream>
#include <vector>
#include <deque>
#include <map>
#include <random>

using namespace std;

// This class represents a directed graph using adjacency list representation
class Graph {

private:
    unsigned V;    // No. of vertices
    vector<vector<int>> adj;    // Pointer to an array containing adjacency lists

public:

    Graph(const unsigned vertices) : V(vertices), adj(vertices, vector<int>()) {}  // Constructor

    Graph(const Graph &) = delete; // copy constructor

    void addEdge(int v, int w) {
        adj[v].push_back(w); // Add w to v’s list.
    } // function to add an edge to graph

    void BFS(int s);  // prints BFS traversal from a given source s
};

void Graph::BFS(int s) {
    // Mark all the vertices as not visited
    vector<bool> visited(V, false);

    // Create a queue for BFS
    deque<int> queue;

    // Mark the current node as visited and enqueue it
    visited[s] = true;
    queue.push_back(s);

    while (!queue.empty()) {
        // Dequeue a vertex from queue and print it
        int f = queue.front();
        cout << f << endl;
        queue.pop_front();

        // Get all adjacent vertices of the dequeued vertex s
        // If a adjacent has not been visited, then mark it visited
        // and enqueue it
        for (auto i : adj[f]) {
            if (!visited[i]) {
                visited[i] = true;
                queue.push_back(i);
            }
        }
    }
}

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
        uniform_real_distribution<double> unif_real(0., 1.);

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

    // Copy constructor of Sequence_dna
    Sequence_dna(const Sequence_dna &) = delete;

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
        time_left -= total_transition_rates;

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

        } else if (time_left < 0.){
            time_left = 0.;
        }
        return time_left;
    }

    // Generate substitutions until equilibrium using the method 'next_substitution' (slow)
    void burn_in(){
        double d = INFINITY;
        for (unsigned m{0}; m < 4 * length; m++){
            next_substitution(d);
        }
    }

    // Generate substitutions until equilibrium by calculating explicitly the equilibrium (fast)
    void burn_in_fast(){
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
};


// Driver program to test methods of graph class
int main() {

    Sequence_dna seq(999);
    cout << seq.get_dna() << endl;
    seq.burn_in_fast();
    string dna = seq.get_dna();
    string protein = seq.translate();
    cout << protein << endl;

    double t = 400;
    t = seq.next_substitution(t);

    cout << dna << endl;
    cout << seq.get_dna() << endl;
    cout << t << endl;

    // Create a graph given in the above diagram
    Graph g(4);
    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(1, 2);
    g.addEdge(2, 0);
    g.addEdge(2, 3);
    g.addEdge(3, 3);

    cout << "Following is Breadth First Traversal "
         << "(starting from vertex 2)\n";
    g.BFS(2);

    return 0;
}
