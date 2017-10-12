// Program to print BFS traversal from a given source vertex. BFS(int s) 
// traverses vertices reachable from s.
#include <iostream>
#include <vector>
#include <deque>
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
        time_left -= 1./total_transition_rates;

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

    void substitute_for(double t){
        while (t>0){
            t = next_substitution(t);
        }
    }
    // Generate substitutions until equilibrium using the method 'next_substitution' (slow)
    void burn_in() {
        double d = INFINITY;
        for (unsigned m{0}; m < 4 * length; m++) {
            next_substitution(d);
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
};

class Node {
private:

public:
    string name;
    double length;
    string subtree;
    Sequence_dna sequence_dna;
    unsigned long index;

    Node(const string & name, const string & length, const string & subtree, const Sequence_dna & seq) :
            name(name), length(stod(length)), subtree(subtree), sequence_dna(seq) {
        static unsigned long static_index{0};
        index = static_index;
        static_index++;
    }  // Constructor

    bool is_leaf(){
        return subtree.length() == 0;
    }
};


// This class represents a directed graph using adjacency list representation
class Graph {

private:
    vector<vector<unsigned long>> adj;    // Pointer to an array containing adjacency lists
    vector<Node> nodes;
public:

    Graph(const string newick, unsigned int length){
        Node root("ROOT", "0", newick, Sequence_dna(length));
        root.sequence_dna.burn_in_fast();
        nodes.push_back(root);
        adj.push_back(vector<unsigned long>());

        parse_newick(newick, root);
    }  // Constructor

    Graph(const Graph &) = delete; // copy constructor

    void addNode(Node node, Node & parent) {
        node.sequence_dna.substitute_for(node.length);
        nodes.push_back(node);
        adj.push_back(vector<unsigned long>());

        addEdge(parent, node);
        if (!node.is_leaf()){
            parse_newick(node.subtree, node);
        }
    }

    void addEdge(const Node & from, Node & to) {
        adj[from.index].push_back(to.index); // Add w to v’s list.
    } // function to add an edge to graph

    void parse_newick(string newick, Node & parent) {
        unsigned long max_point{newick.size()};
        unsigned long point{0};

        while (point < max_point) {
            string subtree{""};
            if (newick[point] == '(') {

                unsigned long postpoint{point};
                unsigned long nbr_open{1};

                for (unsigned long i{point+1}; i < max_point; i++) {
                    if (nbr_open == 0) {
                        postpoint = i;
                        break;
                    } else if (newick[i] == '(') {
                        nbr_open++;
                    } else if (newick[i] == ')') {
                        nbr_open--;
                    };
                }
                subtree = newick.substr(point+1, postpoint - point - 2);
                point = postpoint;
            }

            long next_sep = newick.substr(point).find(",");
            string name_suffix{""};

            if (next_sep == -1){
                name_suffix = newick.substr(point);
                point = max_point;
            } else {
                name_suffix = newick.substr(point, (unsigned long)next_sep);
                point = point + (unsigned long)next_sep + 1;
            }

            long ddot = name_suffix.rfind(':');
            string length{""};
            string name{""};

            if (ddot != -1) {
                length = name_suffix.substr((unsigned long)ddot + 1);
                name = name_suffix.substr(0, (unsigned long)ddot);
            } else {
                name = name_suffix;
                length = "0";
            }

            addNode(Node(name, length, subtree, parent.sequence_dna), parent);
        }
    }

    void BFS(unsigned long s);  // prints BFS traversal from a given source s
};

void Graph::BFS(unsigned long s) {
    // Mark all the vertices as not visited
    vector<bool> visited(nodes.size(), false);

    // Create a queue for BFS
    deque<unsigned long> queue;

    // Mark the current node as visited and enqueue it
    visited[s] = true;
    queue.push_back(s);

    while (!queue.empty()) {
        // Dequeue a vertex from queue and print it
        unsigned long f = queue.front();

        if (nodes[f].is_leaf()) {
            cout << ">" << nodes[f].name << endl;
            cout << nodes[f].sequence_dna.get_dna() << endl;
            //cout << "AA: " << nodes[f].sequence_dna.translate() << endl;
        }

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

// Driver program to test methods of graph class
int main() {
    string mammals{
            "(Homo_sapiens:1.014243,(Canis_lupus_familiaris:0.844021,(Equus_caballus:0.805325,((Bos_taurus:0.182884,((Capra_hircus:0.011884,Capra_aegagrus:0.011884)'0.7-1.9':0.044887,Ovis_aries:0.056771)'4.6-7.0':0.126113)'13.6-25.4':0.455399,Sus_scrofa:0.638283)'63.1-64.8':0.167041)'75.6-88.1':0.038697)'79.1-92.4':0.170222)'92.5-115.2'"};
    string newick{"(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F"};
    string wnewick{"((F:0.3,G:0.4)A:0.5,B:0.2,(C:0.3,D:0.4)E:0.5)F"};


    Graph g(mammals, 999);
    g.BFS(0);

    return 0;
}
