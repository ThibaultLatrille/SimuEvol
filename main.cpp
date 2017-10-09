// Program to print BFS traversal from a given source vertex. BFS(int s) 
// traverses vertices reachable from s.
#include <iostream>
#include <vector>
#include <deque>
#include <map>

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

string nucleotides{"ATGC"};
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

// This class represents a DNA sequence
class Sequence_dna {
private:
    string dna;    // the DNA sequence string (contains only A,T,G and C)

public:

    Sequence_dna(const unsigned length) : dna(length, ' ') {
        srand(0); //seeds the rand() function
        for (int i = 0; i < length; i++) { dna[i] = nucleotides[rand() % 4]; }
    }  // Constructor

    Sequence_dna(const Sequence_dna &) = delete; // copy constructor

    bool is_coding() { return dna.size() % 3 == 0; }

    string translate() {
        if (!is_coding()) { throw logic_error("The DNA sequence is not coding (not a multiple of 3)"); }
        string protein(dna.size() / 3, ' ');
        for (unsigned i = 0; i < protein.size(); i++) {
            protein[i] = codon_table[dna.substr(3 * i, 3)];
        }
        return protein;
    }

    string get_dna() { return dna; }
};


// Driver program to test methods of graph class
int main() {

    Sequence_dna seq(99);
    string dna = seq.get_dna();
    cout << dna << endl;

    string protein = seq.translate();
    cout << protein << endl;

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
    vector<Graph> v{};

    return 0;
}
