// Program to print BFS traversal from a given source vertex. BFS(int s) 
// traverses vertices reachable from s.
#include <iostream>
#include <list>
#include <vector>

using namespace std;

// This class represents a directed graph using adjacency list representation
class Graph {
    unsigned V;    // No. of vertices
    vector<vector<int>> adj;    // Pointer to an array containing adjacency lists
public:
    Graph(const unsigned vertices): V(vertices), adj(vertices, vector<int>()) {}  // Constructor
    Graph(const Graph&) = delete; // copy constructor
    void addEdge(int v, int w) {
        adj[v].push_back(w); // Add w to v’s list.
    } // function to add an edge to graph
    void BFS(int s);  // prints BFS traversal from a given source s
};

void Graph::BFS(int s) {
    // Mark all the vertices as not visited
    vector<bool> visited(V, false);

    // Create a queue for BFS
    list<int> queue;

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

// Driver program to test methods of graph class
int main() {
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
