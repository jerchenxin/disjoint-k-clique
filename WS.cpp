// g++ -O3 -std=c++17 WS.cpp -o WS
// Watts-Strogatz model

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <set>

using namespace std;

// Edge structure
struct Edge {
    int source;
    int target;
};

// Check if an edge already exists in the graph
bool isEdgeExist(const vector<Edge>& graph, int source, int target) {
    for (const Edge& edge : graph) {
        if (edge.source == source && edge.target == target)
            return true;
    }
    return false;
}

// Generate Watts-Strogatz small-world graph
vector<Edge> generateWattsStrogatzGraph(int numNodes, int numNeighbors, double rewiringProbability) {
    vector<Edge> graph;

    // Create initial ring lattice
    for (int i = 0; i < numNodes; ++i) {
        for (int j = 1; j <= numNeighbors / 2; ++j) {
            int target = (i + j) % numNodes;
            Edge newEdge = { i, target };
            graph.push_back(newEdge);
        }
    }

    // Rewire edges with a certain probability
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);

    for (int i = 0; i < numNodes; ++i) {
        for (int j = 0; j < numNeighbors / 2; ++j) {
            if (dis(gen) < rewiringProbability) {
                int target = (i + j) % numNodes;
                int newTarget;

                do {
                    newTarget = static_cast<int>(dis(gen) * numNodes);
                } while (newTarget == i || isEdgeExist(graph, i, newTarget));

                graph[i * numNeighbors / 2 + j].target = newTarget;
            }
        }
    }

    return graph;
}


// Write the graph to a file
void writeGraphToFile(const vector<Edge>& graph, const string& filename) {
    ofstream file(filename);
    if (!file) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    // Write the edge information
    for (const Edge& edge : graph) {
        file << edge.source << " " << edge.target << endl;
    }

    file.close();
    cout << "Graph written to file: " << filename << endl;
}


int main(int argc, char* argv[]) {
    int numNodes = atoi(argv[1]);               // Number of nodes
    int numNeighbors = atoi(argv[2]);            // Number of initial neighbors (even number)
    double rewiringProbability(atof(argv[3])); // Probability of rewiring an edge
    string filename(argv[4]);

    vector<Edge> graph = generateWattsStrogatzGraph(numNodes, numNeighbors, rewiringProbability);

    // Write the graph to a file
    writeGraphToFile(graph, filename);

    return 0;
}