// g++ -O3 -std=c++17 greedyOriginalDegreeDAG.cpp -o greedyOriginalDegreeDAG 

#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <random>
#include <cmath>
#include <chrono>

using namespace std;

int n;

vector<vector<int>> output;
vector<vector<int>> graph;
vector<vector<int>> degree;
vector<vector<int>> nodes;
vector<int> valid;

bool FindClique(int k, vector<int>& prev);
void BuildGraph(int k, int u);
void DeleteNode(int u);

int main(int argc, char* argv[]) {
    auto start = std::chrono::high_resolution_clock::now();

    int k = atoi(argv[2]);
    unsigned long long m = 0;
    n = 0;
    vector<pair<int, int>> tmpG;
    
    ifstream f(argv[1]);
    int u, v;
    while (f >> u >> v) {
        m++;
        tmpG.emplace_back(u, v);
        n = max(max(n, u), v);
    }
    
    n++;

    degree.resize(k+1, vector<int>(n, 0));
    nodes.resize(k+1);

    for (auto [u, v] : tmpG) {
        degree[k][u]++;
        degree[k][v]++;
    }

#ifndef PRINT_CLIQUE
    cout << "n: " << n << endl;
#endif

    valid.resize(n, 1);
    graph.resize(n);

    // rank
    vector<pair<int, int>> rank;
    for (int i=0;i<n;i++) {
        rank.emplace_back(degree[k][i], i);
    }
    sort(rank.begin(), rank.end());

    vector<int> r(n);
    for (int i=0;i<n;i++) {
        r[rank[i].second] = i;
    } 

    for (auto [u, v] : tmpG) {
        if (r[u] > r[v]) {
            graph[u].emplace_back(v);
        } else {
            graph[v].emplace_back(u);
        }
    }
    
    // modified degree
    for (int i=0;i<n;i++) {
        degree[k][i] = graph[i].size();
    }

    for (auto [i, id] : rank) {
        if (valid[id]) {
            if (degree[k][id] < k-1) {
                continue;
            }
            BuildGraph(k, id);
            vector<int> prev = {id};
            if (FindClique(k-1, prev)) {
                DeleteNode(id);
                output.emplace_back(move(prev));
            }
        }
    }
    
#ifndef PRINT_CLIQUE
    auto end = std::chrono::high_resolution_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    cout << "Time(ms): " << diff.count() / 1000000 << endl;
    
    cout << "Solution Number: " << output.size() << endl;
#else
    for (auto& c : output) {
        for (auto i : c) {
            cout << i << " ";
        }
        cout << endl;
    }
#endif
}

bool FindClique(int k, vector<int>& prev) {
    if (k == 2) {
        for (auto i : nodes[k]) {
            for (int j=0;j<degree[k][i];j++) {
                prev.emplace_back(i);
                prev.emplace_back(graph[i][j]);
                DeleteNode(i);
                DeleteNode(graph[i][j]);
                return true;
            }
        }
    } else {
        for (auto i : nodes[k]) {
            if (degree[k][i] < k-1) {
                continue;
            }

            BuildGraph(k, i);

            prev.push_back(i);
            if (FindClique(k-1, prev)) {
                DeleteNode(i);
                return true;
            }
            prev.pop_back();
        }
    }
    return false;
}

void BuildGraph(int k, int u) {
    static vector<int> candidate;
    if (candidate.size() < n) {
        candidate.resize(n, 0);
    }

    nodes[k-1].clear();

    for (int i=0;i<degree[k][u];i++) {
        int v = graph[u][i];
        if (valid[v]) {
            candidate[v] = 1;
            nodes[k-1].push_back(v);
        }
    }

    for (int i=0;i<degree[k][u];i++) {
        int v = graph[u][i];
        int nextIndex = 0;
        for (int j=0;j<degree[k][v];j++) {
            if (candidate[graph[v][j]] == 1) {
                swap(graph[v][j], graph[v][nextIndex]);
                nextIndex++;
            }
        }
        degree[k-1][v] = nextIndex;
    }

    for (int i=0;i<degree[k][u];i++) {
        int v = graph[u][i];
        candidate[v] = 0;
    }
}

void DeleteNode(int u) {
    valid[u] = 0;
}