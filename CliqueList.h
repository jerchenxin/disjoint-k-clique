#pragma once

#include <vector>
#include <unordered_set>
#include <string>

using namespace std;

struct Subgraph {
    vector<vector<int>> graph;
    vector<vector<int>> degree;
    vector<vector<int>> nodes;
    vector<int> candicate;
    vector<int> score;

    Subgraph(vector<vector<int>>& degree, vector<vector<int>>& nodes, int n) : degree(degree), nodes(nodes) {
        graph.resize(n);
        candicate.resize(n, 0);
        score.resize(n, 0);
    }

    Subgraph(vector<vector<int>>& graph, vector<vector<int>>& degree, vector<vector<int>>& nodes, int n) : graph(graph), degree(degree), nodes(nodes) {
        candicate.resize(n, 0);
        score.resize(n, 0);
    }
};

class CliqueList {
public:
    CliqueList(string fileName, int k);
    void ListSingleThread();
    void ListMultiThread(int threadNum);

    void CountSingle();
    void CountMulti(int threadNum);

    void CountScoreSingle();
    void CountScoreMulti(int threadNum);

    // Prune
    inline int Prune(vector<int>& degree, int u);
    inline int Prune(Subgraph& subgraph, vector<int>& degree, int u);

    // color
    void InitColor();
    int UniqueColor(vector<int>& degree, int u);
    int UniqueColor(vector<vector<int>>& graph, vector<int>& degree, int u);

    // ListSingleThread
    void FindClique(int k, vector<int>& prev);
    void BuildGraph(int k, int u);

    // ListMultiThread
    void FindClique(Subgraph& subgraph, int k, vector<int>& prev, vector<vector<int>>& output);
    void BuildGraph(Subgraph& subgraph, int k, int u);
    void MakeSubGraph(Subgraph& subgraph, int k, int u);

    // CountSingle
    void FindCliqueCount(int k);

    // CountMulti
    void FindClique(Subgraph& subgraph, int k, unsigned long long& num);

    // CountScoreSingle
    void FindCliqueCountScore(int k, vector<int>& prev);

    // CountScoreMulti
    void FindCliqueCountScore(Subgraph& subgraph, int k, vector<int>& prev);
    void FindCliqueCountScoreMulti(Subgraph& subgraph, int k, vector<int>& prev, vector<int>& score);

public:
    vector<vector<int>> graph;
    vector<vector<int>> degree;
    vector<vector<int>> nodes;
    vector<int> color;
    int n;
    int k;

    vector<vector<int>> output;
    unsigned long long cliqueNum;
    vector<int> score;
};