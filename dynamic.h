#pragma once

#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <array>
#include <random>
#include <cmath>
#include <chrono>
#include <functional>

#include "lightPruneClass.h"


using namespace std;


struct KSubgraph {
    vector<vector<int>> graph;
    vector<vector<int>> degree;
    vector<set<int>> nodes;
    vector<int> candidate;

    KSubgraph() {}

    KSubgraph(int k, int n) {
        graph.resize(n);
        degree.resize(k+1, vector<int>(n, 0));
        nodes.resize(k+1);
        candidate.resize(n, 0);
    }
};


class DynamicGraph {
public:
    DynamicGraph(string fileName, int k);

    void BuildIndex(string fileName);
    void FindIndexByClique(vector<int>& c);

    void BuildDynamicIndex();
    void MakeSubGraph(KSubgraph& subgraph, int k, vector<int>& c);
    void BuildSubGraph(KSubgraph& subgraph, int k, int u);
    void FindAllClique(KSubgraph& subgraph, int k, vector<int>& prev);
    void FindAllCliqueMulti(KSubgraph& subgraph, int k, vector<int>& prev, vector<set<Clique*>>& tmpIndex);

    void Insertion(int u, int v);
    void AddEdge(int u, int v);
    void MakeSubGraph(KSubgraph& subgraph, int k, int s, int t);
    void FindAllClique(KSubgraph& subgraph, int k, vector<int>& prev, set<int>& modifiedGroup, int s1, int t1);
    void Swap(set<int>& modifiedGroup);
    void TrySwap(int groupID, set<int>& modifiedGroup);
    void MakeSubGraph(KSubgraph& subgraph, int k, int groupID);
    void MakeFreeSubgraph(KSubgraph& subgraph, int k, int s, int t);
    void FindAllFreeClique(KSubgraph& subgraph, int k, vector<int>& prev, set<Clique*>& freeCliqueList);
    void MakeSubgraphByFreeEdge(KSubgraph& subgraph, int k, int s, int t);
    void FindAllCliqueByFreeEdge(KSubgraph& subgraph, int k, vector<int>& prev, set<int>& modifiedGroup, int s1, int t1);

    void Deletion(int u, int v);
    void DeleteEdge(int u, int v);
    void MakeSubgraphByFreeNodes(KSubgraph& subgraph, int k, set<int>& freeNodes, set<int>& newGroupIDSet);
    void FindAllCliqueByFreeNodes(KSubgraph& subgraph, int k, vector<int>& prev, set<int>& modifiedGroup, set<int>& freeNodes);
    void ExtendIndexByFreeNodes(set<int>& newFreeNodes, set<int>& completeGroup, set<int>& modifiedGroup);

    void DeleteIndex(Clique* c);

    vector<vector<int>> DisjointClique(set<Clique*>& cs);

    void Info();

    void IndexCorrectness();

public:
    int n;
    int k;

    KSubgraph subgraph;

    vector<vector<int>> graph;
    vector<vector<int>> reverseGraph;
    vector<int> vertexRank;

    vector<set<Clique*>> index;
    vector<int> group;
    map<int, vector<int>> invGroup;

    vector<vector<int>> output;
};

