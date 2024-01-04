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
#include <random>
#include <cmath>
#include <chrono>
#include <functional>

#include "CliqueList.h"

using namespace std;

struct Clique {
    vector<int> nodes;
};

// auto heapCmp = [](pair<int, Clique*>& a, pair<int, Clique*>& b){return a.first>b.first;};

class MyPQ {
public:
    MyPQ() = default;

    void make_heap() {
        std::make_heap(elements.begin(), elements.end(), [](pair<int, Clique*>& a, pair<int, Clique*>& b){return a.first>b.first;});
    }

    void push_back(pair<int, Clique*> e) {
        elements.push_back(e);
    }

    void push(pair<int, Clique*> e) {
        elements.push_back(e);
        std::push_heap(elements.begin(), elements.end(), [](pair<int, Clique*>& a, pair<int, Clique*>& b){return a.first>b.first;});
    }

    void pop() {
        std::pop_heap(elements.begin(), elements.end(), [](pair<int, Clique*>& a, pair<int, Clique*>& b){return a.first>b.first;});
        elements.pop_back();
    }

    pair<int, Clique*>& top() {
        return elements[0];
    }

    bool empty() {
        return elements.empty();
    }

    unsigned long long size() {
        return elements.size();
    }

    vector<pair<int, Clique*>> elements;
};

using PQ = MyPQ;


class ThirdGreedy {
public:
    ThirdGreedy(string fileName, int k);

    void MakeSubGraph(Subgraph& subgraph, int k, int u);
    void BuildGraph(Subgraph& subgraph, int k, int u);
    void FindCliqueMulti(Subgraph& subgraph, int k, vector<int>& prev, Clique** targetClique, int& targetScore, int curScore);

    void FindAllClique(int k, vector<int>& prev, Clique** targetClique, int& targetScore, int curScore);
    void BuildGraph(int k, int u);

    void DeleteNode(int u);
    inline bool check(vector<int>& c, int k);

    void solve(string fileName);

public:
    int n;
    int k;
    int threadNum = 64;

    vector<vector<int>> graph;
    vector<vector<int>> degree;
    vector<vector<int>> nodes;
    vector<int> candidate;

    vector<int> score;
    vector<int> valid;
    vector<int> orderedVertex;
    vector<int> vertexRank;

    vector<vector<int>> output;
};