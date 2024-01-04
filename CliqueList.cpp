#include "CliqueList.h"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <chrono>
#include <omp.h>
#include <queue>


CliqueList::CliqueList(string fileName, int k) : k(k), cliqueNum(0) {
    ifstream f(fileName);
    int u, v;
    n = 0;
    while (f >> u >> v) {
        n = max(max(n, u), v);

        if (graph.size() < n + 1) {
            graph.resize(n + 1);
        }

        if (u > v) // big : small small small
            graph[u].emplace_back(v);
        else 
            graph[v].emplace_back(u);
    }
    n++;

    degree.resize(k+1, vector<int>(n, 0));
    nodes.resize(k+1);
    score.resize(n, 0);

    for (int i=0;i<n;i++) {
        degree[k][i] = graph[i].size();
    }

#ifdef USE_COLOR
    InitColor();
#endif
}

void CliqueList::ListSingleThread() {
    for (int i=0;i<n;i++) {
        if (Prune(degree[k], i) < k-1) {
            continue;
        }

        BuildGraph(k, i);
        vector<int> prev = {i};
        FindClique(k-1, prev);
    }
}

void CliqueList::ListMultiThread(int threadNum) {
    omp_set_num_threads(threadNum);

    int i;

    #pragma omp parallel private(i)
    {
        int pid = omp_get_thread_num();
        vector<vector<int>> res;
        Subgraph subgraph(degree, nodes, n);
        #pragma omp for schedule(dynamic, 1) nowait
        for (i=0;i<n;i++) {
            if (Prune(subgraph, degree[k], i) < k-1) {
                continue;
            }

            MakeSubGraph(subgraph, k, i);
            vector<int> prev = {i};
            FindClique(subgraph, k-1, prev, res);
        }

        #pragma omp critical
        output.insert(output.end(), res.begin(), res.end());
    }
}

void CliqueList::CountSingle() {
    for (int i=0;i<n;i++) {
        if (Prune(degree[k], i) < k-1) {
            continue;
        }

        BuildGraph(k, i);
        FindCliqueCount(k-1);
    }
}

void CliqueList::CountMulti(int threadNum) {
    omp_set_num_threads(threadNum);

    unsigned long long num = 0;
    int i;

    #pragma omp parallel private(i) reduction(+:num)
    {
        int pid = omp_get_thread_num();
        Subgraph subgraph(degree, nodes, n);
        #pragma omp for schedule(dynamic, 1) nowait
        for (i=0;i<n;i++) {
            if (Prune(subgraph, degree[k], i) < k-1) {
                continue;
            }

            MakeSubGraph(subgraph, k, i);
            FindClique(subgraph, k-1, num);
        }
    }

    cliqueNum = num;
}

void CliqueList::CountScoreSingle() {
    for (int i=0;i<n;i++) {
        if (Prune(degree[k], i) < k-1) {
            continue;
        }

        BuildGraph(k, i);
        vector<int> prev = {i};
        FindCliqueCountScore(k-1, prev);
    }
}

void CliqueList::CountScoreMulti(int threadNum) {
    omp_set_num_threads(threadNum);

    int i;

    #pragma omp parallel private(i)
    {
        int pid = omp_get_thread_num();
        Subgraph subgraph(degree, nodes, n);
        #pragma omp for schedule(dynamic, 1) nowait
        for (i=0;i<n;i++) {
            if (Prune(subgraph, degree[k], i) < k-1) {
                continue;
            }

            MakeSubGraph(subgraph, k, i);
            vector<int> prev = {i};
            FindCliqueCountScore(subgraph, k-1, prev);
            // FindCliqueCountScoreMulti(subgraph, k-1, prev, score);
        }

        #pragma omp critical
        {
            for (int j=0;j<n;j++) {
                score[j] += subgraph.score[j];
            }
        }
    }
}

void CliqueList::FindClique(int k, vector<int>& prev) {
    if (k == 2) {
        for (auto i : nodes[k]) {
            for (int j=0;j<degree[k][i];j++) {
                vector<int> tmp = prev;
                tmp.push_back(i);
                tmp.push_back(graph[i][j]);
                output.emplace_back(move(tmp));
            }
        }
    } else {
        for (auto i : nodes[k]) {
            if (Prune(degree[k], i) < k-1) {
                continue;
            }

            BuildGraph(k, i);
            prev.push_back(i);
            FindClique(k-1, prev);
            prev.pop_back();
        }
    }
}

void CliqueList::BuildGraph(int k, int u) {
    static vector<int> candidate;
    if (candidate.size() < n) {
        candidate.resize(n, 0);
    }

    nodes[k-1].clear();

    for (int i=0;i<degree[k][u];i++) {
        int v = graph[u][i];
        candidate[v] = 1;
        nodes[k-1].push_back(v);
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

void CliqueList::FindClique(Subgraph& subgraph, int k, vector<int>& prev, vector<vector<int>>& output) {
    auto& nodes = subgraph.nodes;
    auto& degree = subgraph.degree;
    auto& graph = subgraph.graph;

    if (k == 2) {
        for (auto i : nodes[k]) {
            for (int j=0;j<degree[k][i];j++) {
                vector<int> tmp = prev;
                tmp.push_back(i);
                tmp.push_back(graph[i][j]);
                output.emplace_back(move(tmp));
            }
        }
    } else {
        for (auto i : nodes[k]) {
            if (Prune(subgraph, degree[k], i) < k-1) {
                continue;
            }

            BuildGraph(subgraph, k, i);
            prev.push_back(i);
            FindClique(subgraph, k-1, prev, output);
            prev.pop_back();
        }
    }
}

void CliqueList::BuildGraph(Subgraph& subgraph, int k, int u) {
    auto& nodes = subgraph.nodes;
    auto& degree = subgraph.degree;
    auto& graph = subgraph.graph;
    auto& candidate = subgraph.candicate;

    nodes[k-1].clear();

    for (int i=0;i<degree[k][u];i++) {
        int v = graph[u][i];
        candidate[v] = 1;
        nodes[k-1].push_back(v);
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

void CliqueList::MakeSubGraph(Subgraph& subgraph, int k, int u) {
    auto& nodes = subgraph.nodes;
    auto& degree = subgraph.degree;
    // auto& graph = subgraph.graph;
    auto& candidate = subgraph.candicate;

    nodes[k-1].clear();

    for (int i=0;i<degree[k][u];i++) {
        int v = graph[u][i];
        candidate[v] = 1;
        nodes[k-1].push_back(v);
    }

    for (int i=0;i<degree[k][u];i++) {
        int v = graph[u][i];

        subgraph.graph[v].clear();

        for (int j=0;j<degree[k][v];j++) {
            if (candidate[graph[v][j]] == 1) {
                subgraph.graph[v].push_back(graph[v][j]);
            }
        }

        degree[k-1][v] = subgraph.graph[v].size();
    }

    for (int i=0;i<degree[k][u];i++) {
        int v = graph[u][i];
        candidate[v] = 0;
    }
}

void CliqueList::InitColor() {
    vector<int> degeneracyOrder;
    degeneracyOrder.reserve(n);

    vector<vector<int>> g(n);
    for (int u=0;u<n;u++) {
        for (auto v : graph[u]) {
            g[u].emplace_back(v);
            g[v].emplace_back(u);
        }
    }

    vector<int> degree(n, 0);
    vector<int> visit(n, 0);
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    for (int i=0;i<n;i++) {
        degree[i] = g[i].size();
        pq.push(make_pair(g[i].size(), i));
    }

    while (!pq.empty()) {
        int id = pq.top().second;
        pq.pop();
        if (visit[id] == 0) {
            visit[id] = 1;
            degeneracyOrder.emplace_back(id);
            for (auto v : g[id]) {
                if (degree[v] > 0) {
                    pq.push(make_pair(--degree[v], v));
                }
            }
        }
    }

    color.resize(n, -1);
    for (auto it=degeneracyOrder.rbegin();it!=degeneracyOrder.rend();it++) {
        int i = *it;
        set<int> usedColor;
        for (auto v : g[i]) {
            if (color[v] != -1) {
                usedColor.insert(color[v]);
            }
            for (int j=0;j<n;j++) {
                if (usedColor.find(j) == usedColor.end()) {
                    color[i] = j;
                    break;
                }
            }
        }
    }
}

int CliqueList::UniqueColor(vector<int>& degree, int u) {
    unordered_set<int> uniqueColor;
    for (int i=0;i<degree[u];i++) {
        uniqueColor.emplace(color[graph[u][i]]);
    }
    return uniqueColor.size();
}

int CliqueList::UniqueColor(vector<vector<int>>& graph, vector<int>& degree, int u) {
    unordered_set<int> uniqueColor;
    for (int i=0;i<degree[u];i++) {
        uniqueColor.emplace(color[graph[u][i]]);
    }
    return uniqueColor.size();
}

int CliqueList::Prune(vector<int>& degree, int u) {
#ifndef USE_COLOR
    return degree[u];
#else
    return UniqueColor(degree, u);
#endif
}

int CliqueList::Prune(Subgraph& subgraph, vector<int>& degree, int u) {
#ifndef USE_COLOR
    return degree[u];
#else
    return UniqueColor(subgraph.graph, degree, u);
#endif
}


void CliqueList::FindCliqueCount(int k) {
    if (k == 2) {
        for (auto i : nodes[k]) {
            cliqueNum += degree[k][i];
        }
    } else {
        for (auto i : nodes[k]) {
            if (Prune(degree[k], i) < k-1) {
                continue;
            }

            BuildGraph(k, i);
            FindCliqueCount(k-1);
        }
    }
}

void CliqueList::FindClique(Subgraph& subgraph, int k, unsigned long long& num) {
    auto& nodes = subgraph.nodes;
    auto& degree = subgraph.degree;
    auto& graph = subgraph.graph;

    if (k == 2) {
        for (auto i : nodes[k]) {
            num += degree[k][i];
        }
    } else {
        for (auto i : nodes[k]) {
            if (Prune(subgraph, degree[k], i) < k-1) {
                continue;
            }

            BuildGraph(subgraph, k, i);
            FindClique(subgraph, k-1, num);
        }
    }
}

void CliqueList::FindCliqueCountScore(int k, vector<int>& prev) {
    if (k == 2) {
        for (auto i : nodes[k]) {
            for (int j=0;j<degree[k][i];j++) {
                for (auto p : prev) {
                    score[p]++;
                }
                score[i]++;
                score[graph[i][j]]++;
            }
        }
    } else {
        for (auto i : nodes[k]) {
            if (Prune(degree[k], i) < k-1) {
                continue;
            }

            BuildGraph(k, i);
            prev.push_back(i);
            FindCliqueCountScore(k-1, prev);
            prev.pop_back();
        }
    }
}

void CliqueList::FindCliqueCountScore(Subgraph& subgraph, int k, vector<int>& prev) {
    auto& nodes = subgraph.nodes;
    auto& degree = subgraph.degree;
    auto& graph = subgraph.graph;
    auto& score = subgraph.score;

    if (k == 2) {
        int edgeNum = 0;
        for (auto i : nodes[k]) {
            edgeNum += degree[k][i];
            score[i] += degree[k][i];
            for (int j=0;j<degree[k][i];j++) {
                score[graph[i][j]]++;
            }
        }
        for (auto p : prev) {
            score[p] += edgeNum;
        }
    } else {
        for (auto i : nodes[k]) {
            if (Prune(subgraph, degree[k], i) < k-1) {
                continue;
            }

            BuildGraph(subgraph, k, i);
            prev.push_back(i);
            FindCliqueCountScore(subgraph, k-1, prev);
            prev.pop_back();
        }
    }
}

void CliqueList::FindCliqueCountScoreMulti(Subgraph& subgraph, int k, vector<int>& prev, vector<int>& score) {
    auto& nodes = subgraph.nodes;
    auto& degree = subgraph.degree;
    auto& graph = subgraph.graph;
    // auto& score = subgraph.score;

    if (k == 2) {
        for (auto i : nodes[k]) {
            for (int j=0;j<degree[k][i];j++) {
                for (auto p : prev) {
                    #pragma omp atomic
                    score[p]++;
                }
                #pragma omp atomic
                score[i]++;

                #pragma omp atomic
                score[graph[i][j]]++;
            }
        }
    } else {
        for (auto i : nodes[k]) {
            if (Prune(subgraph, degree[k], i) < k-1) {
                continue;
            }

            BuildGraph(subgraph, k, i);
            prev.push_back(i);
            FindCliqueCountScore(subgraph, k-1, prev);
            prev.pop_back();
        }
    }
}