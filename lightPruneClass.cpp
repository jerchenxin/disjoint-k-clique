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
#include <omp.h>

#include "lightPruneClass.h"

#define DENSITY 5

ThirdGreedy::ThirdGreedy(string fileName, int k) : k(k) {

}

void ThirdGreedy::solve(string fileName) {

    auto start = std::chrono::high_resolution_clock::now();
    CliqueList* cl = new CliqueList(fileName, k);
    // cl->CountScoreSingle();
    cl->CountScoreMulti(threadNum);
    score = move(cl->score);
    n = cl->n;

#ifndef PRINT_CLIQUE
    auto endScore = std::chrono::high_resolution_clock::now();
    auto diffScore = std::chrono::duration_cast<std::chrono::nanoseconds>(endScore - start);
    cout << "Score Time(ms): " << diffScore.count() / 1000000 << endl;
#endif

    degree.resize(k+1, vector<int>(n, 0));
    nodes.resize(k+1);
    candidate.resize(n, 0);
    valid.resize(n, 1);
    graph.resize(n);

    // calculate rank
    vector<pair<int, int>> rank;
    for (int i=0;i<n;i++) {
        rank.emplace_back(score[i], i);
    }
    sort(rank.begin(), rank.end());

    vertexRank.resize(n);
    orderedVertex.reserve(n);
    for (int i=0;i<n;i++) {
        int id = rank[i].second;
        orderedVertex.push_back(id);
        vertexRank[id] = i;
    }
    
    for (int u=0;u<n;u++) {
        for (int v : cl->graph[u]) {
            if (vertexRank[u] > vertexRank[v]) {
                graph[u].emplace_back(v);
            } else {
                graph[v].emplace_back(u);
            }
        }
    }

#ifndef PRINT_CLIQUE    
    cout << "n: " << n << endl;
#endif

    for (int i=0;i<n;i++) {
        degree[k][i] = graph[i].size();
    }

    delete cl;

    unsigned long long midTime = 0;

    unsigned long long maxNum = 0;

#ifndef PRINT_CLIQUE
    auto endPreprocess = std::chrono::high_resolution_clock::now();
    auto diffPreprocess = std::chrono::duration_cast<std::chrono::nanoseconds>(endPreprocess - start);
    cout << "Preprocess Time(ms): " << diffPreprocess.count() / 1000000 << endl;
#endif
    
    // the third
    {
        PQ pq;
        pq.elements.reserve(n);
        vector<int> redo(n, 0);

        omp_set_num_threads(threadNum);
        int i;
        #pragma omp parallel private(i)
        {
            int pid = omp_get_thread_num();
            Subgraph subgraph(degree, nodes, n);
            #pragma omp for schedule(dynamic, 1) nowait
            for (i=0;i<n;i++) {
                if (score[i] == 0 || subgraph.degree[k][i] < k-1) {
                    continue;
                }

                Clique* targetClique = new Clique();
                int scoreThreshold = INT32_MAX;

                MakeSubGraph(subgraph, k, i);
                vector<int> prev = {i};
                FindCliqueMulti(subgraph, k-1, prev, &targetClique, scoreThreshold, score[i]);
                if (scoreThreshold != INT32_MAX) {
                    #pragma omp critical
                    {
                        pq.push(make_pair(scoreThreshold, targetClique));
                    }
                } else {
                    delete targetClique;
                }
            }
        }

        #ifndef PRINT_CLIQUE
            auto endScore = std::chrono::high_resolution_clock::now();
            auto diffScore = std::chrono::duration_cast<std::chrono::nanoseconds>(endScore - start);
            cout << "Init Time(ms): " << diffScore.count() / 1000000 << endl;
        #endif

        while (!pq.empty()) {
            auto c = pq.top().second;
            if (check(c->nodes, k)) {
                for (auto i : c->nodes) {
                    DeleteNode(i);
                }
                output.emplace_back(c->nodes);
                delete c;   
                pq.pop();
            } else {
                int probeID = c->nodes[0];
                pq.pop();

                if (valid[probeID] == 0 || degree[k][probeID] < k-1) {
                    delete c;
                    continue;
                }

                Clique* targetClique = c;
                int scoreThreshold = INT32_MAX;

                redo[probeID]++;

                BuildGraph(k, probeID);
                vector<int> prev = {probeID};
                FindAllClique(k-1, prev, &targetClique, scoreThreshold, score[probeID]);

                if (scoreThreshold != INT32_MAX) {
                    pq.push(make_pair(scoreThreshold, targetClique));
                } else {
                    delete targetClique;
                }
            }
        }

#ifndef PRINT_CLIQUE
        int sum = 0;
        for (int i=0;i<n;i++) {
            sum += redo[i];
        }
        cout << "redo " << sum << endl;

        // for (auto i : orderedVertex) {
            // cout << redo[i] << endl;
        // }
#endif
    }

#ifdef PRINT_CLIQUE
    for (auto& c : output) {
        for (auto i : c) {
            cout << i << " ";
        }
        cout << endl;
    }
#endif

    auto end = std::chrono::high_resolution_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

#ifndef PRINT_CLIQUE
    cout << "Time(ms): " << diff.count() / 1000000 << endl;
    cout << "Mid Time1(ms): " << midTime / 1000000 << endl;
    cout << "Max Clique: " << maxNum << endl;
    
    cout << "Solution Number: " << output.size() << endl;
#endif
}

void ThirdGreedy::FindAllClique(int k, vector<int>& prev, Clique** targetClique, int& targetScore, int curScore) {
    if (k == 2) {
        for (auto i : nodes[k]) {
            int i_score = score[i] + curScore;

            if (i_score >= targetScore) {
                continue;
            }

            for (int j=0;j<degree[k][i];j++) {
                int j_id = graph[i][j];
                int cliqueScore = i_score + score[j_id];

                if (cliqueScore < targetScore) {
                    (*targetClique)->nodes = prev;
                    (*targetClique)->nodes.reserve(nodes.size()+2);
                    (*targetClique)->nodes.push_back(i);
                    (*targetClique)->nodes.push_back(j_id);
                    targetScore = cliqueScore;
                }
            }
        }
    } else {
        for (auto i : nodes[k]) {
            if (degree[k][i] < k-1 || curScore + score[i] >= targetScore) {
                continue;
            }

            BuildGraph(k, i);
            prev.push_back(i);
            FindAllClique(k-1, prev, targetClique, targetScore, curScore + score[i]);
            prev.pop_back();
        }
    }
}

void ThirdGreedy::BuildGraph(int k, int u) {
    nodes[k-1].clear();

    for (int i=0;i<degree[k][u];i++) {
        int v = graph[u][i];
        if (valid[v]) {
            candidate[v] = 1;
            nodes[k-1].push_back(v);
        } else {
            swap(graph[u][i], graph[u][degree[k][u]-1]); // move the invalid one into the end
            degree[k][u]--;
            i--;
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

    // reset
    for (int v : nodes[k-1]) {
        candidate[v] = 0;
    }
}

void ThirdGreedy::DeleteNode(int u) {
    valid[u] = 0;
}

bool ThirdGreedy::check(vector<int>& c, int k) {
    for (int i=0;i<k;i++) {
        if (!valid[c[i]]) {
            return false;
        }
    }
    return true;
}

void ThirdGreedy::MakeSubGraph(Subgraph& subgraph, int k, int u) {
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

void ThirdGreedy::BuildGraph(Subgraph& subgraph, int k, int u) {
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

void ThirdGreedy::FindCliqueMulti(Subgraph& subgraph, int k, vector<int>& prev, Clique** targetClique, int& targetScore, int curScore) {
    auto& nodes = subgraph.nodes;
    auto& degree = subgraph.degree;
    auto& graph = subgraph.graph;

    if (k == 2) {
        for (auto i : nodes[k]) {
            int i_score = score[i] + curScore;

            if (i_score >= targetScore) {
                continue;
            }

            for (int j=0;j<degree[k][i];j++) {
                int j_id = graph[i][j];
                int cliqueScore = i_score + score[j_id];

                if (cliqueScore < targetScore) {
                    (*targetClique)->nodes = prev;
                    (*targetClique)->nodes.reserve(nodes.size()+2);
                    (*targetClique)->nodes.push_back(i);
                    (*targetClique)->nodes.push_back(j_id);
                    targetScore = cliqueScore;
                }
            }
        }
    } else {
        for (auto i : nodes[k]) {
            if (degree[k][i] < k-1 || curScore + score[i] >= targetScore) {
                continue;
            }

            BuildGraph(subgraph, k, i);
            prev.push_back(i);
            FindCliqueMulti(subgraph, k-1, prev, targetClique, targetScore, curScore + score[i]);
            prev.pop_back();
        }
    }
}
