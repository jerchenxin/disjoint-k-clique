// g++ -O3 -std=c++17 BuildCliqueGraph.cpp -o BuildCliqueGraph

#include <algorithm>
#include <unordered_set>
#include <set>
#include <fstream>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <random>
#include <cstdlib>
#include <iostream>
#include <functional>
#include <chrono>

using namespace std;

struct Clique {
    unsigned long long id;
    vector<int> nodes;
};

bool Overlap(Clique* a, Clique* b) {
    auto it_a = a->nodes.begin();
    auto it_b = b->nodes.begin();

    while (it_a != a->nodes.end() && it_b != b->nodes.end()) {
        if (*it_a < *it_b) {
            it_a++;
        } else if (*it_a > *it_b) {
            it_b++;
        } else {
            return true;
        }
    }

    return false;
}

int main(int argc, char* argv[]) {
    string fileName(argv[1]);
    int k = atoi(argv[2]);
    unsigned long long m = 0;

    vector<vector<Clique*>> cliqueGraph;

    FILE *f = nullptr;
    f = fopen(fileName.c_str(), "r");
    vector<int> tmp(k);
    while (fscanf(f, "%d", &tmp[0]) == 1) {
        for (int i=1;i<k;i++) {
            fscanf(f, "%d", &tmp[i]);
        }
        sort(tmp.begin(), tmp.end());
        if (cliqueGraph.size() < tmp.back() + 1) {
            cliqueGraph.resize(tmp.back() + 1);
        }
        Clique* c = new Clique();
        c->nodes = tmp;
        c->id = m++;
        for (auto i : c->nodes) {
            cliqueGraph[i].push_back(c);
        }
    }
    fclose(f);

    vector<pair<int, int>> g;

    auto start = std::chrono::high_resolution_clock::now();

    for (int i=0;i<cliqueGraph.size();i++) {
        for (auto c : cliqueGraph[i]) {
            if (c->nodes.front() == i) { // calculate once
                unordered_set<Clique*> cand;
                for (auto j : c->nodes) {
                    cand.insert(cliqueGraph[j].begin(), cliqueGraph[j].end());
                }
                cand.erase(c);
                for (auto ne : cand) {
                    if (ne->id < c->id && Overlap(ne, c)) {
                        g.emplace_back(c->id, ne->id);
                        // cout << c->id << " " << ne->id << endl;
                    }
                }
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    cout << "Time(ms): " << diff.count() / 1000000 << endl;

    string outputFileName(argv[3]);
    ofstream output(outputFileName);
    output << m << " " << g.size() << endl;

    for (auto [u, v] : g) {
        output << u << " " << v << endl;
    }
}