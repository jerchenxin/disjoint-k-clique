// g++ -O3 -std=c++17 greedyAllclique.cpp CliqueList.cpp -fopenmp -o greedyAllclique 
// g++ -O3 -std=c++17 greedyAllclique.cpp CliqueList.cpp -fopenmp -o greedyAllclique -DPRINT_CLIQUE


#include "CliqueList.h"

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

int k;
vector<int> valid;

bool check(vector<int>& c) {
    for (int i=0;i<k;i++) {
        if (!valid[c[i]]) return false;
    }
    return true;
}

int main(int argc, char* argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    string fileName(argv[1]);
    k = atoi(argv[2]);

    CliqueList cl(fileName, k);
    cl.ListMultiThread(64);

    auto endList = std::chrono::high_resolution_clock::now();
    auto diffList = std::chrono::duration_cast<std::chrono::nanoseconds>(endList - start);

#ifndef PRINT_CLIQUE
    cout << "List Time(ms): " << diffList.count() / 1000000 << endl;

    cout << "n: " << cl.n << endl;
    cout << "K-Clique Number: " << cl.output.size() << endl;
#endif

    auto& cliqueList = cl.output;
    int n = cl.n;
    valid.resize(n, 1);
    vector<int> degree(n, 0);
    for (auto& c : cliqueList) {
        for (int i=0;i<k;i++) {
            degree[c[i]]++;
        }
    }

    for (auto& c : cliqueList) {
        int score = 0;
        for (int i=0;i<k;i++) {
            score += degree[c[i]];
        }
        c.emplace_back(score);
    }

    sort(cliqueList.begin(), cliqueList.end(), [](auto& a, auto& b){return a[k] < b[k];});

    vector<vector<int>> output;

    for (auto& c : cliqueList) {
        if (check(c)) {
            for (int i=0;i<k;i++) {
                valid[c[i]] = 0;
            }
            c.pop_back();
            output.emplace_back(c);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
#ifndef PRINT_CLIQUE
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

