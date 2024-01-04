// g++ -O3 -std=c++17 bench.cpp dynamic.cpp lightPruneClass.cpp CliqueList.cpp -fopenmp -o bench

#include <chrono>

#include "dynamic.h"

using namespace std;

int main(int argc, char* argv[]) {
    string graphName(argv[1]);

    DynamicGraph dg("/home/xchen/cliqueDataset/" + graphName, atoi(argv[2]));

    vector<pair<int, int>> update;
    {
        ifstream file("/home/xchen/cliqueDataset/update/" + graphName + ".inc");
        int u, v;
        while (file >> u >> v) {
            update.emplace_back(u, v);
        }
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        for (auto [u, v] : update) {
            dg.Deletion(u, v);
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        cout << "delete Time: " << diff.count() / update.size() << endl;
    }

    dg.Info();

    {
        auto start = std::chrono::high_resolution_clock::now();

        for (auto [u, v] : update) {
            dg.Insertion(u, v);
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        cout << "add Time: " << diff.count() / update.size() << endl;
    }

    dg.Info();

    // dg.IndexCorrectness();
    
    // dg.Info();

    // for (auto [k, v] : dg.invGroup) {
    //     for (auto i : v) {
    //         cout << i << " ";
    //     }
    //     cout << endl;
    // }
}