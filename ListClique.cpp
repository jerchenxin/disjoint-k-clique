// g++ -O3 -std=c++17 ListClique.cpp CliqueList.cpp -fopenmp -o ListClique

#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <algorithm>
#include <tuple>
#include <map>
#include <chrono>
#include <random>
#include <omp.h>

#include "CliqueList.h"

using namespace std;

int main(int argc, char *argv[]) {
    string fileName(argv[1]);

    int k = atoi(argv[2]);

    CliqueList cl(fileName, k);

    cl.ListMultiThread(64);
    
    string outputFileName(argv[3]);
    ofstream output(outputFileName);

    for (auto c : cl.output) {
        for (auto i : c) {
            output << i << " ";
        }
        output << endl;
    }
}