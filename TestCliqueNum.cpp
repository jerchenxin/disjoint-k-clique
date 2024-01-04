// g++ -O3 -std=c++17 TestCliqueNum.cpp CliqueList.cpp -fopenmp -o TestCliqueNum

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

    for (int k=3;k<=7;k++) {
        CliqueList cl(fileName, k);

        auto start = std::chrono::high_resolution_clock::now();
        cl.CountMulti(64);
        auto end = std::chrono::high_resolution_clock::now();
        auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        cout << k << "-Count Time(ms): " << diff.count() / 1000000 << endl;
        cout << "Clique Num: " << cl.cliqueNum << endl << endl;
    }
}