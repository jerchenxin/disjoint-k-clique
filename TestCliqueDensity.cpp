// g++ -O3 -std=c++17 TestCliqueDensity.cpp CliqueList.cpp -fopenmp -o TestCliqueDensity

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
        cl.CountScoreMulti(64);
        auto end = std::chrono::high_resolution_clock::now();
        auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        cout << k << "-Count Time(ms): " << diff.count() / 1000000 << endl;

        unsigned long long alpha = 0;
        for (int i=0;i<cl.n;i++) {
            alpha += cl.score[i] * (cl.score[i] - 1) / 2;
        }

        cout << "Alpha: " << alpha << endl << endl;
    }
}