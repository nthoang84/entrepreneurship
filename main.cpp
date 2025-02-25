#include "cagetti_denardi.h"

#include <iostream>
#include <chrono>

using namespace std;

int main() {
    auto start = chrono::high_resolution_clock::now();
    CagettiDeNardi model;
    model.computeAssetGrid();
    model.computeIncomeInvDist();
    model.computePolicy(0.065);
    model.simulate();
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    cout << ">> Elapsed time: " << duration.count() << " seconds" << endl;
    return 0;
}