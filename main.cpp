#include <bits/stdc++.h>
#include "mbff.h"
using namespace std;

int main() {
    ifstream fin("test.in");

    int num_flops;
    vector<float> x, y;

    fin >> num_flops;
    x.resize(num_flops), y.resize(num_flops);

    for (int i = 0; i < num_flops; i++) {
        int idx;
        float a, b;
        fin >> idx >> a >> b;

        x[idx] = a, y[idx] = b;
    }

    int num_paths;
    vector<pair<int, int> > paths;

    fin >> num_paths;
    paths.resize(num_paths);

    for (int i = 0; i < num_paths; i++) {
        int a, b;
        fin >> a >> b;
        paths[i] = make_pair(a, b);
    }

    gpl::MBFF inst = gpl::MBFF(num_flops, num_paths, x, y, paths, 25);

    inst.Run(500, 20.0, 1.0);
}
