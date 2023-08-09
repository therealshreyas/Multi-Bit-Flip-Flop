#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <set>
#include <stack>
#include <queue>
#include <random>
#include <chrono>
#include <map>
#include <ctime>	
using namespace std;

typedef pair<int, int> pii;

#define sz(x) (int)x.size()
#define f first
#define s second
#define mp make_pair
#define pb push_back

auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
std::mt19937 mt(seed);


int TESTS, BBOX, N, M;
vector<pair<float, float> > pointsets;


void generatePoints(int cardinality, int max) {
  int random = rand();
  int pointX = mt() % (BBOX * 100 + 1);
  int randomTwo = rand();
  int pointY = mt() % (BBOX * 100 + 1);
  //cout << pointX << " " << pointY << "\n";
  pair<float, float> coordinate;
  coordinate.first = (float)pointX/100.0 + 70;
  coordinate.second = (float)pointY/100.0 + 70;
  // Checks the see if a pair(coordinate) is in the vector pointsets. If it is,
  // it recurses back to the beginning of the function to generate a new point.
  if (!(find(pointsets.begin(), pointsets.end(), coordinate) !=
        pointsets.end())) {
    pointsets.push_back(coordinate);
  } else {
    return generatePoints(cardinality, max);
  }
}

int main(int argc, char *argv[]) {
  srand((int)time(0)); 
 
  BBOX = stoi(argv[2]);
  N = stoi(argv[4]);
  M = stoi(argv[6]);
  TESTS = stoi(argv[8]);

  // ofstream fout; fout.open("test.in");
  //cout << TESTS << "\n";
  for (int T = 1; T <= TESTS; T++) {
    for (int i = 0; i < N; i++) {
      generatePoints(N, M);
    }


    cout << (int)pointsets.size() << " " << M << "\n";

    for (int i = 0; i < pointsets.size(); i++) {
      cout << pointsets[i].first << " " << pointsets[i].second
           << "\n";
    }
 	
 	vector<pii> paths;
 	while (sz(paths) < M) {
 		int a = mt() % (N), b = mt() % (N);
 		if (a > b) swap(a, b);
 		if (!(find(paths.begin(), paths.end(), mp(a, b)) != paths.end()) && a != b) paths.push_back(mp(a, b));
 	}

 	for (int i = 0; i < sz(paths); i++) {
 		cout << paths[i].f << " " << paths[i].s << "\n";
 	}

    pointsets.clear();
  }


  return 0;
}




