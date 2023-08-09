#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <queue>
#include <stack>
#include <map>
#include <set>
#include <ctime>
#include <chrono>
#include <random>
#include <omp.h>
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>
#include <lemon/maps.h>
#include "ilcplex/ilocplex.h"
using namespace std;
using namespace lemon;


using Graph = ListDigraph;
using Node = ListDigraph::Node;
using Edge = ListDigraph::Arc;

auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
std::mt19937 mt(seed);


/* hard coded values for height, width, power consumption, etc. */
float WIDTH = 2.448, HEIGHT = 1.2;
float RATIOS[6] = {1, 0.875, 0.854, 0.854, 0.844, 0.844}, POWER[6] = {1, 0.875, 0.854, 0.854, 0.844, 0.844};
int BITCNT[6] = {1, 4, 8, 16, 32, 64};


struct Point {
		float x = 0, y = 0;
};

struct Tray {
		Point pt;
		vector<Point> slots;
		int cand[70];
};

struct Flop {
		Point pt;
		int idx = 0;
		float prob = 0;

		bool operator < (const Flop& a) const {
				return prob < a.prob;
		}
};

struct Path {
		int a = 0, b = 0;
};


void setIO(string fin, string fout) {
	ios_base::sync_with_stdio(0); cin.tie(0);
	freopen(fin.c_str(),"r",stdin);
	freopen(fout.c_str(),"w",stdout);
}


float GetDist(Point a, Point b) {
		return (abs(a.x - b.x) + abs(a.y - b.y));
}

int GetRows(int K) {
		if (K == 1 || K == 2 || K == 4) return 1;
		if (K == 8) return 2;
		else return 4;
}


vector<Point> GetSlots(Point tray, int rows, int cols) {
		int bit_idx = 0;
		for (int i = 1; i < 6; i++) {
				if (rows * cols == BITCNT[i]) bit_idx = i;
		}

		float center_x = tray.x, center_y = tray.y;

		vector<Point> slots;
		for (int i = 0; i < rows * cols; i++) {
				int new_col = i % cols, new_row = i / cols;

				Point new_slot;
				new_slot.x = center_x + WIDTH * RATIOS[bit_idx] * ((new_col + 0.5) - (cols / 2.0));
				new_slot.y = center_y + HEIGHT * ((new_row + 0.5) - (rows / 2.0));


				if (new_slot.x >= 0 && new_slot.y >= 0) slots.push_back(new_slot);
		}

		return slots;
}


Flop GetNewFlop(vector<Flop> prob_dist, float tot_dist) {
		float rand_num = (float)(mt() % 101), cum_sum = 0;

		Flop new_flop;
		for (int i = 0; i < (int)prob_dist.size(); i++) {
				cum_sum += prob_dist[i].prob;
				new_flop = prob_dist[i];
				if (cum_sum * 100.0 >= rand_num * tot_dist) break;
		}

		return new_flop;
}

vector<Tray> GetStartTrays(vector<Flop>& flops, int num_trays) {

		int num_flops = (int)flops.size();

		/* pick a random flop */
		int rand_idx = mt() % (num_flops);
		Tray tray_zero; tray_zero.pt = flops[rand_idx].pt;

		set<int> used_flops; used_flops.insert(rand_idx);
		vector<Tray> trays; trays.push_back(tray_zero);
	
		float tot_dist = 0;
		for (int i = 0; i < num_flops; i++) {
				float contr = GetDist(flops[i].pt, tray_zero.pt);
				flops[i].prob = contr, tot_dist += contr;
		}

		while ((int)trays.size() < num_trays) {

				vector<Flop> prob_dist;
				for (int i = 0; i < num_flops; i++) if (!used_flops.count(i)) {
						prob_dist.push_back(flops[i]);
				}
				sort(begin(prob_dist), end(prob_dist));

				Flop new_flop = GetNewFlop(prob_dist, tot_dist);
				used_flops.insert(new_flop.idx);

				Tray new_tray; new_tray.pt = new_flop.pt;
				trays.push_back(new_tray);

				for (int i = 0; i < num_flops; i++) {
						float new_contr = GetDist(flops[i].pt, new_tray.pt);
						flops[i].prob += new_contr, tot_dist += new_contr;
				}

		}

		return trays;
}

void PrintTrays(vector<Tray> trays) {
		cout << "#trays: " << (int)trays.size() << "\n\n";
		for (int i = 0; i < (int)trays.size(); i++) {
				cout << "Center: (" << trays[i].pt.x << ", " << trays[i].pt.y << ")\n\n";

				cout << "#slots: " << (int)trays[i].slots.size() << "\n";
				for (int j = 0; j < (int)trays[i].slots.size(); j++) {
						cout << "(" << trays[i].slots[j].x << ", " << trays[i].slots[j].y << ")\n";
				}

				cout << "\n\n";
		}
}


Tray GetOneBit(Point pt) {
		Tray tray; 
		tray.pt.x = pt.x - WIDTH / 2.0, tray.pt.y = pt.y - HEIGHT / 2.0;
		tray.slots.push_back(pt);
		return tray;
}

vector<pair<int, int> > MinCostFlow(vector<Flop> flops, vector<Tray> &trays, int sz) {
		
		int num_flops = (int)flops.size(), num_trays = (int)trays.size();



		Graph graph; 
		vector<Node> nodes; vector<Edge> edges;
		vector<int> cur_costs, cur_caps;

		Node src = graph.addNode(), sink = graph.addNode();
		for (int i = 0; i < num_flops; i++) {
				Node flop_node = graph.addNode();
				nodes.push_back(flop_node);

				Edge src_to_flop = graph.addArc(src, flop_node);
				edges.push_back(src_to_flop);
				cur_costs.push_back(0), cur_caps.push_back(1);
		}


	 	vector<pair<int, int> > slot_to_tray;

		for (int i = 0; i < num_trays; i++) {
				vector<Point> tray_slots = trays[i].slots; 
				//cout << "#slots in tray: " << (int)tray_slots.size() << "\n";
				for (int j = 0; j < (int)tray_slots.size(); j++) {

						Node slot_node = graph.addNode();
						nodes.push_back(slot_node);

						for (int k = 0; k < num_flops; k++) {
								Edge flop_to_slot = graph.addArc(nodes[k], slot_node);
								edges.push_back(flop_to_slot);

								int edge_cost = 100 * (int)GetDist(flops[k].pt, tray_slots[j]);
								cur_costs.push_back(edge_cost), cur_caps.push_back(1);
						}

						Edge slot_to_sink = graph.addArc(slot_node, sink);
						edges.push_back(slot_to_sink);
						cur_costs.push_back(0), cur_caps.push_back(1);

						slot_to_tray.push_back(make_pair(i, j));
				}

		}

		Graph::ArcMap<int> costs(graph), caps(graph), flow(graph);
		Graph::NodeMap<int> labels(graph);
		NetworkSimplex<Graph, int, int> new_graph(graph);

		for (int i = 0; i < (int)edges.size(); i++) {
				costs[edges[i]] = cur_costs[i], caps[edges[i]] = cur_caps[i];
		}

		labels[src] = -1;
		for (int i = 0; i < (int)nodes.size(); i++) {
			labels[nodes[i]] = i;
		}
		labels[sink] = (int)nodes.size();

		new_graph.costMap(costs);
		new_graph.upperMap(caps);
		new_graph.stSupply(src, sink, num_flops);
		new_graph.run();
		new_graph.flowMap(flow);


		vector<pair<int, int> > clusters(num_flops);
		for (Graph::ArcIt itr(graph); itr != INVALID; ++itr) {
				int u = labels[graph.source(itr)], v = labels[graph.target(itr)];
				if (flow[itr] == 1 && u < num_flops && v >= num_flops) {
						v -= num_flops; 
						int tray_idx = slot_to_tray[v].first, slot_idx = slot_to_tray[v].second;
						clusters[u] = make_pair(tray_idx, slot_idx), trays[tray_idx].cand[slot_idx] = u;
				}
		}

		return clusters;
}

vector<Tray> RunLP(vector<Flop> flops, vector<Tray> trays, vector<pair<int, int> > clusters, int sz) {
		int num_flops = (int)flops.size(), num_trays = (int)trays.size();

		IloEnv my_env;
		IloModel my_model(my_env);

	
		IloNumVar tray_x[num_trays], tray_y[num_trays];
		for (int i = 0; i < num_trays; i++) {
				tray_x[i] = IloNumVar(my_env, 0, 1000000, ILOFLOAT);
				tray_y[i] = IloNumVar(my_env, 0, 1000000, ILOFLOAT);
		}


		IloNumVar d_x[num_flops], d_y[num_flops];
		for (int i = 0; i < num_flops; i++) {	

				d_x[i] = IloNumVar(my_env, 0, 1000000, ILOFLOAT);
				d_y[i] = IloNumVar(my_env, 0, 1000000, ILOFLOAT);

				int tray_idx = clusters[i].first, slot_idx = clusters[i].second;

				Point flop = flops[i].pt;
				Point tray = trays[tray_idx].pt;
				Point slot = trays[tray_idx].slots[slot_idx]; 

				/* 

				trying to minimize the sum of (|slot.x - flop.x| + |slot.y - flop.y|) over all flops
				for this LP, each flop uses the same slot shift as before:
					the new slot.x = new tray.x + shift in x for previous slot position
					the new slot.y = new tray.y + shift in y for previous slot position

				we already have the shifts (can derive from previous slot locations)

				*/

				float shift_x = slot.x - tray.x, shift_y = slot.y - tray.x;



				/* 
			
				let d_x + d_y = d

				|tray.x + shift_x - x| = d_x
				==> 
				|tray.x + shift_x - x| <= d_x
					
				so ...

				tray.x + shift_x - x <= d_x ==> shift_x - x <= d_x - tray.x
				and
				-tray.x - shift_x + x <= d_x ==> x - shift_x <= d_x + tray.x

				repeat for _y ! 

				*/


				my_model.add(shift_x - flop.x <= d_x[i] - tray_x[tray_idx]);
				my_model.add(flop.x - shift_x <= d_x[i] + tray_x[tray_idx]);


				my_model.add(shift_y - flop.y <= d_y[i] - tray_y[tray_idx]);
				my_model.add(flop.y - shift_y <= d_y[i] + tray_y[tray_idx]);
		}


		IloExpr obj(my_env);
		for (int i = 0; i < num_flops; i++) {
				obj += (d_x[i] + d_y[i]);
		}
		my_model.add(IloMinimize(my_env, obj));
		obj.end();

		IloCplex my_cplex(my_env);
		my_cplex.setOut(my_env.getNullStream());
		my_cplex.setParam(IloCplex::Threads, 1);
		my_cplex.setParam(IloCplex::EpGap, 0.01);
		my_cplex.setParam(IloCplex::TiLim, 18000);
		my_cplex.extract(my_model);
		my_cplex.solve();
		my_cplex.getStatus();

		float D = 0;
		for (int i = 0; i < num_flops; i++) {
				D += (my_cplex.getValue(d_x[i]) + my_cplex.getValue(d_y[i]));
		}
		//cout << "Objective Value: " << D << "\n";

		vector<Tray> new_trays;
		for (int i = 0; i < num_trays; i++) {
				Tray new_tray;
				new_tray.pt.x = my_cplex.getValue(tray_x[i]);
				new_tray.pt.y = my_cplex.getValue(tray_y[i]);
				new_trays.push_back(new_tray);
		}
		
		return new_trays;
}

void RunILP(vector<Flop> flops, vector<Path> paths, vector<vector<Tray> > all_trays, float ALPHA, float BETA) {
        /* convert all_trays into a single vectors */
        vector<Tray> trays;
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < (int)all_trays[i].size(); j++) trays.push_back(all_trays[i][j]);
        }

		int num_flops = (int)flops.size(), num_paths = (int)paths.size(), num_trays = (int)trays.size();
		
		IloEnv my_env;
		IloModel my_model(my_env);

		vector<IloNumVar> B[num_flops][num_trays];
		vector<IloNumVar> d_x[num_flops][num_trays], d_y[num_flops][num_trays];

		for (int i = 0; i < num_flops; i++) {
				for (int j = 0; j < num_trays; j++) {
					d_x[i][j].resize((int)trays[j].slots.size()), d_y[i][j].resize((int)trays[j].slots.size()), B[i][j].resize((int)trays[j].slots.size());
						for (int k = 0; k < (int)trays[j].slots.size(); k++) if (trays[j].cand[k] == i) {
								IloNumVar a = IloNumVar(my_env, 0, 10000000, ILOFLOAT);
								IloNumVar b = IloNumVar(my_env, 0, 10000000, ILOFLOAT);
								IloNumVar c = IloNumVar(my_env, 0, 1, ILOINT);
								d_x[i][j][k] = a, d_y[i][j][k] = b;
								B[i][j][k] = c;
						}
				}
		}

		/* add constraints for D */
		float maxD = 0;
		for (int i = 0; i < num_flops; i++) {
				for (int j = 0; j < num_trays; j++) {

						vector<Point> slots = trays[j].slots;


						for (int k = 0; k < (int)slots.size(); k++) if (trays[j].cand[k] == i) {
								
								float shift_x = slots[k].x - flops[i].pt.x, shift_y = slots[k].y - flops[i].pt.y;
								maxD = max(maxD, max(shift_x, -shift_x) + max(shift_y, -shift_y));


								my_model.add(B[i][j][k] >= 0);
								my_model.add(B[i][j][k] <= 1);

								my_model.add(d_x[i][j][k] - shift_x * B[i][j][k] >= 0);
								my_model.add(d_x[i][j][k] + shift_x * B[i][j][k] >= 0);

	
								my_model.add(d_y[i][j][k] - shift_y * B[i][j][k] >= 0);
								my_model.add(d_y[i][j][k] + shift_y * B[i][j][k] >= 0);
								
						}
				}
		} 


		IloExpr D(my_env);
		for (int i = 0; i < num_flops; i++) {
				for (int j = 0; j < num_trays; j++) {
						for (int k = 0; k < (int)trays[j].slots.size(); k++) if (trays[j].cand[k] == i) {
								my_model.add(d_x[i][j][k] + d_y[i][j][k] <= maxD);
								D += d_x[i][j][k] + d_y[i][j][k];
						}
				}
		}


	
		vector<IloNumVar> z_x(num_paths), z_y(num_paths);

		for (int i = 0; i < num_paths; i++) {
				int flop_a_idx = paths[i].a, flop_b_idx = paths[i].b;


				IloNumVar a = IloNumVar(my_env, 0, 1000000000, ILOFLOAT);
				IloNumVar b = IloNumVar(my_env, 0, 1000000000, ILOFLOAT);
				z_x[i] = a, z_y[i] = b;

				IloExpr d_x_a(my_env), d_y_a(my_env);
				IloExpr d_x_b(my_env), d_y_b(my_env);

				for (int j = 0; j < num_trays; j++) {
                        vector<Point> slots = trays[j].slots;
						for (int k = 0; k < (int)slots.size(); k++) {
								if (trays[j].cand[k] == flop_a_idx) {
                                        float shift_x = slots[k].x - flops[flop_a_idx].pt.x, shift_y = slots[k].y - flops[flop_a_idx].pt.y;
										d_x_a += (shift_x * B[flop_a_idx][j][k]);
										d_y_a += (shift_y * B[flop_a_idx][j][k]);
								}
								if (trays[j].cand[k] == flop_b_idx) {
                                        float shift_x = slots[k].x - flops[flop_b_idx].pt.x, shift_y = slots[k].y - flops[flop_b_idx].pt.y;
										d_x_b += (shift_x * B[flop_b_idx][j][k]);
										d_y_b += (shift_y * B[flop_b_idx][j][k]);
								}
						}
				}


				my_model.add(z_x[i] + d_x_a - d_x_b >= 0);
				my_model.add(z_x[i] + d_x_b - d_x_a >= 0);

				my_model.add(z_y[i] + d_y_a - d_y_b >= 0);
				my_model.add(z_y[i] + d_y_b - d_y_a >= 0);
		}

		IloExpr Z(my_env);
		for (int i = 0; i < num_paths; i++) {
				my_model.add(z_x[i] + z_y[i] <= maxD);
				Z += (z_x[i] + z_y[i]);
		}

		


		/* add constraints for B */
		for (int i = 0; i < num_trays; i++) {
				for (int j = 0; j < (int)trays[i].slots.size(); j++) {
						IloExpr b_slot(my_env);
						for (int k = 0; k < num_flops; k++) if (trays[i].cand[j] == k) {
								b_slot += B[k][i][j];
						}
						my_model.add(b_slot >= 0);
						my_model.add(b_slot <= 1);
						b_slot.clear();
				}
		}

		for (int i = 0; i < num_flops; i++) {
				IloExpr b_flop(my_env);
				for (int j = 0; j < num_trays; j++) {
						for (int k = 0; k < (int)trays[j].slots.size(); k++) if (trays[j].cand[k] == i) {
								b_flop += B[i][j][k];
						}
				}
				my_model.add(b_flop == 1);
				b_flop.clear();
		}

		/* add constraints for W */
		vector<IloNumVar> e(num_trays);
		for (int i = 0; i < num_trays; i++) {
				e[i] = IloNumVar(my_env, 0, 1, ILOINT);
				IloExpr slots_used(my_env);
				for (int j = 0; j < (int)trays[i].slots.size(); j++) {
						for (int k = 0; k < num_flops; k++) if (trays[i].cand[j] == k) {
								my_model.add(e[i] - B[k][i][j] >= 0);
								slots_used += B[k][i][j];
						}
				}
				my_model.add(slots_used - e[i] >= 0);
				my_model.add(e[i] >= 0);
				my_model.add(e[i] <= 1);
				slots_used.clear();
		}

		vector<float> w(num_trays);
		for (int i = 0; i < num_trays; i++) {
				int bit_idx = 0;
				for (int j = 0; j < 6; j++) {
						if (BITCNT[j] == (int)trays[i].slots.size()) bit_idx = j;
				}
				if (BITCNT[bit_idx] == 1) w[i] = 1.00;
				else w[i] = ((float)BITCNT[bit_idx]) * (POWER[bit_idx] - BITCNT[bit_idx] * 0.0015);
		}

		IloExpr W(my_env);
		for (int i = 0; i < num_trays; i++) {
				W += (w[i] * e[i]);
		}



		IloExpr obj(my_env);
		obj = D + ALPHA * W + BETA * Z;
		my_model.add(IloMinimize(my_env, obj));

		IloCplex my_cplex(my_env);
		my_cplex.setParam(IloCplex::Threads, 20);
		my_cplex.setParam(IloCplex::EpGap, 0.01);
		my_cplex.setParam(IloCplex::TiLim, 18000);
		my_cplex.setOut(my_env.getNullStream());

		my_cplex.extract(my_model);
		my_cplex.solve();
		my_cplex.getStatus();


		cout << "Total: " << my_cplex.getObjValue() << "\n";
		cout << "D: " << my_cplex.getValue(D) << "\n";
		cout << "Z: " << my_cplex.getValue(Z) << "\n";
		cout << "W: " << my_cplex.getValue(W) << "\n";

		
		set<int> tray_idx;
		for (int i = 0; i < num_flops; i++) {
				for (int j = 0; j < num_trays; j++) {
						for (int k = 0; k < (int)trays[j].slots.size(); k++) if (trays[j].cand[k] == i) {
								if (my_cplex.getValue(B[i][j][k]) == 1) {
										tray_idx.insert(j);
								} 
						}
				}
		}


		map<int, int> tray_sz;
		for (auto x: tray_idx) {
				tray_sz[(int)trays[x].slots.size()]++;
		}

		for (auto x: tray_sz) {
				cout << "Tray size = " << x.first << ": " << x.second << "\n";
		}


		D.clear(), Z.clear(), W.clear();

}


int main(int argc, char *argv[]) {
		float ALPHA = stof(argv[1]), BETA = stof(argv[2]);
		setIO(argv[3], argv[4]);

		int num_flops, num_paths;
		cin >> num_flops >> num_paths;

		vector<Flop> flops(num_flops);
		for (int i = 0; i < num_flops; i++) {
				float x, y; cin >> x >> y; 

				Point new_pt; 
				new_pt.x = x, new_pt.y = y;

				Flop new_flop; 
				new_flop.pt = new_pt, new_flop.idx = i;
				flops[i] = new_flop;
		}

		vector<Path> paths(num_paths);
		for (int i = 0; i < num_paths; i++) {
				int a, b; cin >> a >> b;

				Path new_path;
				new_path.a = a, new_path.b = b;
				paths[i] = new_path;
		}


		time_t start, end; time(&start);

		vector<vector<Tray> > trays(6);
		vector<vector<pair<int, int> > > clusters(6);

        for (int i = 0; i < 6; i++) {
                int num_trays = (num_flops + (BITCNT[i] - 1)) / BITCNT[i];
                trays[i].resize(num_trays);
                clusters[i].resize(num_flops);
        }

		for (int i = 0; i < num_flops; i++) {
				Tray one_bit = GetOneBit(flops[i].pt);
				one_bit.cand[0] = i;
				trays[0][i] = one_bit;
		}


		omp_set_num_threads(5);

		#pragma omp parallel for
		for (int i = 1; i < 6; i++) {

				int rows = GetRows(BITCNT[i]), cols = BITCNT[i] / rows;
				float AR = (cols * WIDTH * RATIOS[i]) / (rows * HEIGHT);

				int num_trays = (num_flops + (BITCNT[i] - 1)) / BITCNT[i];
				trays[i] = GetStartTrays(flops, num_trays);
				for (int j = 0; j < num_trays; j++) trays[i][j].slots = GetSlots(trays[i][j].pt, rows, cols);

				// clusters[i] = (tray that flop i belongs to, slot that flop i belongs to)
				for (int j = 0; j < 15; j++) {
						clusters[i] = MinCostFlow(flops, trays[i], BITCNT[i]);
						#pragma omp critical 
						{
							trays[i] = RunLP(flops, trays[i], clusters[i], BITCNT[i]);
						}
						for (int k = 0; k < num_trays; k++) {
							trays[i][k].slots = GetSlots(trays[i][k].pt, rows, cols);
						}
				}
				clusters[i] = MinCostFlow(flops, trays[i], BITCNT[i]);

				for (int j = 0; j < num_trays; j++) {
					trays[i][j].slots = GetSlots(trays[i][j].pt, rows, cols);
				}
		}


        clock_t ilp_start, ilp_end; ilp_start = clock();

		RunILP(flops, paths, trays, ALPHA, BETA);

        ilp_end = clock();
        cout << "ILP time: " << fixed <<  double(ilp_end - ilp_start) / double(CLOCKS_PER_SEC) << setprecision(5) << " seconds\n";


		time(&end);
		cout << "Total time: " << fixed <<  double(end - start) << setprecision(5) << " seconds\n";
}










