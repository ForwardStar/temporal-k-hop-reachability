#include "advanced_two_hop.h"

bool cmp(std::pair<int, long long> i, std::pair<int, long long> j) {
    return i.second > j.second;
}

bool check(TemporalGraph* G, int u, int w, int ts, int te, int d, std::unordered_map<int, std::vector<std::vector<int>>> &L) {
    int t = G->tmax + 1;
    int ts_w = ts + 1;
    while (ts_w <= G->tmax + 1) {
        if (L[w][ts_w].size() > 0) {
            int d_w = d;
            while (d_w > 0) {
                t = std::min(t, L[w][ts_w][d_w]);
                if (t <= te) {
                    return true;
                }
                d_w -= (d_w & (-d_w));
            }
        }
        ts_w += (ts_w & (-ts_w));
    }
    return false;
}

bool AdvancedTwoHopIndex::reachable(TemporalGraph* G, int u, int v, int ts, int te, int d) {
    if (u == v) {
        return true;
    }
    if (order[v] < order[u]) {
        for (auto it = L_out[u].begin(); it != L_out[u].end(); it++) {
            int w = it->first;
            if (w != v && L_in[v].find(w) == L_in[v].end()) {
                continue;
            }
            int l = 1;
            int r = d + 1;
            while (l < r) {
                int mid = l + r >> 1;
                if (check(G, u, w, ts, te, mid, L_out[u])) {
                    r = mid;
                }
                else {
                    l = mid + 1;
                }
            }
            if (r == d + 1) {
                continue;
            }
            if (w == v) {
                return true;
            }
            if (d == l) {
                continue;
            }
            if (check(G, v, w, ts, te, d - l, L_in[v])) {
                return true;
            }
        }
    }
    else {
        for (auto it = L_in[v].begin(); it != L_in[v].end(); it++) {
            int w = it->first;
            if (w != u && L_out[u].find(w) == L_out[u].end()) {
                continue;
            }
            int l = 1;
            int r = d + 1;
            while (l < r) {
                int mid = l + r >> 1;
                if (check(G, v, w, ts, te, mid, L_in[v])) {
                    r = mid;
                }
                else {
                    l = mid + 1;
                }
            }
            if (r == d + 1) {
                continue;
            }
            if (w == u) {
                return true;
            }
            if (d == l) {
                continue;
            }
            if (check(G, u, w, ts, te, d - l, L_out[u])) {
                return true;
            }
        }
    }
    return false;
}

void AdvancedTwoHopIndex::construct_for_a_vertex(TemporalGraph* G, int u, bool reverse, std::vector<std::unordered_map<int, std::vector<std::vector<int>>>> &L, std::string algorithm) {
    std::vector<std::vector<std::pair<int, int>>> binary_indexed_tree;
    binary_indexed_tree.resize(G->n);
    std::queue<std::vector<int>> Q_BFS;
    auto smaller_interval_first = [](std::vector<int> i, std::vector<int> j) {
        return (i[2] - i[1] > j[2] - j[1]) || (i[2] - i[1] == j[2] - j[1] && i[3] > j[3]);
    };
    std::priority_queue<std::vector<int>, std::vector<std::vector<int>>, decltype(smaller_interval_first)> Q_PR(smaller_interval_first);
    std::vector<int> start;
    start.push_back(u);
    start.push_back(G->tmax + 1);
    start.push_back(-1);
    start.push_back(0);
    if (algorithm == "BFS") {
        Q_BFS.push(start);
    }
    else {
        Q_PR.push(start);
    }

    while (1) {
        std::vector<int> current;
        if (algorithm == "BFS") {
            if (Q_BFS.empty()) {
                break;
            }
            current = Q_BFS.front();
            Q_BFS.pop();
        }
        else {
            if (Q_PR.empty()) {
                break;
            }
            current = Q_PR.top();
            Q_PR.pop();
        }
        int v = current[0];
        int ts = current[1];
        int te = current[2];
        int d = current[3];
        if (u == v && d > 0) {
            continue;
        }

        // Dequeue checking
        if (algorithm == "BFS") {
            bool flag = false;
            if (binary_indexed_tree[v].size() > 0) {
                int t = ts + 1;
                while (t <= G->tmax + 1) {
                    if (binary_indexed_tree[v][t].first < te && binary_indexed_tree[v][t].second == d) {
                        flag = true;
                        break;
                    }
                    t += (t & (-t));
                }
                if (flag) {
                    continue;
                }
            }
            if (flag) {
                continue;
            }
        }

        // Check whether u can reach v in ([ts, te], d)
        if (u != v && ((!reverse && reachable(G, u, v, ts, te, d)) || (reverse && reachable(G, v, u, ts, te, d)))) {
            continue;
        }

        // Insert the path into the index
        if (u != v) {
            if (L[v].find(u) == L[v].end()) {
                L[v][u] = std::vector<std::vector<int>>();
                L[v][u].resize(G->tmax + 2);
            }
            int t1 = ts + 1;
            while (t1 > 0) {
                if (L[v][u][t1].size() == 0) {
                    L[v][u][t1].assign(k + 1, G->tmax + 1);
                }
                int t2 = d;
                while (t2 <= k) {
                    L[v][u][t1][t2] = std::min(te, L[v][u][t1][t2]);
                    t2 += (t2 & (-t2));
                }
                t1 -= (t1 & (-t1));
            }
        }

        if (d == k) {
            continue;
        }
        TemporalGraph::Edge* e = G->getHeadEdge(v);
        if (reverse) {
            e = G->getHeadInEdge(v);
        }
        while (e) {
            if (order[e->to] < order[u]) {
                e = e->next;
                continue;
            }
            int ts_new = std::min(ts, e->interaction_time);
            int te_new = std::max(te, e->interaction_time);
            bool flag = false;
            // enqueue checking
            if (algorithm == "BFS") {
                if (binary_indexed_tree[e->to].size() > 0) {
                    int t = ts_new + 1;
                    while (t <= G->tmax + 1) {
                        if (binary_indexed_tree[e->to][t].first <= te_new) {
                            flag = true;
                            break;
                        }
                        t += (t & (-t));
                    }
                }
            }
            if (!flag) {
                if (algorithm == "BFS") {
                    if (binary_indexed_tree[e->to].size() == 0) {
                        binary_indexed_tree[e->to].assign(G->tmax + 2, std::make_pair(G->tmax + 1, 0));
                    }
                    int t = ts_new + 1;
                    while (t > 0) {
                        if (te < binary_indexed_tree[e->to][t].first) {
                            binary_indexed_tree[e->to][t].first = te_new;
                            binary_indexed_tree[e->to][t].second = d + 1;
                        }
                        t -= (t & (-t));
                    }
                }
                std::vector<int> next;
                next.push_back(e->to);
                next.push_back(ts_new);
                next.push_back(te_new);
                next.push_back(d + 1);
                if (algorithm == "BFS") {
                    Q_BFS.push(next);
                }
                else {
                    Q_PR.push(next);
                }
            }
            e = e->next;
        }
    }
}

AdvancedTwoHopIndex::AdvancedTwoHopIndex(TemporalGraph* G, int k_input, int t_threshold, std::string algorithm) {
    k = k_input;
    index_construct_algorithm = algorithm;
    L_in.resize(G->n);
    L_out.resize(G->n);

    std::vector<std::pair<int, long long>> vertex_set;
    for (int u = 0; u < G->n; u++) {
        vertex_set.push_back(std::make_pair(u, ((long long)G->in_degree[u] + 1) * (G->degree[u] + 1)));
    }
    std::sort(vertex_set.begin(), vertex_set.end(), cmp);
    order = new int[G->n];
    for (int i = 0; i < vertex_set.size(); i++) {
        order[vertex_set[i].first] = i;
    }
    std::cout << "Vertex ordering completed. Start constructing index..." << std::endl;

    unsigned long long start_time = currentTime();
    int i = 0;
    for (auto it = vertex_set.begin(); it != vertex_set.end(); it++) {
        int u = it->first;
        construct_for_a_vertex(G, u, false, L_in, algorithm);
        construct_for_a_vertex(G, u, true, L_out, algorithm);
        putProcess(double(++i) / G->n, currentTime() - start_time);
    }
}

AdvancedTwoHopIndex::~AdvancedTwoHopIndex() {
    delete [] order;
}

void AdvancedTwoHopIndex::solve(TemporalGraph* G, char* query_file, char* output_file, int k) {
    int s, t, ts, te;
    int query_num = 0;
    std::ifstream fin(query_file);
    std::ofstream fout(output_file);

    while (fin >> s >> t >> ts >> te) {
        ++query_num;
    }

    fin = std::ifstream(query_file);

    int i = 0;
    unsigned long long start_time = currentTime();
    while (fin >> s >> t >> ts >> te) {
        // Perform online BFS Search
        if (reachable(G, s, t, ts, te, k)) {
            fout << "Reachable" << std::endl;
        }
        else {
            fout << "Not reachable" << std::endl;
        }
        putProcess(double(++i) / query_num, currentTime() - start_time);
    }

    std::cout << "Average: " << timeFormatting((currentTime() - start_time) / query_num).str() << std::endl;
}