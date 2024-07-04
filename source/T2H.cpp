#include "T2H.h"

bool cmp_t2_increasing(std::pair<std::pair<int, int>, int> i, std::pair<std::pair<int, int>, int> j) {
    return i.second < j.second;
}

bool cmp_t2_decreasing(std::pair<std::pair<int, int>, int> i, std::pair<std::pair<int, int>, int> j) {
    return i.second > j.second;
}

bool cmp_degree(std::pair<int, long long> i, std::pair<int, long long> j) {
    return i.second > j.second;
}

unsigned long long T2HIndex::size() {
    unsigned long long num_intervals = 0;
    for (int u = 0; u < L_out.size(); u++) {
        for (auto s1 : L_out[u]) {
            for (int i = 0; i <= k; i++) {
                num_intervals += s1[i].size();
            }
        }
    }
    for (int u = 0; u < L_in.size(); u++) {
        for (auto s1 : L_in[u]) {
            for (int i = 0; i <= k; i++) {
                num_intervals += s1[i].size();
            }
        }
    }
    return num_intervals;
}

double T2HIndex::max_number_of_paths() {
    unsigned long long res = 0;
    for (int u = 0; u < L_out.size(); u++) {
        unsigned long long tmp1 = 0;
        for (auto s1 : L_out[u]) {
            unsigned long long tmp2 = 0;
            for (int i = 0; i <= k; i++) {
                tmp2 += s1[i].size();
            }
            tmp1 = std::max(tmp1, tmp2);
        }
        for (auto s1 : L_in[u]) {
            unsigned long long tmp2 = 0;
            for (int i = 0; i <= k; i++) {
                tmp2 += s1[i].size();
            }
            tmp1 = std::max(tmp1, tmp2);
        }
        res += tmp1;
    }
    return double(res) / L_out.size();
}

int T2HIndex::find_index(std::vector<int> &L_neighbours, int u) {
    int l = 0;
    int r = L_neighbours.size() - 1;
    while (l < r) {
        int mid = l + r >> 1;
        if (order[L_neighbours[mid]] < order[u]) {
            l = mid + 1;
        }
        else {
            r = mid;
        }
    }
    if (L_neighbours.size() > 0 && L_neighbours[l] == u) {
        return l;
    }
    else {
        return -1;
    }
}

std::pair<int, int>& T2HIndex::binary_search_ts_Lin(int u, int v_idx, int k, int ts) {
    int l = 0;
    int r = L_in[u][v_idx][k].size() - 1;
    if (r == -1) {
        return null_interval;
    }
    while (l < r) {
        int mid = (l + r) / 2;
        if (L_in[u][v_idx][k][mid].first < ts) {
            l = mid + 1;
        }
        else {
            r = mid;
        }
    }
    return L_in[u][v_idx][k][l];
}

std::pair<int, int>& T2HIndex::binary_search_te_Lin(int u, int v_idx, int k, int te) {
    int l = 0;
    int r = L_in[u][v_idx][k].size() - 1;
    if (r == -1) {
        return null_interval;
    }
    while (l < r) {
        int mid = (l + r + 1) / 2;
        if (L_in[u][v_idx][k][mid].second > te) {
            r = mid - 1;
        }
        else {
            l = mid;
        }
    }
    return L_in[u][v_idx][k][l];
}

std::pair<int, int>& T2HIndex::binary_search_ts_Lout(int u, int v_idx, int k, int ts) {
    int l = 0;
    int r = L_out[u][v_idx][k].size() - 1;
    if (r == -1) {
        return null_interval;
    }
    while (l < r) {
        int mid = (l + r + 1) / 2;
        if (L_out[u][v_idx][k][mid].first < ts) {
            r = mid - 1;
        }
        else {
            l = mid;
        }
    }
    return L_out[u][v_idx][k][l];
}

std::pair<int, int>& T2HIndex::binary_search_te_Lout(int u, int v_idx, int k, int te) {
    int l = 0;
    int r = L_out[u][v_idx][k].size() - 1;
    if (r == -1) {
        return null_interval;
    }
    while (l < r) {
        int mid = (l + r) / 2;
        if (L_out[u][v_idx][k][mid].second > te) {
            l = mid + 1;
        }
        else {
            r = mid;
        }
    }
    return L_out[u][v_idx][k][l];
}

bool T2HIndex::reachable(int u, int v, int ts, int te, int k) {
    if (u == v) {
        return true;
    }

    int v_idx = find_index(L_out_neighbours[u], v);
    if (v_idx != -1) {
        for (int i = 1; i <= k; i++) {
            auto interval = binary_search_ts_Lout(u, v_idx, i, ts);
            if (interval.first >= ts && interval.second <= te) {
                return true;
            }
        }
    }

    int u_idx = find_index(L_in_neighbours[v], u);
    if (u_idx != -1) {
        for (int i = 1; i <= k; i++) {
            auto interval = binary_search_ts_Lin(v, u_idx, i, ts);
            if (interval.first >= ts && interval.second <= te) {
                return true;
            }
        }
    }

    int i = 0, j = 0;
    while (i < L_out_neighbours[u].size() && j < L_in_neighbours[v].size()) {
        while (j < L_in_neighbours[v].size() && order[L_in_neighbours[v][j]] < order[L_out_neighbours[u][i]]) {
            j++;
        }
        if (j < L_in_neighbours[v].size() && L_in_neighbours[v][j] == L_out_neighbours[u][i]) {
            int d1 = k - 1, d2 = 1, ts_max = -1;
            while (d1 > 0) {
                auto interval1 = binary_search_te_Lin(v, j, d2, te);
                if (interval1.first >= ts && interval1.second <= te) {
                    ts_max = std::max(ts_max, interval1.first);
                }

                auto interval2 = binary_search_ts_Lout(u, i, d1, ts);
                if (interval2.first >= ts && interval2.second <= ts_max) {
                    return true;
                }
                d1--, d2++;
            }
        }
        i++;
    }

    return false;
}

void T2HIndex::construct_for_a_vertex(TemporalGraph* G, int u, bool reverse) {
    auto& L = reverse ? L_out : L_in;
    auto& L_neighbours = reverse ? L_out_neighbours : L_in_neighbours;

    std::vector<std::pair<std::pair<int, int>, int>> edge_set;
    std::unordered_set<int> Q[2];
    Q[0].insert(u);

    for (int d = 0; d < k; d++) {
        edge_set.clear();
        for (int v : Q[d % 2]) {
            auto e = G->getHeadEdge(v);
            if (reverse) {
                e = G->getHeadInEdge(v);
            }
            while (e) {
                if (order[e->to] > order[u]) {
                    edge_set.push_back(std::make_pair(std::make_pair(v, e->to), e->interaction_time));
                }
                e = G->getNextEdge(e);
            }
        }

        std::sort(edge_set.begin(), edge_set.end(), reverse ? cmp_t2_decreasing : cmp_t2_increasing);

        for (auto e : edge_set) {
            int v = e.first.first, w = e.first.second, t = e.second;
            if (u == v) {
                if (L_neighbours[w].size() == 0 || L_neighbours[w][L[w].size() - 1] != u) {
                    L[w].push_back(std::vector<std::vector<std::pair<int, int>>>());
                    L[w][L[w].size() - 1].resize(k + 1);
                    L_neighbours[w].push_back(u);
                    affected_vertices.insert(w);
                }
                auto it = minimal_intervals[w].begin();
                if (!reverse) {
                    it = minimal_intervals[w].lower_bound(std::make_pair(-t, 0));
                    if (it != minimal_intervals[w].end() && it->first == -t && it->second == t) {
                        continue;
                    }
                }
                else {
                    it = minimal_intervals[w].lower_bound(std::make_pair(t, 0)); 
                    if (it != minimal_intervals[w].end() && it->first == t && it->second == t) {
                        continue;
                    }
                }
                Q[(d + 1) % 2].insert(w);
                L[w][L[w].size() - 1][1].push_back(std::make_pair(t, t));
                if (!reverse) {
                    minimal_intervals[w].insert(std::make_pair(-t, t));
                }
                else {
                    minimal_intervals[w].insert(std::make_pair(t, t));
                }
                continue;
            }
            auto interval = reverse ? binary_search_ts_Lout(v, L[v].size() - 1, d, t) : binary_search_te_Lin(v, L[v].size() - 1, d, t);
            if (interval.first >= 0 && interval.second >= 0 && ((!reverse && interval.second <= t) || (reverse && interval.first >= t))) {
                int ts = interval.first, te = t;
                if (reverse) {
                    ts = t, te = interval.second;
                }
                auto it = minimal_intervals[w].begin();
                if (!reverse) {
                    it = minimal_intervals[w].lower_bound(std::make_pair(-te, 0));
                    if (it != minimal_intervals[w].end() && it->second >= ts) {
                        continue;
                    }
                }
                else {
                    it = minimal_intervals[w].lower_bound(std::make_pair(ts, 0));
                    if (it != minimal_intervals[w].end() && it->second <= te) {
                        continue;
                    }
                }
                if (it != minimal_intervals[w].end()) {
                    affected_vertices.insert(w);
                    if (!reverse) {
                        if (it->second == te) {
                            minimal_intervals[w].erase(it);
                        }
                        minimal_intervals[w].insert(std::make_pair(-te, ts));
                    }
                    else {
                        if (it->first == ts) {
                            minimal_intervals[w].erase(it);
                        }
                        minimal_intervals[w].insert(std::make_pair(ts, te));
                    }
                }
                if ((!reverse && !reachable(u, w, ts, te, d + 1)) || (reverse && !reachable(w, u, ts, te, d + 1))) {
                    Q[(d + 1) % 2].insert(w);
                    if (L_neighbours[w].size() == 0 || L_neighbours[w][L[w].size() - 1] != u) {
                        L[w].push_back(std::vector<std::vector<std::pair<int, int>>>());
                        L[w][L[w].size() - 1].resize(k + 1);
                        L_neighbours[w].push_back(u);
                    }
                    if (!reverse) {
                        auto& interval1 = binary_search_te_Lin(w, L[w].size() - 1, d + 1, te);
                        if (interval1.second == te) {
                            interval1.first = ts;
                            continue;
                        }
                    }
                    else {
                        auto& interval1 = binary_search_ts_Lout(w, L[w].size() - 1, d + 1, ts);
                        if (interval1.first == ts) {
                            interval1.second = te;
                            continue;
                        }
                    }
                    L[w][L[w].size() - 1][d + 1].push_back(std::make_pair(ts, te));
                }
            }
        }

        Q[d % 2].clear();
    }

    for (int v : affected_vertices) {
        minimal_intervals[v].clear();
    }
    affected_vertices.clear();
}

T2HIndex::T2HIndex(TemporalGraph* G, int k_input) {
    k = k_input;
    L_in.resize(G->n);
    L_out.resize(G->n);
    L_in_neighbours.resize(G->n);
    L_out_neighbours.resize(G->n);
    minimal_intervals.resize(G->n);

    std::vector<std::pair<int, long long>> vertex_set;
    for (int u = 0; u < G->n; u++) {
        vertex_set.push_back(std::make_pair(u, ((long long)G->in_degree[u] + 1) * (G->degree[u] + 1)));
    }
    std::sort(vertex_set.begin(), vertex_set.end(), cmp_degree);
    order = new int[G->n];
    for (int i = 0; i < vertex_set.size(); i++) {
        order[vertex_set[i].first] = i;
    }
    std::cout << "Vertex ordering completed. Start constructing the index..." << std::endl;

    unsigned long long start_time = currentTime();
    int i = 0;
    for (auto it = vertex_set.begin(); it != vertex_set.end(); it++) {
        int u = it->first;
        construct_for_a_vertex(G, u, false);
        construct_for_a_vertex(G, u, true);
        putProcess(double(++i) / G->n, currentTime() - start_time);
    }
}

T2HIndex::~T2HIndex() {
    delete [] order;
}

void T2HIndex::solve(TemporalGraph* G, char* query_file, char* output_file) {
    int s, t, ts, te, k;
    int query_num = 0;
    std::ifstream fin(query_file);
    std::ofstream fout(output_file);

    while (fin >> s >> t >> ts >> te >> k) {
        ++query_num;
    }

    fin = std::ifstream(query_file);

    int i = 0;
    unsigned long long start_time = currentTime();
    while (fin >> s >> t >> ts >> te >> k) {
        if (reachable(s, t, ts, te, k)) {
            fout << "Reachable" << std::endl;
        }
        else {
            fout << "Not reachable" << std::endl;
        }
        putProcess(double(++i) / query_num, currentTime() - start_time);
    }

    std::cout << "Average: " << timeFormatting((currentTime() - start_time) / query_num).str() << std::endl;
}