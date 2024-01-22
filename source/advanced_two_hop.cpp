#include "advanced_two_hop.h"

bool cmp(std::pair<int, long long> i, std::pair<int, long long> j) {
    return i.second > j.second;
}

bool cmp1(std::vector<int> i, std::vector<int> j) {
    return i[1] > j[1] || (i[1] == j[1] && i[2] < j[2]);
}

int AdvancedTwoHopIndex::size() {
    int num_intervals = 0;
    for (auto it = L_out.begin(); it != L_out.end(); it++) {
        num_intervals += it->size();
    }
    for (auto it = L_in.begin(); it != L_in.end(); it++) {
        num_intervals += it->size();
    }
    return num_intervals;
}

int AdvancedTwoHopIndex::find_lowest_index(std::vector<std::vector<int>> &L, int v) {
    int l = 0;
    int r = L.size() - 1;
    while (l < r) {
        int mid = l + r >> 1;
        if (order[L[mid][0]] < order[v]) {
            l = mid + 1;
        }
        else {
            r = mid;
        }
    }
    return r;
}

bool AdvancedTwoHopIndex::reachable(TemporalGraph* G, int u, int v, int ts, int te, int d) {
    if (u == v) {
        return true;
    }

    int idx = find_lowest_index(L_out[u], v);
    if (idx != -1) {
        while (L_out[u][idx][0] == v) {
            if (L_out[u][idx][3] > d) {
                break;
            }
            int l = idx;
            int r = cut_out[u][idx] - 1;
            while (l < r) {
                int mid = l + r >> 1;
                if (L_out[u][mid][1] >= ts) {
                    r = mid;
                }
                else {
                    l = mid + 1;
                }
            }
            if (L_out[u][l][1] >= ts && L_out[u][l][2] <= te) {
                return true;
            }
            idx = cut_out[u][idx];
        }
    }

    idx = find_lowest_index(L_in[v], u);
    if (idx != -1) {
        while (L_in[v][idx][0] == u) {
            if (L_in[v][idx][3] > d) {
                break;
            }
            int l = idx;
            int r = cut_in[v][idx] - 1;
            while (l < r) {
                int mid = l + r >> 1;
                if (L_in[v][mid][1] >= ts) {
                    r = mid;
                }
                else {
                    l = mid + 1;
                }
            }
            if (L_in[v][l][1] >= ts && L_in[v][l][2] <= te) {
                return true;
            }
            idx = cut_in[v][idx];
        }
    }

    int i = 0;
    int j = 0;
    while (i < L_out[u].size() && j < L_in[v].size()) {
        if (order[L_out[u][i][0]] < order[L_in[v][j][0]]) {
            if (i >= next_out[u].size()) {
                return false;
            }
            i = next_out[u][i];
        }
        else if (order[L_out[u][i][0]] > order[L_in[v][j][0]]) {
            if (j >= next_in[v].size()) {
                return false;
            }
            j = next_in[v][j];
        }
        else {
            int w = L_out[u][i][0];
            int d1 = d + 1;
            int t1 = te + 1;
            while (i < L_out[u].size() && L_out[u][i][0] == w) {
                if (L_out[u][i][3] > d) {
                    if (i >= next_out[u].size()) {
                        return false;
                    }
                    i = next_out[u][i];
                    break;
                }
                int l = i;
                int r = cut_out[u][i] - 1;
                while (l < r) {
                    int mid = l + r + 1 >> 1;
                    if (L_out[u][mid][1] < ts) {
                        r = mid - 1;
                    }
                    else {
                        l = mid;
                    }
                }
                if (L_out[u][l][1] >= ts && L_out[u][l][2] <= te) {
                    d1 = L_out[u][l][3];
                    t1 = L_out[u][l][2];
                    break;
                }
                i = cut_out[u][i];
            }

            if (d1 <= d) {
                while (j < L_in[v].size() && L_in[v][j][0] == w) {
                    if (L_in[v][j][3] > d - d1) {
                        if (j >= next_in[v].size()) {
                            return false;
                        }
                        j = next_in[v][j];
                        break;
                    }
                    int l = j;
                    int r = cut_in[v][j] - 1;
                    while (l < r) {
                        int mid = l + r + 1 >> 1;
                        if ((!is_temporal_path && L_in[v][mid][1] < ts) || (is_temporal_path && L_in[v][mid][1] < t1)) {
                            r = mid - 1;
                        }
                        else {
                            l = mid;
                        }
                    }
                    if (((!is_temporal_path && L_in[v][l][1] >= ts) || (is_temporal_path && L_in[v][l][1] >= t1)) && L_in[v][l][2] <= te) {
                        return true;
                    }
                    j = cut_in[v][j];
                }
            }

            if (i >= next_out[u].size()) {
                i = L_out[u].size();
            }
            else {
                i = next_out[u][i];
            }
        }
    }

    return false;
}

void AdvancedTwoHopIndex::construct_for_a_vertex(TemporalGraph* G, int u, bool reverse, std::vector<std::vector<std::vector<int>>> &L, 
                                                std::vector<std::vector<int>> &next_idx, 
                                                std::vector<std::vector<int>> &cut_idx,
                                                int t_threshold) {
    for (auto u : affected_vertices) {
        inc_index[u].clear();
        binary_indexed_tree[u] = std::vector<std::pair<int, int>>();
    }
    affected_vertices.clear();
    Q.push(std::vector<int>{u, G->tmax + 1, -1, 0});

    while (!Q.empty()) {
        std::vector<int> current = Q.front();
        Q.pop();
        int v = current[0], ts = current[1], te = current[2], d = current[3];
        if (u == v && d > 0) {
            continue;
        }

        // Check minimality to avoid expanding non-minimal paths
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
        }
        if (flag) {
            continue;
        }

        // Check whether u can reach v in ([ts, te], d)
        if (u != v && ((!reverse && reachable(G, u, v, ts, te, d)) || (reverse && reachable(G, v, u, ts, te, d)))) {
            continue;
        }

        // Insert the path into the index
        if (u != v) {
            affected_vertices.insert(v);
            while (next_idx[v].size() < L[v].size()) {
                next_idx[v].push_back(L[v].size());
            }
            inc_index[v].push_back(std::vector<int>{u, ts, te, d});
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
            if (!is_temporal_path || (!reverse && te <= e->interaction_time) || (reverse && ts >= e->interaction_time)) {
                int ts_new = std::min(ts, e->interaction_time);
                int te_new = std::max(te, e->interaction_time);
                if (t_threshold != -1 && te_new - ts_new + 1 > t_threshold) {
                    e = e->next;
                    continue;
                }
                visited_paths++;
                bool flag = false;
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
                if (!flag) {
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
                    Q.push(std::vector<int>{e->to, ts_new, te_new, d + 1});
                }
            }
            e = e->next;
        }
    }

    // Organize all paths in an increasing ts order
    for (auto v : affected_vertices) {
        int i = 0, j = 0;
        for (i = 0; i <= inc_index[v].size(); i++) {
            if (i == inc_index[v].size() || (i > 0 && inc_index[v][i][3] != inc_index[v][i - 1][3])) {
                std::sort(inc_index[v].begin() + j, inc_index[v].begin() + i, cmp1);
                int cur_tmin = G->tmax + 1;
                while (j < i) {
                    if (inc_index[v][j][2] < cur_tmin) {
                        cur_tmin = inc_index[v][j][2];
                        L[v].push_back(inc_index[v][j]);
                    }
                    j++;
                }
                while (cut_idx[v].size() < L[v].size()) {
                    cut_idx[v].push_back(L[v].size());
                }
            }
        }
    }
}

AdvancedTwoHopIndex::AdvancedTwoHopIndex(TemporalGraph* G, int k_input, int t_threshold, std::string path_type) {
    k = k_input;
    L_in.resize(G->n);
    L_out.resize(G->n);
    next_in.resize(G->n);
    next_out.resize(G->n);
    cut_in.resize(G->n);
    cut_out.resize(G->n);
    inc_index.resize(G->n);
    binary_indexed_tree.resize(G->n);

    if (path_type == "Temporal") {
        is_temporal_path = true;
    }
    else {
        is_temporal_path = false;
    }

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
        construct_for_a_vertex(G, u, false, L_in, next_in, cut_in, t_threshold);
        if (G->is_directed) {
            construct_for_a_vertex(G, u, true, L_out, next_out, cut_out, t_threshold);
        }
        else {
            L_out = L_in;
            next_out = next_in;
            cut_out = cut_in;
        }
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