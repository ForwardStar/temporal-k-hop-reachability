#include "advanced_two_hop.h"

bool cmp(std::pair<int, long long> i, std::pair<int, long long> j) {
    return i.second > j.second;
}

bool cmp1(std::vector<int> i, std::vector<int> j) {
    return i[1] > j[1] || (i[1] == j[1] && i[2] < j[2]);
}

int AdvancedTwoHopIndex::size() {
    int num_intervals = 0;
    for (int u = 0; u < L_out.size(); u++) {
        for (auto s : L_out[u]) {
            for (int d = 0; d < k; d++) {
                num_intervals += s.second[d].size();
            }
        }
    }
    for (int u = 0; u < L_in.size(); u++) {
        for (auto s : L_in[u]) {
            for (int d = 0; d < k; d++) {
                num_intervals += s.second[d].size();
            }
        }
    }
    return num_intervals;
}

bool AdvancedTwoHopIndex::reachable(TemporalGraph* G, int u, int v, int ts, int te, int k) {
    if (u == v) {
        return true;
    }

    if (L_out[u].find(v) != L_out[u].end()) {
        for (int d = 0; d < k; d++) {
            int l = 0;
            int r = L_out[u][v][d].size() - 1;
            if (l > r) {
                continue;
            }
            while (l < r) {
                int mid = l + r + 1 >> 1;
                if (L_out[u][v][d][mid][1] < ts) {
                    r = mid - 1;
                }
                else {
                    l = mid;
                }
            }
            if (L_out[u][v][d][l][1] >= ts && L_out[u][v][d][l][2] <= te) {
                return true;
            }
        }
    }

    if (L_in[v].find(u) != L_in[v].end()) {
        for (int d = 0; d < k; d++) {
            int l = 0;
            int r = L_in[v][u][d].size() - 1;
            if (l > r) {
                continue;
            }
            while (l < r) {
                int mid = l + r + 1 >> 1;
                if (L_in[v][u][d][mid][1] < ts) {
                    r = mid - 1;
                }
                else {
                    l = mid;
                }
            }
            if (L_in[v][u][d][l][1] >= ts && L_in[v][u][d][l][2] <= te) {
                return true;
            }
        }
    }

    for (auto s : L_out[u]) {
        int w = s.first;
        if (L_in[v].find(w) != L_in[v].end()) {
            for (int d1 = 0; d1 < k - 1; d1++) {
                int l = 0;
                int r = L_out[u][w][d1].size() - 1;
                if (l > r) {
                    continue;
                }
                while (l < r) {
                    int mid = l + r + 1 >> 1;
                    if (L_out[u][w][d1][mid][1] < ts) {
                        r = mid - 1;
                    }
                    else {
                        l = mid;
                    }
                }
                if (L_out[u][w][d1][l][1] >= ts && L_out[u][w][d1][l][2] <= te) {
                    int t1 = L_out[u][w][d1][l][2];
                    for (int d2 = 0; d2 <= k - d1 - 2; d2++) {
                        int l = 0;
                        int r = L_in[v][w][d2].size() - 1;
                        if (l > r) {
                            continue;
                        }
                        while (l < r) {
                            int mid = l + r + 1 >> 1;
                            if ((!is_temporal_path && L_in[v][w][d2][mid][1] < ts) || (is_temporal_path && L_in[v][w][d2][mid][1] < t1)) {
                                r = mid - 1;
                            }
                            else {
                                l = mid;
                            }
                        }
                        if (((!is_temporal_path && L_in[v][w][d2][l][1] >= ts) || (is_temporal_path && L_in[v][w][d2][l][1] >= t1)) && L_in[v][w][d2][l][2] <= te) {
                            return true;
                        }
                    }
                }
            }
        }
    }

    return false;
}

void AdvancedTwoHopIndex::construct_for_a_vertex(TemporalGraph* G, int u, bool reverse, std::vector<std::unordered_map<int, std::vector<std::vector<std::vector<int>>>>> &L, int t_threshold) {
    for (auto u : affected_vertices) {
        temp_paths[u].clear();
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
            temp_paths[v].push_back(std::vector<int>{u, ts, te, d});
            affected_vertices.insert(v);
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
        L[v][u] = std::vector<std::vector<std::vector<int>>>();
        L[v][u].resize(k);
        int i = 0, j = 0;
        for (i = 0; i <= temp_paths[v].size(); i++) {
            if (i == temp_paths[v].size() || (i > 0 && temp_paths[v][i][3] != temp_paths[v][i - 1][3])) {
                std::sort(temp_paths[v].begin() + j, temp_paths[v].begin() + i, cmp1);
                int d = temp_paths[v][j][3] - 1;
                int cur_tmin = G->tmax + 1;
                while (j < i) {
                    if (temp_paths[v][j][2] < cur_tmin) {
                        cur_tmin = temp_paths[v][j][2];
                        L[v][u][d].push_back(temp_paths[v][j]);
                    }
                    j++;
                }
            }
        }
    }
}

AdvancedTwoHopIndex::AdvancedTwoHopIndex(TemporalGraph* G, int k_input, int t_threshold, std::string path_type) {
    k = k_input;
    L_in.resize(G->n);
    L_out.resize(G->n);
    temp_paths.resize(G->n);
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
        construct_for_a_vertex(G, u, false, L_in, t_threshold);
        construct_for_a_vertex(G, u, true, L_out, t_threshold);
        putProcess(double(++i) / G->n, currentTime() - start_time);
    }
}

AdvancedTwoHopIndex::~AdvancedTwoHopIndex() {
    delete [] order;
}

void AdvancedTwoHopIndex::solve(TemporalGraph* G, char* query_file, char* output_file) {
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