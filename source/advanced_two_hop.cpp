#include "advanced_two_hop.h"

bool cmp(std::pair<int, long long> i, std::pair<int, long long> j) {
    return i.second > j.second;
}

bool cmp1(std::vector<int> i, std::vector<int> j) {
    return i[0] > j[0] || (i[0] == j[0] && i[1] < j[1]);
}

unsigned long long AdvancedTwoHopIndex::size() {
    unsigned long long num_intervals = 0;
    for (int u = 0; u < L_out.size(); u++) {
        for (auto s1 : L_out[u]) {
            for (auto s2 : s1) {
                num_intervals += s2.second.size();
            }
        }
    }
    for (int u = 0; u < L_in.size(); u++) {
        for (auto s1 : L_in[u]) {
            for (auto s2 : s1) {
                num_intervals += s2.second.size();
            }
        }
    }
    return num_intervals;
}

int AdvancedTwoHopIndex::find_index(std::vector<int> &L_neighbours, int u) {
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

bool AdvancedTwoHopIndex::reachable(TemporalGraph* G, int u, int v, int ts, int te, int k) {
    if (u == v) {
        return true;
    }

    int idx = find_index(L_out_neighbours[u], v);
    if (idx != -1) {
        for (auto L : L_out[u][idx]) {
            if (L.first > k) {
                break;
            }
            int l = 0;
            int r = L.second.size() - 1;
            while (l < r) {
                int mid = l + r + 1 >> 1;
                if (L.second[mid].first < ts) {
                    r = mid - 1;
                }
                else {
                    l = mid;
                }
            }
            if (L.second[l].first >= ts && L.second[l].second <= te) {
                return true;
            }
        }
    }

    idx = find_index(L_in_neighbours[v], u);
    if (idx != -1) {
        for (auto L : L_in[v][idx]) {
            if (L.first > k) {
                break;
            }
            int l = 0;
            int r = L.second.size() - 1;
            while (l < r) {
                int mid = l + r + 1 >> 1;
                if (L.second[mid].first < ts) {
                    r = mid - 1;
                }
                else {
                    l = mid;
                }
            }
            if (L.second[l].first >= ts && L.second[l].second <= te) {
                return true;
            }
        }
    }

    int j = 0;
    for (int i = 0; i < L_out_neighbours[u].size(); i++) {
        int w = L_out_neighbours[u][i];
        while (j < L_in_neighbours[v].size() && order[L_in_neighbours[v][j]] < order[w]) {
            j++;
        }
        if (j >= L_in_neighbours[v].size()) {
            break;
        }
        if (w == L_in_neighbours[v][j]) {
            int ts_max = -1;
            int te_min = G->tmax + 1;
            int idx = 0;
            int d2 = L_in[v][j][idx].first;
            for (int idx2 = L_out[u][i].size() - 1; idx2 >= 0; idx2--) {
                int l = 0;
                int r = L_out[u][i][idx2].second.size() - 1;
                int d1 = L_out[u][i][idx2].first;
                if (d1 >= k) {
                    continue;
                }
                while (l < r) {
                    int mid = l + r + 1 >> 1;
                    if (L_out[u][i][idx2].second[mid].first < ts) {
                        r = mid - 1;
                    }
                    else {
                        l = mid;
                    }
                }
                if (L_out[u][i][idx2].second[l].first >= ts && L_out[u][i][idx2].second[l].second <= te) {
                    int t1 = L_out[u][i][idx2].second[l].second;
                    while (d1 + d2 <= k) {
                        int l = 0;
                        int r = L_in[v][j][idx].second.size() - 1;
                        while (l < r) {
                            int mid = l + r >> 1;
                            if (L_in[v][j][idx].second[mid].second > te) {
                                l = mid + 1;
                            }
                            else {
                                r = mid;
                            }
                        }
                        if (L_in[v][j][idx].second[l].first >= ts && L_in[v][j][idx].second[l].second <= te && L_in[v][j][idx].second[l].first >= ts_max) {
                            ts_max = L_in[v][j][idx].second[l].first;
                            te_min = L_in[v][j][idx].second[l].second;
                        }
                        idx++;
                        if (idx >= L_in[v][j].size()) {
                            d2 = k + 1;
                        }
                        else {
                            d2 = L_in[v][j][idx].first;
                        }
                    }
                    if (((!is_temporal_path && ts_max >= ts) || (is_temporal_path && ts_max >= t1)) && te_min <= te) {
                        return true;
                    }
                }
            }
        }
    }

    return false;
}

void AdvancedTwoHopIndex::construct_for_a_vertex(TemporalGraph* G, int u, bool reverse, 
                                                std::vector<std::vector<std::vector<std::pair<int, std::vector<std::pair<int, int>>>>>> &L,
                                                std::vector<std::vector<int>> &L_neighbours,
                                                int t_threshold) {
    for (auto u : affected_vertices) {
        temp_paths[u].clear();
        binary_indexed_tree[u].clear();
        temp_binary_indexed_tree[u].clear();
    }
    affected_vertices.clear();
    Q.push(std::vector<int>{u, G->tmax + 1, -1, 0});
    int cur_len = 0;
    std::vector<std::vector<int>> cur_paths;

    while (!Q.empty()) {
        std::vector<int> current = Q.front();
        Q.pop();
        int v = current[0], ts = current[1], te = current[2], d = current[3];
        if (u == v && d > 0) {
            continue;
        }

        if (d > cur_len) {
            cur_len = d;
            for (auto e : cur_paths) {
                int v = e[0];
                temp_binary_indexed_tree[v].clear();
                int ts_e = e[1];
                int te_e = e[2];
                if (binary_indexed_tree[v].size() == 0) {
                    binary_indexed_tree[v].assign(G->tmax + 2, G->tmax + 1);
                }
                int t = ts_e + 1;
                while (t) {
                    binary_indexed_tree[v][t] = std::min(binary_indexed_tree[v][t], te_e);
                    t -= (t & (-t));
                }
            }
            cur_paths.clear();
        }

        // Check minimality to avoid expanding non-minimal paths
        bool flag = false;
        if (binary_indexed_tree[v].size() > 0) {
            int t = ts + 1;
            while (t <= G->tmax + 1) {
                if (binary_indexed_tree[v][t] < te) {
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
            temp_paths[v].push_back(std::vector<int>{ts, te, d});
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
                deleted_edges.push_back(std::make_pair(e, reverse));
                e = G->deleteEdge(e);
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
                        if (binary_indexed_tree[e->to][t] <= te_new) {
                            flag = true;
                            break;
                        }
                        t += (t & (-t));
                    }
                }
                if (!flag) {
                    if (temp_binary_indexed_tree[e->to].size() > 0) {
                        int t = ts_new + 1;
                        while (t <= G->tmax + 1) {
                            if (temp_binary_indexed_tree[e->to][t] <= te_new) {
                                flag = true;
                                break;
                            }
                            t += (t & (-t));
                        }
                    }
                    if (!flag) {
                        if (temp_binary_indexed_tree[e->to].size() == 0) {
                            temp_binary_indexed_tree[e->to].assign(G->tmax + 2, G->tmax + 1);
                        }
                        int t = ts_new + 1;
                        while (t) {
                            temp_binary_indexed_tree[e->to][t] = std::min(temp_binary_indexed_tree[e->to][t], te_new);
                            t -= (t & (-t));
                        }
                        cur_paths.push_back(std::vector<int>{e->to, ts_new, te_new});
                        Q.push(std::vector<int>{e->to, ts_new, te_new, d + 1});
                    }
                }
            }
            e = e->next;
        }
    }

    // Organize all paths in an increasing ts order
    for (auto v : affected_vertices) {
        int idx = L[v].size();
        L_neighbours[v].push_back(u);
        L[v].push_back(std::vector<std::pair<int, std::vector<std::pair<int, int>>>>());
        int i = 0, j = 0;
        for (i = 0; i <= temp_paths[v].size(); i++) {
            if (i == temp_paths[v].size() || (i > 0 && temp_paths[v][i][2] != temp_paths[v][i - 1][2])) {
                std::sort(temp_paths[v].begin() + j, temp_paths[v].begin() + i, cmp1);
                int d = temp_paths[v][j][2];
                int idx2 = L[v][idx].size();
                L[v][idx].push_back(std::pair<int, std::vector<std::pair<int, int>>>());
                L[v][idx][idx2].first = d;
                int cur_tmin = G->tmax + 1;
                while (j < i) {
                    if (temp_paths[v][j][1] < cur_tmin) {
                        cur_tmin = temp_paths[v][j][1];
                        L[v][idx][idx2].second.push_back(std::make_pair(temp_paths[v][j][0], temp_paths[v][j][1]));
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
    L_in_neighbours.resize(G->n);
    L_out_neighbours.resize(G->n);
    temp_paths.resize(G->n);
    binary_indexed_tree.resize(G->n);
    temp_binary_indexed_tree.resize(G->n);

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
    std::cout << "Vertex ordering completed. Start constructing the index..." << std::endl;

    unsigned long long start_time = currentTime();
    int i = 0;
    for (auto it = vertex_set.begin(); it != vertex_set.end(); it++) {
        int u = it->first;
        construct_for_a_vertex(G, u, false, L_in, L_in_neighbours, t_threshold);
        construct_for_a_vertex(G, u, true, L_out, L_out_neighbours, t_threshold);
        for (int i = deleted_edges.size() - 1; i >= 0; i--) {
            auto e = deleted_edges[i].first;
            bool in_edge = deleted_edges[i].second; 
            if (e->last) {
                e->last->next = e;
            }
            else {
                if (in_edge) {
                    G->head_in_edge[e->from] = e;
                }
                else {
                    G->head_edge[e->from] = e;
                }
            }
            if (e->next) {
                e->next->last = e;
            }
        }
        deleted_edges.clear();
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