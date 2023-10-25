#include "advanced_two_hop.h"

bool cmp(std::pair<int, long long> i, std::pair<int, long long> j) {
    return i.second > j.second;
}

bool cmp1(std::vector<int> i, std::vector<int> j) {
    return i[1] > j[1];
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

    if (index_construct_algorithm == "PrioritySearch" || index_construct_algorithm == "BFS-naive") {
        int idx = find_lowest_index(L_out[u], v);
        for (int i = idx; i < L_out[u].size(); i++) {
            if (L_out[u][i][0] != v) {
                break;
            }
            if (L_out[u][i][1] >= ts && L_out[u][i][2] <= te && L_out[u][i][3] <= d) {
                return true;
            }
        }

        idx = find_lowest_index(L_in[v], u);
        for (int i = idx; i < L_in[v].size(); i++) {
            if (L_in[v][i][0] != u) {
                break;
            }
            if (L_in[v][i][1] >= ts && L_in[v][i][2] <= te && L_in[v][i][3] <= d) {
                return true;
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
            else if (L_out[u][i][1] >= ts && L_out[u][i][2] <= te) {
                if (L_in[v][j][1] >= ts && L_in[v][j][2] <= te) {
                    if (L_out[u][i][3] + L_in[v][j][3] <= d) {
                        return true;
                    }
                    if (i >= next_out[u].size() || j >= next_in[v].size()) {
                        break;
                    }
                    i = next_out[u][i];
                    j = next_in[v][j];
                }
                else {
                    j++;
                }
            }
            else {
                i++;
            }
        }
    }
    else {
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
                int d1 = d;
                while (i < L_out[u].size() && L_out[u][i][0] == w) {
                    if (L_out[u][i][3] >= d) {
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
                        if (i >= next_out[u].size()) {
                            i = L_out[u].size();
                        }
                        else {
                            i = next_out[u][i];
                        }
                        break;
                    }
                    i = cut_out[u][i];
                }

                if (d1 < d) {
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
                            if (L_in[v][mid][1] < ts) {
                                r = mid - 1;
                            }
                            else {
                                l = mid;
                            }
                        }
                        if (L_in[v][l][1] >= ts && L_in[v][l][2] <= te) {
                            return true;
                        }
                        j = cut_in[v][j];
                    }
                }
            }
        }
    }

    return false;
}

void AdvancedTwoHopIndex::construct_for_a_vertex(TemporalGraph* G, int u, bool reverse, std::vector<std::vector<std::vector<int>>> &L, 
                                                std::vector<std::vector<int>> &next_idx, 
                                                std::vector<std::vector<int>> &cut_idx,
                                                int t_threshold,
                                                std::string algorithm) {
    std::unordered_set<int> affected_vertices;
    std::vector<std::vector<std::vector<int>>> inc_index;
    inc_index.resize(G->n);
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
    if (algorithm.rfind("BFS", 0) == 0) {
        Q_BFS.push(start);
    }
    else {
        Q_PR.push(start);
    }

    while (1) {
        std::vector<int> current;
        if (algorithm.rfind("BFS", 0) == 0) {
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
        if (algorithm.rfind("BFS", 0) == 0) {
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
            affected_vertices.insert(v);
            while (next_idx[v].size() < L[v].size()) {
                next_idx[v].push_back(L[v].size());
            }
            std::vector<int> next;
            next.push_back(u);
            next.push_back(ts);
            next.push_back(te);
            next.push_back(d);
            inc_index[v].push_back(next);
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
            if (t_threshold != -1 && te_new - ts_new + 1 > t_threshold) {
                e = e->next;
                continue;
            }
            bool flag = false;
            // enqueue checking
            if (algorithm.rfind("BFS", 0) == 0) {
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
                if (algorithm.rfind("BFS", 0) == 0) {
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
                if (algorithm.rfind("BFS", 0) == 0) {
                    Q_BFS.push(next);
                }
                else {
                    Q_PR.push(next);
                }
            }
            e = e->next;
        }
    }

    if (algorithm == "PrioritySearch" || algorithm == "BFS-naive") {
        for (auto v : affected_vertices) {
            for (auto path : inc_index[v]) {
                L[v].push_back(path);
            }
        }
    }
    else {
        // Compress the index to form a minimal set of intervals
        for (auto v : affected_vertices) {
            int i = 0;
            int j = 0;
            std::vector<int> BIT;
            BIT.assign(G->tmax + 2, G->tmax + 1);
            for (i = 0; i <= inc_index[v].size(); i++) {
                if (i == inc_index[v].size() || (i > 0 && inc_index[v][i][3] != inc_index[v][i - 1][3])) {
                    std::sort(inc_index[v].begin() + j, inc_index[v].begin() + i, cmp1);
                    while (j < i) {
                        int t = inc_index[v][j][1] + 1;
                        bool flag = false;
                        while (t <= G->tmax + 1) {
                            if (BIT[t] <= inc_index[v][j][2]) {
                                flag = true;
                                break;
                            }
                            t += (t & (-t));
                        }
                        if (!flag) {
                            L[v].push_back(inc_index[v][j]);
                            int t = inc_index[v][j][1] + 1;
                            while (t > 0) {
                                BIT[t] = std::min(BIT[t], inc_index[v][j][2]);
                                t -= (t & (-t));
                            }
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
}

AdvancedTwoHopIndex::AdvancedTwoHopIndex(TemporalGraph* G, int k_input, int t_threshold, std::string algorithm) {
    k = k_input;
    index_construct_algorithm = algorithm;
    L_in.resize(G->n);
    L_out.resize(G->n);
    next_in.resize(G->n);
    next_out.resize(G->n);
    cut_in.resize(G->n);
    cut_out.resize(G->n);

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
        construct_for_a_vertex(G, u, false, L_in, next_in, cut_in, t_threshold, algorithm);
        construct_for_a_vertex(G, u, true, L_out, next_out, cut_out, t_threshold, algorithm);
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