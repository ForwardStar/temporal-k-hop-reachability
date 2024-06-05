#include "T2H.h"

bool cmp_t2(std::pair<std::pair<int, int>, int> i, std::pair<std::pair<int, int>, int> j) {
    return i.second < j.second;
}

bool cmp_degree(std::pair<int, long long> i, std::pair<int, long long> j) {
    return i.second > j.second;
}

unsigned long long T2HIndex::size() {
    unsigned long long num_intervals = 0;
    for (int u = 0; u < L_out.size(); u++) {
        for (auto s1 : L_out[u]) {
            for (int i = 0; i <= k; i++) {
                num_intervals += s1.second[i].size();
            }
        }
    }
    for (int u = 0; u < L_in.size(); u++) {
        for (auto s1 : L_in[u]) {
            for (int i = 0; i <= k; i++) {
                num_intervals += s1.second[i].size();
            }
        }
    }
    return num_intervals;
}

unsigned long long T2HIndex::max_number_of_paths() {
    unsigned long long res = 0;
    for (int u = 0; u < L_out.size(); u++) {
        for (auto s1 : L_out[u]) {
            unsigned long long tmp = 0;
            for (int i = 0; i <= k; i++) {
                tmp += s1.second[i].size();
            }
            res = std::max(res, tmp);
        }
    }
    for (int u = 0; u < L_in.size(); u++) {
        for (auto s1 : L_in[u]) {
            unsigned long long tmp = 0;
            for (int i = 0; i <= k; i++) {
                tmp += s1.second[i].size();
            }
            res = std::max(res, tmp);
        }
    }
    return res;
}

std::pair<int, int>& T2HIndex::binary_search1(std::vector<std::unordered_map<int, std::vector<std::vector<std::pair<int, int>>>>> &L, int u, int v, int k, int ts) {
    int l = 0;
    int r = L[u][v][k].size() - 1;
    if (r == -1) {
        return null_interval;
    }
    while (l < r) {
        int mid = (l + r) / 2;
        if (L[u][v][k][mid].first < ts) {
            l = mid + 1;
        }
        else {
            r = mid;
        }
    }
    return L[u][v][k][l];
}

std::pair<int, int>& T2HIndex::binary_search2(std::vector<std::unordered_map<int, std::vector<std::vector<std::pair<int, int>>>>> &L, int u, int v, int k, int te) {
    int l = 0;
    int r = L[u][v][k].size() - 1;
    if (r == -1) {
        return null_interval;
    }
    while (l < r) {
        int mid = (l + r + 1) / 2;
        if (L[u][v][k][mid].first > te) {
            r = mid - 1;
        }
        else {
            l = mid;
        }
    }
    return L[u][v][k][l];
}

bool T2HIndex::reachable(int u, int v, int ts, int te, int k) {
    if (u == v) {
        return true;
    }

    if (L_out[u].find(v) != L_out[u].end()) {
        for (int i = 0; i <= k; i++) {
            auto interval = binary_search1(L_out, u, v, i, ts);
            if (interval.first >= ts && interval.second <= te) {
                return true;
            }
        }
    }

    if (L_in[v].find(u) != L_in[v].end()) {
        for (int i = 0; i <= k; i++) {
            auto interval = binary_search1(L_in, v, u, i, ts);
            if (interval.first >= ts && interval.second <= te) {
                return true;
            }
        }
    }

    for (auto p : L_out[u]) {
        int w = p.first;
        if (L_in[v].find(w) != L_in[v].end()) {
            int d1 = k - 1, d2 = 1, ts_max = -1;
            while (d1 >= 0) {
                auto interval1 = binary_search2(L_in, v, w, d2, te);
                if (interval1.first >= ts && interval1.second <= te) {
                    ts_max = std::max(ts_max, interval1.first);
                }

                auto interval2 = binary_search1(L_out, u, w, d1, ts);
                if (interval2.first >= ts && interval2.second <= ts_max) {
                    return true;
                }
                d1--, d2++;
            }
        }
    }

    return false;
}

void T2HIndex::construct_for_a_vertex(TemporalGraph* G, int u, bool reverse, std::vector<std::unordered_map<int, std::vector<std::vector<std::pair<int, int>>>>> &L) {
    std::vector<std::vector<int>> current;
    std::vector<std::vector<int>> next;
    current.push_back(std::vector<int>{u, G->tmax + 1, -1, 0});

    // Find the edges in the k-hop subgraph of u
    std::vector<int> f;
    f.assign(G->n, G->n);
    f[u] = 0;
    std::vector<std::pair<std::pair<int, int>, int>> edges;
    std::queue<int> Q;
    Q.push(u);
    while (!Q.empty()) {
        int v = Q.front();
        Q.pop();
        if (f[v] >= k) {
            break;
        }
        auto e = G->getHeadEdge(v);
        if (reverse) {
            e = G->getHeadInEdge(v);
        }
        while (e) {
            if (order[e->to] > order[u]) {
                edges.push_back(std::make_pair(std::make_pair(v, e->to), e->interaction_time));
                if (f[v] + 1 < f[e->to]) {
                    Q.push(e->to);
                    f[e->to] = f[v] + 1;
                }
            }
            e = G->getNextEdge(e);
        }
    }
    std::sort(edges.begin(), edges.end(), cmp_t2);

    // Index construction
    for (auto e : edges) {
        int v = e.first.first, w = e.first.second, t = e.second;
        if (u == w) {
            continue;
        }
        else if (u == v) {
            if (L[w].find(u) == L[w].end()) {
                L[w][u] = std::vector<std::vector<std::pair<int, int>>>();
                L[w][u].resize(k + 1);
            }
            if ((!reverse && !reachable(u, w, t, t, 1)) || (reverse && !reachable(w, u, t, t, 1))) {
                L[w][u][1].push_back(std::make_pair(t, t));
            }
        }
        else {
            if (L[v].find(u) != L[v].end()) {
                for (int j = 1; j < k; j++) {
                    auto interval = binary_search2(L, v, u, j, t);
                    if (interval.first >= 0 && interval.second <= t) {
                        int ts = interval.first, te = t;
                        if ((!reverse && !reachable(u, w, ts, te, j + 1)) || (reverse && !reachable(w, u, ts, te, j + 1))) {
                            if (L[w].find(u) == L[w].end()) {
                                L[w][u] = std::vector<std::vector<std::pair<int, int>>>();
                                L[w][u].resize(k + 1);
                            }
                            auto interval1 = binary_search2(L, w, u, j + 1, te);
                            if (interval1.second == te) {
                                interval1.first = ts;
                            }
                            else {
                                L[w][u][j + 1].push_back(std::make_pair(ts, te));
                            }
                        }
                    }
                }
            }
        }
    }
}

T2HIndex::T2HIndex(TemporalGraph* G, int k_input) {
    k = k_input;
    L_in.resize(G->n);
    L_out.resize(G->n);

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
        construct_for_a_vertex(G, u, false, L_in);
        construct_for_a_vertex(G, u, true, L_out);
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