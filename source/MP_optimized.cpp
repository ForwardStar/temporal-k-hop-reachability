#include "MP_optimized.h"

bool cmp_t1(std::pair<std::pair<int, int>, int> i, std::pair<std::pair<int, int>, int> j) {
    return i.second < j.second;
}

bool cmp1(std::pair<int, int> i, std::pair<int, int> j) {
    return i.first < j.first;
}

unsigned long long MPIndexO::size() {
    unsigned long long num_intervals = 0;
    for (int i = 0; i < L.size(); i++) {
        for (auto it = L[i].begin(); it != L[i].end(); it++) {
            for (int j = 0; j <= k; j++) {
                num_intervals += it->second[j].size();
            }
        }
    }
    return num_intervals;
}

bool MPIndexO::reachable(TemporalGraph* G, int u, int v, int ts, int te, int k_input) {
    if (u == v) {
        return true;
    }

    if (vertex_cover.find(u) != vertex_cover.end()) {
        if (vertex_cover.find(v) != vertex_cover.end()) {
            int i = inv_vertex_cover[u];
            if (L[i].find(v) == L[i].end()) {
                return false;
            }
            for (int j = 1; j <= k_input; j++) {
                int l = 0, r = L[i][v][j].size() - 1;
                if (r == -1) {
                    continue;
                }
                while (l < r) {
                    int mid = (l + r) / 2;
                    if (L[i][v][j][mid].first < ts) {
                        l = mid + 1;
                    }
                    else {
                        r = mid;
                    }
                }
                if (L[i][v][j][l].first >= ts && L[i][v][j][l].second <= te) {
                    return true;
                }
            }
        }
        else {
            auto e = G->getHeadInEdge(v);
            while (e) {
                if (e->interaction_time >= ts && e->interaction_time <= te && reachable(G, u, e->to, ts, e->interaction_time, k_input - 1)) {
                    return true;
                }
                e = G->getNextEdge(e);
            }
        }
    }
    else {
        auto e = G->getHeadEdge(u);
        while (e) {
            if (e->interaction_time >= ts && e->interaction_time <= te && reachable(G, e->to, v, e->interaction_time, te, k_input - 1)) {
                return true;
            }
            e = G->getNextEdge(e);
        }
    }
    
    return false;
}

MPIndexO::MPIndexO(TemporalGraph* G, int k_input) {
    k = k_input;

    // Generate vertex cover
    auto large_degree_first = [](std::pair<int, long> i, std::pair<int, long> j) {
        return i.second < j.second;
    };
    std::priority_queue<std::pair<int, long>, std::vector<std::pair<int, long>>, decltype(large_degree_first)> heap(large_degree_first);
    std::vector<bool> covered;
    covered.resize(G->n);

    for (int i = 0; i < G->edge_set.size(); i++) {
        int u = G->edge_set[i].first.first;
        int v = G->edge_set[i].first.second;
        heap.push(std::make_pair(i, ((long long)G->degree[u] + 1) * (G->in_degree[u] + 1) + ((long long)G->degree[v] + 1) * (G->in_degree[v] + 1)));
    }
    while (heap.size() > 0) {
        auto e = heap.top();
        heap.pop();
        int u = G->edge_set[e.first].first.first;
        int v = G->edge_set[e.first].first.second;
        if (covered[u] || covered[v]) {
            continue;
        }
        covered[u] = true;
        covered[v] = true;
        vertex_cover.insert(u);
        vertex_cover.insert(v);
    }
    std::cout << "Vertex cover size: " << vertex_cover.size() << std::endl;
    L.resize(vertex_cover.size());
    
    // Construct the index by vertex cover
    int i = 0;
    f.assign(G->n, G->n);
    unsigned long long start_time = currentTime();
    for (auto u : vertex_cover) {
        inv_vertex_cover[u] = i;

        // Find the edges in the k-hop subgraph of u
        for (auto v : visited_vertices) {
            f[v] = G->n;
        }
        visited_vertices.clear();
        f[u] = 0;
        visited_vertices.insert(u);
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
            while (e) {
                edges.push_back(std::make_pair(std::make_pair(v, e->to), e->interaction_time));
                if (f[v] + 1 < f[e->to]) {
                    visited_vertices.insert(e->to);
                    Q.push(e->to);
                    f[e->to] = f[v] + 1;
                }
                e = G->getNextEdge(e);
            }
        }
        std::sort(edges.begin(), edges.end(), cmp_t1);

        // Index construction
        for (auto e : edges) {
            int v = e.first.first, w = e.first.second, t = e.second;
            if (u == w) {
                continue;
            }
            else if (u == v) {
                if (L[i].find(w) == L[i].end()) {
                    L[i][w] = std::vector<std::vector<std::pair<int, int>>>();
                    L[i][w].resize(k + 1);
                }
                bool is_minimal = true;
                for (auto path : L[i][w][1]) {
                    if (path.first == t && path.second == t) {
                        is_minimal = false;
                        break;
                    }
                }
                if (is_minimal) {
                    L[i][w][1].push_back(std::make_pair(t, t));
                }
            }
            else {
                if (L[i].find(v) != L[i].end()) {
                    for (int j = 1; j < k; j++) {
                        int max_ts = -1;
                        for (auto e : L[i][v][j]) {
                            if (e.second <= t) {
                                max_ts = std::max(max_ts, e.first);
                            }
                        }
                        if (max_ts != -1) {
                            if (L[i].find(w) == L[i].end()) {
                                L[i][w] = std::vector<std::vector<std::pair<int, int>>>();
                                L[i][w].resize(k + 1);
                            }
                            bool is_minimal = true;
                            for (int j2 = 1; j2 <= j + 1; j2++) {
                                for (auto e : L[i][w][j2]) {
                                    if (e.first >= max_ts && e.second <= t) {
                                        is_minimal = false;
                                        break;
                                    }
                                }
                                if (!is_minimal) {
                                    break;
                                }
                            }
                            if (is_minimal) {
                                for (int j2 = j + 1; j2 <= k; j2++) {
                                    for (auto it = L[i][w][j2].begin(); it != L[i][w][j2].end();) {
                                        if (it->first <= max_ts && it->second >= t) {
                                            it = L[i][w][j2].erase(it);
                                            continue;
                                        }
                                        it++;
                                    }
                                }
                                L[i][w][j + 1].push_back(std::make_pair(max_ts, t));
                            }
                        }
                    }
                }
            }
        }

        // Prune index to keep only vertices in the vertex cover
        for (auto it = L[i].begin(); it != L[i].end();) {
            if (vertex_cover.find(it->first) == vertex_cover.end()) {
                it = L[i].erase(it);
                continue;
            }
            unsigned long long num_paths = 0;
            for (int j = 1; j <= k; j++) {
                num_paths += (unsigned long long)it->second[j].size();
            }
            alpha = std::max(num_paths, alpha);
            it++;
        }
        
        putProcess(double(++i) / vertex_cover.size(), currentTime() - start_time);
    }
}

void MPIndexO::solve(TemporalGraph* G, char* query_file, char* output_file) {
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