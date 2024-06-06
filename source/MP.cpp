#include "MP.h"

bool cmp(std::pair<int, int> i, std::pair<int, int> j) {
    return i.first > j.first || (i.first == j.first && i.second < j.second);
}

bool cmp1(std::pair<std::pair<int, int>, int> i, std::pair<std::pair<int, int>, int> j) {
    return i.second < j.second;
}

unsigned long long MPIndex::size() {
    unsigned long long num_intervals = 0;
    for (auto it = L.begin(); it != L.end(); it++) {
        for (auto it1 = it->begin(); it1 != it->end(); it1++) {
            for (int i = 0; i <= k; i++) {
                num_intervals += it1->second[i].size();
            }
        }
    }
    return num_intervals;
}

bool MPIndex::reachable(TemporalGraph* G, int u, int v, int ts, int te, int k_input) {
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
                for (auto e : L[i][v][j]) {
                    if (e.first >= ts && e.second <= te) {
                        return true;
                    }
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

MPIndex::MPIndex(TemporalGraph* G, int k_input) {
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
    
    // Index construction
    int i = 0;
    unsigned long long start_time = currentTime();
    for (auto u : vertex_cover) {
        inv_vertex_cover[u] = i;
        std::queue<std::vector<int>> Q;
        Q.push(std::vector<int>{u, G->tmax + 1, -1, 0});
        while (!Q.empty()) {
            std::vector<int> current = Q.front();
            Q.pop();
            int v = current[0], ts = current[1], te = current[2], d = current[3];
            if (d == k) {
                break;
            }
            bool is_minimal = false;
            for (int j = 1; j <= d; j++) {
                for (auto e : L[i][v][j]) {
                    if (e.first == ts && e.second == te) {
                        is_minimal = true;
                        break;
                    }
                }
                if (is_minimal) {
                    break;
                }
            }
            if (d > 0 && !is_minimal) {
                // Only minimal k-hop paths are expanded
                continue;
            }
            TemporalGraph::Edge* e = G->getHeadEdge(v);
            while (e) {
                int w = e->to, t = e->interaction_time;
                if (w != u && te <= t) {
                    int ts_new = std::min(ts, t), te_new = std::max(te, t);
                    if (L[i].find(w) == L[i].end()) {
                        L[i][w] = std::vector<std::vector<std::pair<int, int>>>();
                        L[i][w].resize(k + 1);
                    }
                    bool is_minimal = true;
                    for (int j = 1; j <= d + 1; j++) {
                        for (auto e : L[i][w][j]) {
                            if (e.first >= ts_new && e.second <= te_new) {
                                is_minimal = false;
                                break;
                            }
                        }
                        if (!is_minimal) {
                            break;
                        }
                    }
                    if (is_minimal) {
                        for (auto it = L[i][w][d + 1].begin(); it != L[i][w][d + 1].end();) {
                            if (it->first <= ts_new && it->second >= te_new) {
                                it = L[i][w][d + 1].erase(it);
                                continue;
                            }
                            it++;
                        }
                        L[i][w][d + 1].push_back(std::make_pair(ts_new, te_new));
                        Q.push(std::vector<int>{w, ts_new, te_new, d + 1});
                    }
                }
                e = G->getNextEdge(e);
            }
        }

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

void MPIndex::solve(TemporalGraph* G, char* query_file, char* output_file) {
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