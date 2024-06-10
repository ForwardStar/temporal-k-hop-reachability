#include "MP.h"

bool MP_cmp(std::pair<int, int> i, std::pair<int, int> j) {
    return i.first < j.first;
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

    if (L[u].find(v) == L[u].end()) {
        return false;
    }
    for (int j = 1; j <= k_input; j++) {
        int l = 0, r = L[u][v][j].size() - 1;
        if (r == -1) {
            continue;
        }
        while (l < r) {
            int mid = (l + r) / 2;
            if (L[u][v][j][mid].first < ts) {
                l = mid + 1;
            }
            else {
                r = mid;
            }
        }
        if (L[u][v][j][l].first >= ts && L[u][v][j][l].second <= te) {
            return true;
        }
    }
    
    return false;
}

MPIndex::MPIndex(TemporalGraph* G, int k_input) {
    k = k_input;
    L.resize(G->n);
    
    // Index construction
    unsigned long long start_time = currentTime();
    for (int u = 0; u < G->n; u++) {
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
                for (auto e : L[u][v][j]) {
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
                    if (L[u].find(w) == L[u].end()) {
                        L[u][w] = std::vector<std::vector<std::pair<int, int>>>();
                        L[u][w].resize(k + 1);
                    }
                    bool is_minimal = true;
                    for (int j = 1; j <= d + 1; j++) {
                        for (auto e : L[u][w][j]) {
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
                        for (auto it = L[u][w][d + 1].begin(); it != L[u][w][d + 1].end();) {
                            if (it->first <= ts_new && it->second >= te_new) {
                                it = L[u][w][d + 1].erase(it);
                                continue;
                            }
                            it++;
                        }
                        L[u][w][d + 1].push_back(std::make_pair(ts_new, te_new));
                        Q.push(std::vector<int>{w, ts_new, te_new, d + 1});
                    }
                }
                e = G->getNextEdge(e);
            }
        }

        for (auto it = L[u].begin(); it != L[u].end();) {
            unsigned long long num_paths = 0;
            for (int j = 1; j <= k; j++) {
                num_paths += (unsigned long long)it->second[j].size();
                std::sort(it->second[j].begin(), it->second[j].end(), MP_cmp);
            }
            alpha = std::max(num_paths, alpha);
            it++;
        }
        
        putProcess(double(u + 1) / G->n, currentTime() - start_time);
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