#include "baseline.h"

bool cmp(std::pair<int, int> i, std::pair<int, int> j) {
    return i.first > j.first || (i.first == j.first && i.second < j.second);
}

bool cmp1(std::pair<std::pair<int, int>, int> i, std::pair<std::pair<int, int>, int> j) {
    return i.second < j.second;
}

unsigned long long BaselineIndex::size() {
    unsigned long long num_intervals = 0;
    for (auto it = L.begin(); it != L.end(); it++) {
        for (auto it1 = it->begin(); it1 != it->end(); it1++) {
            num_intervals += it1->second.size();
        }
    }
    return num_intervals;
}

bool BaselineIndex::reachable(TemporalGraph* G, int u, int v, int ts, int te, int k_input) {
    if (k_input == 1) {
        if (u == v) {
            return true;
        }
        TemporalGraph::Edge* e = G->getHeadEdge(u);
        while (e) {
            if (e->interaction_time >= ts && e->interaction_time <= te && e->to == v) {
                return true;
            }
            e = e->next;
        }
        return false;
    }
    if (u == v) {
        return true;
    }
    if (vertex_cover.find(u) == vertex_cover.end()) {
        if (vertex_cover.find(v) != vertex_cover.end()) {
            TemporalGraph::Edge* e = G->getHeadEdge(u);
            while (e) {
                if (e->interaction_time >= ts && e->interaction_time <= te) {
                    if ((!is_temporal_path && reachable(G, e->to, v, ts, te, k_input - 1)) || 
                        (is_temporal_path && reachable(G, e->to, v, e->interaction_time, te, k_input - 1))) {
                        return true;
                    }
                }
                e = e->next;
            }
            return false;
        }
        else {
            TemporalGraph::Edge* e1 = G->getHeadEdge(u);
            while (e1) {
                if (e1->interaction_time >= ts && e1->interaction_time <= te) {
                    TemporalGraph::Edge* e2 = G->getHeadInEdge(v);
                    while (e2) {
                        if (e2->interaction_time >= ts && e2->interaction_time <= te) {
                            if ((!is_temporal_path && reachable(G, e1->to, e2->to, ts, te, k_input - 2)) || 
                                (is_temporal_path && e2->interaction_time >= e1->interaction_time && reachable(G, e1->to, e2->to, e1->interaction_time, e2->interaction_time, k_input - 2))) {
                                return true;
                            }
                        }
                        e2 = e2->next;
                    }
                }
                e1 = e1->next;
            }
            return false;
        }
    }
    else {
        if (vertex_cover.find(v) != vertex_cover.end()) {
            int i = inv_vertex_cover[u];
            if (L[i].find(v) == L[i].end()) {
                return false;
            }
            for (int j = 1; j <= k_input; j++) {
                int l = 0;
                if (j > 1) {
                    l = cut[i][v][j - 2];
                }
                int r = cut[i][v][j - 1] - 1;
                if (l > r) {
                    continue;
                }
                while (l < r) {
                    int mid = l + r + 1 >> 1;
                    if (L[i][v][mid].first >= ts) {
                        l = mid;
                    }
                    else {
                        r = mid - 1;
                    }
                }
                if (L[i][v][l].first >= ts && L[i][v][l].second <= te) {
                    return true;
                }
            }
            return false;
        }
        else {
            TemporalGraph::Edge* e = G->getHeadInEdge(v);
            while (e) {
                if (e->interaction_time >= ts && e->interaction_time <= te) {
                    if ((!is_temporal_path && reachable(G, u, e->to, ts, te, k_input - 1)) ||
                        (is_temporal_path && reachable(G, u, e->to, ts, e->interaction_time, k_input - 1))) {
                        return true;
                    }
                }
                e = e->next;
            }
            return false;
        }
    }
}

BaselineIndex::BaselineIndex(TemporalGraph* G, int k_input, int t_threshold, std::string path_type) {
    k = k_input;
    if (path_type == "Temporal") {
        is_temporal_path = true;
    }
    else {
        is_temporal_path = false;
    }

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
    cut.resize(vertex_cover.size());
    L.resize(vertex_cover.size());
    std::vector<std::vector<std::pair<std::pair<int, int>, int>>> T;
    T.resize(G->n);
    std::unordered_set<int> Vn, Vs;
    std::vector<std::pair<int, int>> intervals;
    
    // Construct the index by vertex cover
    int i = 0;
    unsigned long long start_time = currentTime();
    for (auto u : vertex_cover) {
        inv_vertex_cover[u] = i++;
        std::queue<std::vector<int>> Q;
        for (auto v : Vn) {
            T[v].clear();
        }
        T.clear();
        Vn.clear();
        Vs.clear();
        Q.push(std::vector<int>{u, G->tmax + 1, -1, 0});

        // BFS to find minimal paths
        while (!Q.empty()) {
            std::vector<int> current = Q.front();
            Q.pop();
            int v = current[0], ts = current[1], te = current[2], d = current[3];
            // Check minimality to avoid expanding non-minimal paths
            bool flag = false;
            for (auto it1 = T[v].begin(); it1 != T[v].end(); it1++) {
                if (it1->first.first >= ts && it1->first.second <= te && it1->second == d && (it1->first.first != ts || it1->first.second != te)) {
                    flag = true;
                    break;
                }
            }
            if (flag) {
                continue;
            }
            if ((u == v && d != 0)) {
                continue;
            }
            TemporalGraph::Edge* e = G->getHeadEdge(v);
            while (e) {
                int w = e->to, t = e->interaction_time;
                // Note that for temporal paths, edges should follow an increasing time order
                if (!is_temporal_path || te <= t) {
                    int ts_new = std::min(ts, t);
                    int te_new = std::max(te, t);
                    if (t_threshold == -1 || te_new - ts_new + 1 <= t_threshold) {
                        visited_paths++;
                        bool flag = false;
                        for (auto it1 = T[w].begin(); it1 != T[w].end();) {
                            if (it1->first.first >= ts_new && it1->first.second <= te_new) {
                                flag = true;
                                break;
                            }
                            if (it1->first.first <= ts_new && it1->first.second >= te_new && it1->second >= d + 1) {
                                it1 = T[w].erase(it1);
                                continue;
                            }
                            it1++;
                        }
                        if (!flag) {
                            if (w != u && Vs.find(w) == Vs.end() && vertex_cover.find(w) != vertex_cover.end()) {
                                Vs.insert(w);
                            }
                            if (Vn.find(w) == Vn.end()) {
                                Vn.insert(w);
                            }
                            T[w].push_back(std::make_pair(std::make_pair(ts_new, te_new), d + 1));
                            if (d + 1 < k) {
                                Q.push(std::vector<int>{w, ts_new, te_new, d + 1});
                            }
                        }
                    }
                }
                e = e->next;
            }
        }

        for (auto v : Vn) {
            max_number_of_paths = std::max(max_number_of_paths, (unsigned long long)T[v].size());
        }

        // Compress all enqueued paths to generate a minimal path set
        for (auto v : Vs) {
            L[i - 1][v] = std::vector<std::pair<int, int>>();
            cut[i - 1][v] = std::vector<int>();
            for (int j = 1; j <= k; j++) {
                intervals.clear();
                for (auto it1 = T[v].begin(); it1 != T[v].end(); it1++) {
                    if (it1->second > j) {
                        break;
                    }
                    intervals.push_back(it1->first);
                }
                std::sort(intervals.begin(), intervals.end(), cmp);
                int tmin = G->tmax + 1;
                for (auto it2 = intervals.begin(); it2 != intervals.end(); it2++) {
                    if (tmin > it2->second) {
                        tmin = it2->second;
                        L[i - 1][v].push_back(*it2);
                    }
                }
                cut[i - 1][v].push_back(L[i - 1][v].size());
            }
        }
        
        putProcess(double(i) / vertex_cover.size(), currentTime() - start_time);
    }
}

void BaselineIndex::solve(TemporalGraph* G, char* query_file, char* output_file) {
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