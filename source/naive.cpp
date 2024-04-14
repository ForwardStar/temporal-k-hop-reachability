#include "naive.h"

unsigned long long NaiveIndex::size() {
    unsigned long long num_vertices = 0;
    for (int ts = 0; ts < L.size(); ts++) {
        for (int te = ts; te < L[ts].size(); te++) {
            for (int u = 0; u < L[ts][te].size(); u++) {
                num_vertices += L[ts][te][u].size();
            }
        }
    }
    return num_vertices;
}

bool NaiveIndex::reachable(TemporalGraph* G, int u, int v, int ts, int te, int k_input) {
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
    if (k_input == 0) {
        return false;
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
            if (L[ts][te][i].find(v) != L[ts][te][i].end() && L[ts][te][i][v] <= k) {
                return true;
            }
            else {
                return false;
            }
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

NaiveIndex::NaiveIndex(TemporalGraph* G, int k_input, int t_threshold, std::string path_type) {
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
    // std::cout << "Vertex cover size: " << vertex_cover.size() << std::endl;
    int i = 0;
    for (auto u : vertex_cover) {
        inv_vertex_cover[u] = i++;
    }

    unsigned long long start_time = currentTime();
    L.resize(G->tmax + 1);
    for (int ts = 0; ts <= G->tmax; ts++) {
        L[ts].resize(G->tmax + 1);
        for (int te = ts; te <= G->tmax; te++) {
            TemporalGraph* Gp = G->projectedGraph(ts, te);
            L[ts][te].resize(vertex_cover.size());
            // Construct the index by vertex cover
            int i = 0;
            for (auto u : vertex_cover) {
                L[ts][te][i] = std::unordered_map<int, int>();
                std::vector<int> T;
                T.assign(Gp->n, 2147483647);
                T[u] = -1;
                std::queue<std::vector<int>> Q;
                Q.push(std::vector<int>{u, 0, -1});
                while (!Q.empty()) {
                    int v = Q.front()[0];
                    int dis = Q.front()[1];
                    int t_end = Q.front()[2];
                    if (dis >= k) {
                        break;
                    }
                    if (t_end > T[v]) {
                        continue;
                    }
                    Q.pop();
                    TemporalGraph::Edge* edge = Gp->getHeadEdge(v);
                    while (edge) {
                        if (path_type == "Temporal" && (edge->interaction_time >= T[edge->to] || t_end > edge->interaction_time)) {
                            edge = G->getNextEdge(edge);
                            continue;
                        }
                        if (path_type == "Temporal" || T[edge->to] == 2147483647) {
                            if (vertex_cover.find(edge->to) != vertex_cover.end() && L[ts][te][i].find(edge->to) == L[ts][te][i].end()) {
                                L[ts][te][i][edge->to] = dis + 1;
                            }
                            T[edge->to] = std::min(T[edge->to], edge->interaction_time);
                            Q.push(std::vector<int>{edge->to, dis + 1, edge->interaction_time});
                            edge = Gp->getNextEdge(edge);
                        }
                    }
                }
                i++;
            }
            delete Gp;
        }
        putProcess(double(ts + 1) / (G->tmax + 1), currentTime() - start_time);
    }
}

void NaiveIndex::solve(TemporalGraph* G, char* query_file, char* output_file) {
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