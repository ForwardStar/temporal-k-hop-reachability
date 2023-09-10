#include "index.h"

bool cmp(std::pair<int, int> i, std::pair<int, int> j) {
    return i.second > j.second;
}

bool cmp1(std::pair<int, int> i, std::pair<int, int> j) {
    return i.first > j.first;
}

bool Index::reachable(TemporalGraph* G, int u, int v, int ts, int te, int k_input) {
    if (u == v) {
        return true;
    }
    if (vertex_cover.find(u) == vertex_cover.end()) {
        if (vertex_cover.find(v) != vertex_cover.end()) {
            TemporalGraph::Edge* e = G->getHeadEdge(u);
            while (e) {
                if (e->interaction_time >= ts && e->interaction_time <= te) {
                    if (reachable(G, e->to, v, ts, te, k_input - 1)) {
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
                            if (reachable(G, e1->to, e2->to, ts, te, k_input - 2)) {
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
            if (index[i].find(v) == index[i].end()) {
                return false;
            }
            for (int j = 0; j <= k_input - (k - 2); j++) {
                for (auto it = index[i][v][j].begin(); it != index[i][v][j].end(); it++) {
                    int ts_index = -*it / (G->tmax + 1);
                    int te_index = -*it % (G->tmax + 1);
                    if (ts_index >= ts) {
                        if (te_index <= te) {
                            return true;
                        }
                    }
                    else {
                        break;
                    }
                }
            }
            return false;
        }
        else {
            TemporalGraph::Edge* e = G->getHeadInEdge(v);
            while (e) {
                if (e->interaction_time >= ts && e->interaction_time <= te) {
                    if (reachable(G, u, e->to, ts, te, k_input - 1)) {
                        return true;
                    }
                }
                e = e->next;
            }
            return false;
        }
    }
}

Index::Index(TemporalGraph* G, int k_input) {
    k = k_input;

    // Generate vertex cover
    TemporalGraph* Graph = new TemporalGraph(G, 0, G->tmax);
    Heap* heap = new Heap();
    std::vector<bool> covered;
    covered.resize(Graph->n);

    for (int i = 0; i < Graph->edge_set.size(); i++) {
        int u = Graph->edge_set[i].first.first;
        int v = Graph->edge_set[i].first.second;
        heap->insert(i, -(Graph->degree[u] + 1) * (Graph->in_degree[u] + 1) - (Graph->degree[v] + 1) * (Graph->in_degree[v] + 1));
    }
    while (heap->size() > 0) {
        int i = heap->pop();
        int u = Graph->edge_set[i].first.first;
        int v = Graph->edge_set[i].first.second;
        if (covered[u] || covered[v]) {
            continue;
        }
        covered[u] = true;
        covered[v] = true;
        vertex_cover.insert(u);
        vertex_cover.insert(v);
    }
    std::cout << "Vertex cover size: " << vertex_cover.size() << std::endl;
    
    // Construct the index by vertex cover
    int i = 0;
    unsigned long long start_time = currentTime();
    for (auto it = vertex_cover.begin(); it != vertex_cover.end(); it++) {
        int u = *it;
        inv_vertex_cover[u] = i++;
        index.push_back(std::unordered_map<int, std::vector<std::set<long long>>>());
        auto cmp2 = [](std::vector<int> i, std::vector<int> j) {
            return i[2] - i[1] > j[2] - j[1];
        };
        std::priority_queue<std::vector<int>, std::vector<std::vector<int>>, decltype(cmp2)> Q(cmp2);
        std::vector<int> start;
        std::vector<std::vector<std::vector<std::pair<int, int>>>> visited;
        visited.resize(Graph->n);
        for (int v = 0; v < Graph->n; v++) {
            visited[v].resize(k + 1);
        }
        start.push_back(u);
        start.push_back(Graph->tmax + 1);
        start.push_back(-1);
        start.push_back(0);
        Q.push(start);
        while (!Q.empty()) {
            std::vector<int> current = Q.top();
            Q.pop();
            int v = current[0];
            if ((u == v && current[3] != 0) || current[3] == k) {
                continue;
            }
            TemporalGraph::Edge* e = Graph->getHeadEdge(v);
            while (e) {
                int ts = std::min(current[1], e->interaction_time);
                int te = std::max(current[2], e->interaction_time);
                bool flag = false;
                for (int j = 0; j <= current[3] + 1; j++) {
                    for (auto it1 = visited[e->to][j].begin(); it1 != visited[e->to][j].end(); it1++) {
                        if (it1->first >= ts && it1->second <= te) {
                            flag = true;
                            break;
                        }
                    }
                    if (flag) {
                        break;
                    }
                }
                if (!flag) {
                    std::vector<int> next;
                    next.push_back(e->to);
                    next.push_back(ts);
                    next.push_back(te);
                    next.push_back(current[3] + 1);
                    visited[e->to][current[3] + 1].push_back(std::pair<int, int>(ts, te));
                    Q.push(next);
                }
                e = e->next;
            }
        }
        for (auto it1 = vertex_cover.begin(); it1 != vertex_cover.end(); it1++) {
            int v = *it1;
            if (u == v) {
                continue;
            }
            for (int j = 0; j <= k; j++) {
                for (auto it2 = visited[v][j].begin(); it2 != visited[v][j].end(); it2++) {
                    if (index[i - 1].find(v) == index[i - 1].end()) {
                        index[i - 1][v] = std::vector<std::set<long long>>();
                        for (int l = 0; l < 3; l++) {
                            index[i - 1][v].push_back(std::set<long long>());
                        }
                    }
                    index[i - 1][v][std::max(0, j - (k - 2))].insert(-(long long)it2->first * (Graph->tmax + 1) - it2->second);
                }
            }
        }
        putProcess(double(i) / vertex_cover.size(), currentTime() - start_time);
    }

    delete Graph;
}

void Index::solve(TemporalGraph* Graph, char* query_file, char* output_file, int k) {
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
        if (reachable(Graph, s, t, ts, te, k)) {
            fout << "Reachable" << std::endl;
        }
        else {
            fout << "Not reachable" << std::endl;
        }
        putProcess(double(++i) / query_num, currentTime() - start_time);
    }

    std::cout << "Average: " << timeFormatting((currentTime() - start_time) / query_num).str() << std::endl;
}