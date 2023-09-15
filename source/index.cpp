#include "index.h"

bool cmp(std::pair<int, int> i, std::pair<int, int> j) {
    return i.first > j.first || (i.first == j.first && i.second < j.second);
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
            int l = 0;
            int r = cut[i][v][k_input - (k - 2)] - 1;
            if (l > r) {
                return false;
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
            if (L[i][v][l].second <= te) {
                return true;
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

Index::Index(TemporalGraph* G, int k_input, int t_threshold) {
    k = k_input;

    // Generate vertex cover
    Heap* heap = new Heap();
    std::vector<bool> covered;
    covered.resize(G->n);

    for (int i = 0; i < G->edge_set.size(); i++) {
        int u = G->edge_set[i].first.first;
        int v = G->edge_set[i].first.second;
        heap->insert(i, -(G->degree[u] + 1) * (G->in_degree[u] + 1) - (G->degree[v] + 1) * (G->in_degree[v] + 1));
    }
    while (heap->size() > 0) {
        int i = heap->pop();
        int u = G->edge_set[i].first.first;
        int v = G->edge_set[i].first.second;
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
    
    // Construct the index by vertex cover
    int i = 0;
    unsigned long long start_time = currentTime();
    for (auto it = vertex_cover.begin(); it != vertex_cover.end(); it++) {
        int u = *it;
        inv_vertex_cover[u] = i++;
        L.push_back(std::unordered_map<int, std::vector<std::pair<int, int>>>());
        std::queue<std::vector<int>> Q;
        std::vector<int> start;
        std::unordered_map<int, std::vector<int>> binary_indexed_tree;
        std::unordered_map<int, std::vector<std::pair<std::pair<int, int>, int>>> T;
        std::unordered_set<int> Vp;
        start.push_back(u);
        start.push_back(G->tmax + 1);
        start.push_back(-1);
        start.push_back(0);
        Q.push(start);

        while (!Q.empty()) {
            std::vector<int> current = Q.front();
            Q.pop();
            int v = current[0];
            if ((u == v && current[3] != 0)) {
                continue;
            }
            TemporalGraph::Edge* e = G->getHeadEdge(v);
            while (e) {
                int ts = std::min(current[1], e->interaction_time);
                int te = std::max(current[2], e->interaction_time);
                if (t_threshold == -1 || te - ts + 1 <= t_threshold) {
                    bool flag = false;
                    if (binary_indexed_tree.find(e->to) != binary_indexed_tree.end()) {
                        int t = ts + 1;
                        while (t <= G->tmax + 1) {
                            if (binary_indexed_tree[e->to][t] <= te) {
                                flag = true;
                                break;
                            }
                            t += (t & (-t));
                        }
                    }
                    if (!flag) {
                        if (Vp.find(e->to) == Vp.end() && vertex_cover.find(e->to) != vertex_cover.end()) {
                            Vp.insert(e->to);
                        }
                        if (binary_indexed_tree.find(e->to) == binary_indexed_tree.end()) {
                            T[e->to] = std::vector<std::pair<std::pair<int, int>, int>>();
                            binary_indexed_tree[e->to] = std::vector<int>();
                            binary_indexed_tree[e->to].assign(G->tmax + 2, G->tmax + 1);
                        }
                        T[e->to].push_back(std::make_pair(std::make_pair(ts, te), current[3] + 1));
                        int t = ts + 1;
                        while (t > 0) {
                            binary_indexed_tree[e->to][t] = std::min(binary_indexed_tree[e->to][t], te);
                            t -= (t & (-t));
                        }
                        if (current[3] + 1 < k) {
                            std::vector<int> next;
                            next.push_back(e->to);
                            next.push_back(ts);
                            next.push_back(te);
                            next.push_back(current[3] + 1);
                            Q.push(next);
                        }
                    }
                }
                e = e->next;
            }
        }

        for (auto it1 = Vp.begin(); it1 != Vp.end(); it1++) {
            L[i - 1][*it1] = std::vector<std::pair<int, int>>();
            cut[i - 1][*it1] = std::vector<int>();
            std::vector<std::pair<int, int>> intervals;
            std::vector<std::pair<std::pair<int, int>, int>>::iterator it2;
            for (it2 = T[*it1].begin(); it2 != T[*it1].end(); it2++) {
                if (it2->second <= k - 2) {
                    intervals.push_back(it2->first);
                }
                else {
                    break;
                }
            }
            std::sort(intervals.begin(), intervals.end(), cmp);
            int tmin = G->tmax + 1;
            for (auto it3 = intervals.begin(); it3 != intervals.end(); it3++) {
                if (tmin > it3->second) {
                    tmin = it3->second;
                    L[i - 1][*it1].push_back(*it3);
                }
            }
            cut[i - 1][*it1].push_back(L[i - 1][*it1].size());
            for (int j = k - 1; j <= k; j++) {
                intervals.clear();
                for (it2; it2 != T[*it1].end(); it2++) {
                    if (it2->second > j) {
                        break;
                    }
                    intervals.push_back(it2->first);
                    std::sort(intervals.begin(), intervals.end(), cmp);
                    int tmin = G->tmax + 1;
                    for (auto it3 = intervals.begin(); it3 != intervals.end(); it3++) {
                        if (tmin > it3->second && !reachable(G, u, *it1, it3->first, it3->second, j - 1)) {
                            tmin = it3->second;
                            L[i - 1][*it1].push_back(*it3);
                        }
                    }
                }
                cut[i - 1][*it1].push_back(L[i - 1][*it1].size());
            }
        }
        
        putProcess(double(i) / vertex_cover.size(), currentTime() - start_time);
    }
}

void Index::solve(TemporalGraph* G, char* query_file, char* output_file, int k) {
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