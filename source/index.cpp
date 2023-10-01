#include "index.h"

bool cmp(std::pair<int, int> i, std::pair<int, int> j) {
    return i.first > j.first || (i.first == j.first && i.second < j.second);
}

bool cmp1(std::pair<std::pair<int, int>, int> i, std::pair<std::pair<int, int>, int> j) {
    return i.second < j.second;
}

int Index::size() {
    int num_intervals = 0;
    for (auto it = L.begin(); it != L.end(); it++) {
        for (auto it1 = it->begin(); it1 != it->end(); it1++) {
            num_intervals += it1->second.size();
        }
    }
    return num_intervals;
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
            if (L[i].find(v) == L[i].end()) {
                return false;
            }
            if (index_construct_algorithm.rfind("PrioritySearch", 0) == 0) {
                for (int j = k - 2; j <= k_input; j++) {
                    int l = 0;
                    if (j > k - 2) {
                        l = cut[i][v][j - 1 - (k - 2)];
                    }
                    int r = cut[i][v][j - (k - 2)];
                    for (l; l < r; l++) {
                        if (L[i][v][l].first >= ts && L[i][v][l].second <= te) {
                            return true;
                        }
                    }
                }
            }
            else {
                for (int j = k - 2; j <= k_input; j++) {
                    int l = 0;
                    if (j > k - 2) {
                        l = cut[i][v][j - 1 - (k - 2)];
                    }
                    int r = cut[i][v][j - (k - 2)] - 1;
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

Index::Index(TemporalGraph* G, int k_input, int t_threshold, std::string algorithm) {
    k = k_input;
    index_construct_algorithm = algorithm;

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
    
    // Construct the index by vertex cover
    if (algorithm.rfind("PrioritySearch", 0) == 0) {
        int i = 0;
        unsigned long long start_time = currentTime();
        for (auto it = vertex_cover.begin(); it != vertex_cover.end(); it++) {
            int u = *it;
            inv_vertex_cover[u] = i++;
            L.push_back(std::unordered_map<int, std::vector<std::pair<int, int>>>());
            auto smaller_interval_first = [](std::vector<int> i, std::vector<int> j) {
                return (i[2] - i[1] > j[2] - j[1]) || (i[2] - i[1] == j[2] - j[1] && i[3] > j[3]);
            };
            std::priority_queue<std::vector<int>, std::vector<std::vector<int>>, decltype(smaller_interval_first)> Q(smaller_interval_first);
            std::vector<std::vector<std::vector<int>>> binary_indexed_tree;
            std::vector<std::vector<std::vector<std::pair<int, int>>>> T;
            binary_indexed_tree.resize(G->n);
            T.resize(G->n);
            for (int u = 0; u < G->n; u++) {
                T[u].resize(k + 1);
            }
            std::vector<int> start;
            std::unordered_set<int> Vp;
            start.push_back(u);
            start.push_back(G->tmax + 1);
            start.push_back(-1);
            start.push_back(0);
            Q.push(start);

            while (!Q.empty()) {
                std::vector<int> current = Q.top();
                Q.pop();
                int v = current[0];
                if ((u == v && current[3] != 0)) {
                    continue;
                }
                bool flag = false;
                for (int j = 0; j <= current[3]; j++) {
                    for (auto it1 = T[v][j].begin(); it1 != T[v][j].end(); it1++) {
                        if (it1->first >= current[1] && it1->second <= current[2]) {
                            flag = true;
                            break;
                        }
                    }
                    if (flag) {
                        break;
                    }
                }
                if (flag) {
                    continue;
                }
                if (u != v && Vp.find(v) == Vp.end() && vertex_cover.find(v) != vertex_cover.end()) {
                    Vp.insert(v);
                }
                T[v][current[3]].push_back(std::make_pair(current[1], current[2]));
                if (current[3] == k) {
                    continue;
                }
                TemporalGraph::Edge* e = G->getHeadEdge(v);
                while (e) {
                    int ts = std::min(current[1], e->interaction_time);
                    int te = std::max(current[2], e->interaction_time);
                    if (t_threshold == -1 || te - ts + 1 <= t_threshold) {
                        if (algorithm == "PrioritySearch-naive") {
                            bool flag = false;
                            for (int j = 0; j <= current[3] + 1; j++) {
                                for (auto it1 = T[e->to][j].begin(); it1 != T[e->to][j].end(); it1++) {
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
                                Q.push(next);
                            }
                        }
                        else {
                            bool flag = false;
                            if (binary_indexed_tree[e->to].size() > 0) {
                                int t = ts + 1;
                                while (t <= G->tmax + 1) {
                                    if (binary_indexed_tree[e->to][t].size() > 0) {
                                        int d = current[3] + 1;
                                        while (d > 0) {
                                            if (binary_indexed_tree[e->to][t][d] <= te) {
                                                flag = true;
                                                break;
                                            }
                                            d -= (d & (-d));
                                        }
                                    }
                                    if (flag) {
                                        break;
                                    }
                                    t += (t & (-t));
                                }
                            }
                            if (!flag) {
                                if (binary_indexed_tree[e->to].size() == 0) {
                                    binary_indexed_tree[e->to].assign(G->tmax + 2, std::vector<int>());
                                }
                                int t = ts + 1;
                                while (t > 0) {
                                    if (binary_indexed_tree[e->to][t].size() == 0) {
                                        binary_indexed_tree[e->to][t].assign(k + 1, G->tmax + 1);
                                    }
                                    int d = current[3] + 1;
                                    while (d <= k) {
                                        binary_indexed_tree[e->to][t][d] = std::min(binary_indexed_tree[e->to][t][d], te);
                                        d += (d & (-d));
                                    }
                                    t -= (t & (-t));
                                }
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

            if (algorithm == "PrioritySearch-naive") {
                for (auto it1 = Vp.begin(); it1 != Vp.end(); it1++) {
                    L[i - 1][*it1] = std::vector<std::pair<int, int>>();
                    cut[i - 1][*it1] = std::vector<int>();
                    for (int j = k - 2; j >= 0; --j) {
                        for (auto it2 = T[*it1][j].begin(); it2 != T[*it1][j].end(); it2++) {
                            bool flag = false;
                            for (auto it3 = L[i - 1][*it1].begin(); it3 != L[i - 1][*it1].end(); it3++) {
                                if (it3->first >= it2->first && it3->second <= it2->second) {
                                    flag = true;
                                    break;
                                }
                            }
                            if (!flag) {
                                L[i - 1][*it1].push_back(*it2);
                            }
                        }
                    }
                    cut[i - 1][*it1].push_back(L[i - 1][*it1].size());
                    for (int j = k - 1; j <= k; j++) {
                        for (auto it2 = T[*it1][j].begin(); it2 != T[*it1][j].end(); it2++) {
                            L[i - 1][*it1].push_back(*it2);
                        }
                        cut[i - 1][*it1].push_back(L[i - 1][*it1].size());
                    }
                }
            }

            putProcess(double(i) / vertex_cover.size(), currentTime() - start_time);
        }
    }
    else if (algorithm.rfind("BFS", 0) == 0) {
        int i = 0;
        unsigned long long start_time = currentTime();
        for (auto it = vertex_cover.begin(); it != vertex_cover.end(); it++) {
            int u = *it;
            inv_vertex_cover[u] = i++;
            L.push_back(std::unordered_map<int, std::vector<std::pair<int, int>>>());
            std::queue<std::vector<int>> Q;
            std::vector<int> start;
            std::vector<std::vector<std::pair<int, int>>> binary_indexed_tree;
            std::vector<std::vector<std::pair<std::pair<int, int>, int>>> T;
            binary_indexed_tree.resize(G->n);
            T.resize(G->n);
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
                // Adopts double-checking to avoid expanding non-minimal paths
                if (algorithm == "BFS-naive") {
                    bool flag = false;
                    for (auto it1 = T[v].begin(); it1 != T[v].end(); it1++) {
                        if (it1->first.first >= current[1] && it1->first.second <= current[2] && it1->second == current[3] && (it1->first.first != current[1] || it1->first.second != current[2])) {
                            flag = true;
                            break;
                        }
                    }
                    if (flag) {
                        continue;
                    }
                }
                else {
                    if (binary_indexed_tree[v].size() > 0) {
                        bool flag = false;
                        int t = current[1] + 1;
                        while (t <= G->tmax + 1) {
                            if (binary_indexed_tree[v][t].first < current[2] && binary_indexed_tree[v][t].second == current[3]) {
                                flag = true;
                                break;
                            }
                            t += (t & (-t));
                        }
                        if (flag) {
                            continue;
                        }
                    }
                }
                if ((u == v && current[3] != 0)) {
                    continue;
                }
                TemporalGraph::Edge* e = G->getHeadEdge(v);
                while (e) {
                    int ts = std::min(current[1], e->interaction_time);
                    int te = std::max(current[2], e->interaction_time);
                    if (t_threshold == -1 || te - ts + 1 <= t_threshold) {
                        // BFS-naive: check minimality by enumerating all previous paths
                        if (algorithm == "BFS-naive") {
                            bool flag = false;
                            for (auto it1 = T[e->to].begin(); it1 != T[e->to].end(); it1++) {
                                if (it1->first.first >= ts && it1->first.second <= te) {
                                    flag = true;
                                    break;
                                }
                            }
                            if (!flag) {
                                if (e->to != u && Vp.find(e->to) == Vp.end() && vertex_cover.find(e->to) != vertex_cover.end()) {
                                    Vp.insert(e->to);
                                }
                                T[e->to].push_back(std::make_pair(std::make_pair(ts, te), current[3] + 1));
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
                        // BFS-full: check minimality by binary indexed tree
                        else {
                            bool flag = false;
                            if (binary_indexed_tree[e->to].size() > 0) {
                                int t = ts + 1;
                                while (t <= G->tmax + 1) {
                                    if (binary_indexed_tree[e->to][t].first <= te) {
                                        flag = true;
                                        break;
                                    }
                                    t += (t & (-t));
                                }
                            }
                            if (!flag) {
                                if (e->to != u && Vp.find(e->to) == Vp.end() && vertex_cover.find(e->to) != vertex_cover.end()) {
                                    Vp.insert(e->to);
                                }
                                if (binary_indexed_tree[e->to].size() == 0) {
                                    binary_indexed_tree[e->to].assign(G->tmax + 2, std::make_pair(G->tmax + 1, 0));
                                }
                                T[e->to].push_back(std::make_pair(std::make_pair(ts, te), current[3] + 1));
                                int t = ts + 1;
                                while (t > 0) {
                                    if (te < binary_indexed_tree[e->to][t].first) {
                                        binary_indexed_tree[e->to][t].first = te;
                                        binary_indexed_tree[e->to][t].second = current[3] + 1;
                                    }
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
                    }
                    std::sort(intervals.begin(), intervals.end(), cmp);
                    int tmin = G->tmax + 1;
                    for (auto it3 = intervals.begin(); it3 != intervals.end(); it3++) {
                        if (tmin > it3->second && !reachable(G, u, *it1, it3->first, it3->second, j - 1)) {
                            tmin = it3->second;
                            L[i - 1][*it1].push_back(*it3);
                        }
                    }
                    cut[i - 1][*it1].push_back(L[i - 1][*it1].size());
                }
            }
            
            putProcess(double(i) / vertex_cover.size(), currentTime() - start_time);
        }
    }
    else if (algorithm == "Experimental") {
        int i = 0;
        unsigned long long start_time = currentTime();
        // Reconstruct the graph by vertex cover
        TemporalGraph* Gs = G->inducedSubgraph(vertex_cover);
        std::vector<std::vector<int>> SCCs = Gs->findSCC();
        std::vector<std::vector<int>> last_intervals;
        std::vector<std::vector<int>> cur_intervals;
        last_intervals.resize(G->n);
        cur_intervals.resize(G->n);
        for (auto SCC : SCCs) {
            if (SCC.size() == 1 && vertex_cover.find(SCC[0]) == vertex_cover.end()) {
                continue;
            }
            for (auto u : SCC) {
                inv_vertex_cover[u] = i++;
                L.push_back(std::unordered_map<int, std::vector<std::pair<int, int>>>());
                std::queue<std::vector<int>> Q;
                std::vector<int> start;
                std::vector<std::vector<std::pair<int, int>>> binary_indexed_tree;
                std::vector<std::vector<std::pair<std::pair<int, int>, int>>> T;
                binary_indexed_tree.resize(G->n);
                T.resize(G->n);
                std::unordered_set<int> Vp;
                start.push_back(u);
                start.push_back(G->tmax + 1);
                start.push_back(-1);
                start.push_back(0);
                Q.push(start);
                int last_dis = 0;
                cur_intervals.clear();
                cur_intervals.resize(G->n);
                auto it = last_intervals.begin();

                while (!Q.empty()) {
                    std::vector<int> current = Q.front();
                    Q.pop();
                    int v = current[0];
                    int ts = current[1];
                    int te = current[2];
                    // Adopts double-checking to avoid expanding non-minimal paths
                    if (binary_indexed_tree[v].size() > 0) {
                        bool flag = false;
                        int t = current[1] + 1;
                        while (t > 0) {
                            if (binary_indexed_tree[v][t].first < current[2] && binary_indexed_tree[v][t].second == current[3]) {
                                flag = true;
                                break;
                            }
                            t -= (t & (-t));
                        }
                        if (flag) {
                            continue;
                        }
                    }
                    std::vector<int> temp;
                    temp.push_back(v);
                    temp.push_back(ts);
                    temp.push_back(te);
                    temp.push_back(current[3]);
                    cur_intervals.push_back(temp);
                }
            }
        }
        delete Gs;
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