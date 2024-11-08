#include "online2.h"

bool cmp_t0(std::pair<std::pair<int, int>, int> i, std::pair<std::pair<int, int>, int> j) {
    return i.second < j.second;
}

std::string onlineSearch2(TemporalGraph* Graph, int s, int t, int ts, int te, int k) {
    if (s == t) {
        return "Reachable";
    }

    std::vector<int> f;
    std::unordered_set<int> sources;
    std::vector<std::pair<std::pair<int, int>, int>> edges;
    auto shortest_first = [](std::pair<int, int> i, std::pair<int, int> j) {
        return i.second > j.second;
    };
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, decltype(shortest_first)> Q(shortest_first);
    TemporalGraph* G = new TemporalGraph(Graph->n);
    f.assign(Graph->n, Graph->n);
    f[s] = 0;
    for (int i = ts; i <= te; i++) {
        for (auto e : Graph->temporal_edge[i]) {
            sources.insert(e.first);
            G->addEdge(e.first, e.second, i);
        }
        for (int u : sources) {
            Q.push(std::make_pair(u, f[u]));
        }
        sources.clear();
        while (!Q.empty()) {
            int u = Q.top().first;
            Q.pop();
            for (auto e : G->neighbors[u]) {
                if (f[u] + 1 < f[e.first]) {
                    f[e.first] = f[u] + 1;
                    Q.push(std::make_pair(e.first, f[e.first]));
                }
            }
        }
        for (auto e : Graph->temporal_edge[i]) {
            G->neighbors[e.first].clear();
            G->in_neighbors[e.second].clear();
            G->degree[e.first] = 0;
            G->in_degree[e.second] = 0;
        }
    }
    delete G;

    if (f[t] <= k) {
        return "Reachable";
    }
    return "Not reachable";
}

void online2(TemporalGraph* Graph, char* query_file, char* output_file) {
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
        fout << onlineSearch2(Graph, s, t, ts, te, k) << std::endl;
        putProcess(double(++i) / query_num, currentTime() - start_time);
    }

    std::cout << "Average: " << timeFormatting((currentTime() - start_time) / query_num).str() << std::endl;
}