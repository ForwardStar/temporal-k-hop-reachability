#include "online.h"

bool cmp_t0(std::pair<std::pair<int, int>, int> i, std::pair<std::pair<int, int>, int> j) {
    return i.second < j.second;
}

std::string onlineSearch(TemporalGraph* Graph, int s, int t, int ts, int te, int k) {
    if (s == t) {
        return "Reachable";
    }

    auto G = Graph->projectedGraph(ts, te);

    std::vector<int> f;
    f.assign(G->n, G->n);
    f[s] = 0;
    std::vector<std::pair<std::pair<int, int>, int>> edges;
    std::queue<int> Q;
    Q.push(s);
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
                Q.push(e->to);
                f[e->to] = f[v] + 1;
            }
            e = G->getNextEdge(e);
        }
    }
    std::sort(edges.begin(), edges.end(), cmp_t0);
    
    f.assign(G->n, G->n);
    f[s] = 0;
    for (auto e : edges) {
        int u = e.first.first, v = e.first.second;
        f[v] = std::min(f[v], f[u] + 1);
    }

    if (f[t] <= k) {
        return "Reachable";
    }
    return "Not reachable";
}

void online(TemporalGraph* Graph, char* query_file, char* output_file) {
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
        fout << onlineSearch(Graph, s, t, ts, te, k) << std::endl;
        putProcess(double(++i) / query_num, currentTime() - start_time);
    }

    std::cout << "Average: " << timeFormatting((currentTime() - start_time) / query_num).str() << std::endl;
}