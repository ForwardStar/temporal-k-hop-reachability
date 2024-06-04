#include "online_search.h"

std::string onlineSearch(TemporalGraph* Graph, int s, int t, int ts, int te, int k) {
    if (s == t) {
        return "Reachable";
    }
    
    std::vector<int> f;
    f.assign(Graph->n, Graph->n);
    f[s] = 0;
    for (int i = ts; i <= te; i++) {
        for (auto e : Graph->temporal_edge[i]) {
            int u = e.first, v = e.second;
            f[v] = std::min(f[v], f[u] + 1);
        }
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
        // Perform online BFS Search
        fout << onlineSearch(Graph, s, t, ts, te, k) << std::endl;
        putProcess(double(++i) / query_num, currentTime() - start_time);
    }

    std::cout << "Average: " << timeFormatting((currentTime() - start_time) / query_num).str() << std::endl;
}