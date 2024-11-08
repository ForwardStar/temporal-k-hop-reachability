#include "online1.h"

std::string onlineSearch1(TemporalGraph* Graph, int s, int t, int ts, int te, int k) {
    if (s == t) {
        return "Reachable";
    }

    std::vector<int> T;
    T.assign(Graph->n, 2147483647);
    T[s] = -1;
    std::queue<std::vector<int>> Q;
    TemporalGraph* G = Graph->projectedGraph(ts, te);

    Q.push(std::vector<int>{s, 0, -1});
    while (!Q.empty()) {
        int u = Q.front()[0];
        int dis = Q.front()[1];
        int te = Q.front()[2];
        Q.pop();
        if (dis >= k) {
            break;
        }
        for (auto e : G->neighbors[u]) {
            if (e.second >= T[e.first] || te > e.second) {
                continue;
            }
            if (e.first == t) {
                delete G;
                return "Reachable";
            }
            T[e.first] = std::min(T[e.first], e.second);
            Q.push(std::vector<int>{e.first, dis + 1, e.second});
        }
    }

    delete G;
    return "Not reachable";
}

void online1(TemporalGraph* Graph, char* query_file, char* output_file) {
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
        fout << onlineSearch1(Graph, s, t, ts, te, k) << std::endl;
        putProcess(double(++i) / query_num, currentTime() - start_time);
    }

    std::cout << "Average: " << timeFormatting((currentTime() - start_time) / query_num).str() << std::endl;
}