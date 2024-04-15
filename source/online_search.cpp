#include "online_search.h"

std::string onlineSearch(TemporalGraph* Graph, int s, int t, int ts, int te, int k, std::string path_type) {
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
        TemporalGraph::Edge* edge = G->getHeadEdge(u);
        while (edge) {
            if (path_type == "Temporal" && (edge->interaction_time >= T[edge->to] || te > edge->interaction_time)) {
                edge = G->getNextEdge(edge);
                continue;
            }
            if (path_type == "Temporal" || T[edge->to] == 2147483647) {
                if (edge->to == t) {
                    delete G;
                    return "Reachable";
                }
                T[edge->to] = std::min(T[edge->to], edge->interaction_time);
                Q.push(std::vector<int>{edge->to, dis + 1, edge->interaction_time});
            }
            edge = G->getNextEdge(edge);
        }
    }

    delete G;
    return "Not reachable";
}

void online(TemporalGraph* Graph, char* query_file, char* output_file, std::string path_type) {
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
        fout << onlineSearch(Graph, s, t, ts, te, k, path_type) << std::endl;
        putProcess(double(++i) / query_num, currentTime() - start_time);
    }

    std::cout << "Average: " << timeFormatting((currentTime() - start_time) / query_num).str() << std::endl;
}