#include "online_search.h"

std::string onlineSearch(TemporalGraph* Graph, int s, int t, int ts, int te, int k) {
    if (s == t) {
        return "Reachable";
    }

    std::vector<bool> Vis(Graph->numOfVertices());
    for (int i = 0; i < Graph->numOfVertices(); i++) {
        Vis[i] = 0;
    }
    std::queue<std::pair<int, int>> Q;
    TemporalGraph* G = new TemporalGraph(Graph, ts, te);

    int dis = 0;
    Q.push(std::make_pair(s, 0));
    Vis[s] = true;
    while (!Q.empty()) {
        int u = Q.front().first;
        int dis = Q.front().second;
        if (dis >= k) {
            break;
        }
        Q.pop();
        TemporalGraph::Edge* edge = G->getHeadEdge(u);
        while (edge) {
            if (!Vis[edge->to]) {
                if (edge->to == t) {
                    return "Reachable";
                }
                Q.push(std::make_pair(edge->to, dis + 1));
                Vis[edge->to] = 1;
            }
            edge = G->getNextEdge(edge);
        }
    }

    delete G;
    return "Not reachable";
}

void online(TemporalGraph* Graph, char* query_file, char* output_file, int k) {
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
        fout << onlineSearch(Graph, s, t, ts, te, k) << std::endl;
        putProcess(double(++i) / query_num, currentTime() - start_time);
    }

    std::cout << "Average: " << timeFormatting((currentTime() - start_time) / query_num).str() << std::endl;
}