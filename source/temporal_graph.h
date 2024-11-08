#ifndef TEMPORALGRAPH
#define TEMPORALGRAPH

#include "commonfunctions.h"

class TemporalGraph {
    public:
        // n: the number of vertices; m: the number of edges.
        int n = 0;
        int m = 0;

        // tmax: the maximum time of all temporal edges.
        int tmax;

        // degree: the out-degree of all vertices.
        std::vector<int> degree;

        // in_degree: the in-degree of all vertices.
        std::vector<int> in_degree;

        // edge_set: all edges.
        std::vector<std::pair<std::pair<int, int>, int>> edge_set;

        std::vector<std::vector<std::pair<int, int>>> neighbors, in_neighbors;

        // temporal_edge[t] --> the edge set at time t.
        std::vector<std::vector<std::pair<int, int>>> temporal_edge;

        // findSCC(): find all SCCs in the graph.
        std::vector<std::vector<int>> findSCC();

        // numOfVertices(): get the number of the vertices in the graph.
        int numOfVertices();

        // numOfEdges(): get the number of the edges in the graph.  
        int numOfEdges();

        // updateInfo(): update the number of edges.
        void updateInfo();

        // addInEdge(u, v, t): add only an in-edge (v, u, t) to the graph.
        void addInEdge(int u, int v, int t);

        // addOutEdge(u, v, t): add only an out-edge (u, v, t) to the graph.
        void addOutEdge(int u, int v, int t);

        // addEdge(u, v, t): add an edge (u, v, t) to the graph.
        void addEdge(int u, int v, int t);
        
        // inducedSubgraph(S): return the induced subgraph of S.
        template <typename T>
        TemporalGraph* inducedSubgraph(T S) {
            TemporalGraph* G = new TemporalGraph();
            G->n = n;
            G->m = 0;
            G->tmax = tmax;

            for (auto it = S.begin(); it != S.end(); it++) {
                for (auto e : *it) {
                    temporal_edge[e.second].push_back(std::make_pair(*it, e.first));
                    edge_set.push_back(std::make_pair(std::make_pair(*it, e.first), e.second));
                    G->addEdge(*it, e.first, e.second);
                }
            }

            return G;
        }

        // projectedGraph(ts, te): return the projected graph of [ts, te].
        TemporalGraph* projectedGraph(int ts, int te);

        TemporalGraph() {}
        TemporalGraph(int n_input);
        TemporalGraph(char* graph_file, double fraction);
        ~TemporalGraph() {}
};

#endif