#ifndef TEMPORALGRAPH
#define TEMPORALGRAPH

#include "commonfunctions.h"

class TemporalGraph {
    public:
        // Edge: the structure of the edges in the temporal graph.
        struct Edge {
            // to: the destination of the edge.
            int from, to;
            // interaction_time: time of the interaction.
            int interaction_time;
            // next, last: the next/last edge in linked list structure.
            Edge* next;
            Edge* last;
            Edge(int u, int v, int t, Edge* nextptr): from(u), to(v), interaction_time(t), next(nextptr), last(nullptr) {}
            ~Edge() {}
        };

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

        // head_edge[vertex] --> the head out-edge from this vertex.
        std::vector<Edge*> head_edge;

        // inv_head_edge[vertex] --> the head in-edge from this vertex.
        std::vector<Edge*> head_in_edge;

        // temporal_edge[t] --> the edge set at time t.
        std::vector<std::vector<std::pair<int, int>>> temporal_edge;

        // findSCC(): find all SCCs in the graph.
        std::vector<std::vector<int>> findSCC();

        // numOfVertices(): get the number of the vertices in the graph.
        int numOfVertices();

        // numOfEdges(): get the number of the edges in the graph.  
        int numOfEdges();

        // getHeadEdge(u): get the head out-edge of vertex u.
        Edge* getHeadEdge(int u);

        // getHeadInEdge(u): get the head in-edge of vertex u.
        Edge* getHeadInEdge(int u);

        // getNextEdge(e): get the next edge of edge e.
        Edge* getNextEdge(Edge* e);

        // deleteEdge(e): delete e and return the next edge of e.
        Edge* deleteEdge(Edge* e);

        // getDestination(e): get the destination of the edge e.
        int getDestination(Edge* e);

        // getInteractionTime(e): get the time of the interaction.
        int getInteractionTime(Edge* e);

        // updateInfo(): update the number of edges.
        void updateInfo();

        // addEdge(u, v, t): add an edge (u, v, t) to the graph.
        void addEdge(int u, int v, int t, bool repeat=true);
        
        // inducedSubgraph(S): return the induced subgraph of S.
        template <typename T>
        TemporalGraph* inducedSubgraph(T S) {
            TemporalGraph* G = new TemporalGraph();
            G->n = n;
            G->m = 0;
            G->tmax = tmax;

            while (G->head_edge.size() < G->n) {
                head_edge.push_back(nullptr);
                degree.push_back(0);
            }

            while (G->head_in_edge.size() < G->n) {
                head_in_edge.push_back(nullptr);
                in_degree.push_back(0);
            }

            for (auto it = S.begin(); it != S.end(); it++) {
                Edge* e = getHeadEdge(*it);
                while (e) {
                    temporal_edge[e->interaction_time].push_back(std::make_pair(*it, e->to));
                    edge_set.push_back(std::make_pair(std::make_pair(*it, e->to), e->interaction_time));
                    G->addEdge(*it, e->to, e->interaction_time);
                    e = e->next;
                }
            }

            return G;
        }

        // projectedGraph(ts, te): return the projected graph of [ts, te].
        TemporalGraph* projectedGraph(int ts, int te);

        TemporalGraph() {}
        TemporalGraph(int n_input);
        TemporalGraph(char* graph_file, double fraction);
        ~TemporalGraph();
};

#endif