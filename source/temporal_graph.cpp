#include "temporal_graph.h"

int TemporalGraph::numOfVertices() {
    return n;
}

int TemporalGraph::numOfEdges() {
    return m;
}

void TemporalGraph::addInEdge(int u, int v, int t) {
    in_degree[u]++;
    in_neighbors[u].push_back({v, t});
}

void TemporalGraph::addOutEdge(int u, int v, int t) {
    degree[u]++;
    neighbors[u].push_back({v, t});
}

void TemporalGraph::addEdge(int u, int v, int t) {
    m++;
    addInEdge(v, u, t);
    addOutEdge(u, v, t);
}

TemporalGraph* TemporalGraph::projectedGraph(int ts, int te) {
    TemporalGraph* G = new TemporalGraph();
    G->n = n;
    G->m = 0;
    G->tmax = tmax;
    G->temporal_edge.resize(tmax + 1);
    G->degree.assign(n, 0);
    G->in_degree.assign(n, 0);
    G->neighbors.resize(n);
    G->in_neighbors.resize(n);

    for (int t = ts; t <= te; ++t) {
        for (auto it = temporal_edge[t].begin(); it != temporal_edge[t].end(); it++) {
            int u = it->first;
            int v = it->second;
            G->temporal_edge[t].emplace_back(std::make_pair(u, v));
            G->edge_set.emplace_back(std::make_pair(std::make_pair(u, v), t));
            G->addEdge(u, v, t);
        }
    }

    return G;
}

TemporalGraph::TemporalGraph(int n_input) {
    n = n_input;
    degree.assign(n, 0);
    in_degree.assign(n, 0);
    neighbors.resize(n);
    in_neighbors.resize(n);
}

TemporalGraph::TemporalGraph(char* graph_file, double fraction) {
    int u, v, t;
    std::ifstream fin(graph_file);

    while (fin >> u >> v >> t) {
        n = std::max(n, std::max(u, v) + 1);

        while (neighbors.size() < n) {
            neighbors.emplace_back(std::vector<std::pair<int, int>>());
            degree.emplace_back(0);
            in_neighbors.emplace_back(std::vector<std::pair<int, int>>());
            in_degree.emplace_back(0);
        }

        tmax = std::max(tmax, t);

        while (temporal_edge.size() < t + 1) {
            temporal_edge.emplace_back(std::vector<std::pair<int, int>>());
        }

        temporal_edge[t].emplace_back(std::make_pair(u, v));
        edge_set.emplace_back(std::make_pair(std::make_pair(u, v), t));
    }

    for (int t = 0; t <= tmax * fraction; t++) {
        for (auto e : temporal_edge[t]) {
            addEdge(e.first, e.second, t);
        }
    }
}