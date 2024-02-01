#include "temporal_graph.h"

int TemporalGraph::numOfVertices() {
    return n;
}

int TemporalGraph::numOfEdges() {
    return m;
}

TemporalGraph::Edge* TemporalGraph::getHeadEdge(int u) {
    return head_edge[u];
}

TemporalGraph::Edge* TemporalGraph::getHeadInEdge(int u) {
    return head_in_edge[u];
}

TemporalGraph::Edge* TemporalGraph::getNextEdge(Edge* e) {
    return e->next;
}

TemporalGraph::Edge* TemporalGraph::deleteEdge(Edge* e) {
    int u = e->from;
    if (head_edge[u] == e) {
        head_edge[u] = e->next;
    }
    if (head_in_edge[u] == e) {
        head_in_edge[u] = e->next;
    }
    if (e->last) {
        e->last->next = e->next;
    }
    if (e->next) {
        e->next->last = e->last;
    }
    return e->next;
}

int TemporalGraph::getDestination(Edge* e) {
    return e->to;
}

int TemporalGraph::getInteractionTime(Edge* e) {
    return e->interaction_time;
}

void TemporalGraph::updateInfo() {
    m = 0;
    for (int u = 0; u < n; u++) {
        Edge* e = getHeadEdge(u);
        while (e) {
            m++;
            e = e->next;
        }
    }
}

void TemporalGraph::addEdge(int u, int v, int t, bool repeat) {
    m++;
    degree[u]++;
    if (head_edge[u]) {
        head_edge[u]->last = new Edge(u, v, t, head_edge[u]);
        head_edge[u] = head_edge[u]->last;
    }
    else {
        head_edge[u] = new Edge(u, v, t, nullptr);
    }
    in_degree[v]++;
    if (head_in_edge[v]) {
        head_in_edge[v]->last = new Edge(v, u, t, head_in_edge[v]);
        head_in_edge[v] = head_in_edge[v]->last;
    }
    else {
        head_in_edge[v] = new Edge(v, u, t, nullptr);
    }
    if (!is_directed && repeat) {
        m--;
        addEdge(v, u, t, false);
    }
}

TemporalGraph* TemporalGraph::projectedGraph(int ts, int te) {
    TemporalGraph* G = new TemporalGraph();
    G->n = n;
    G->m = 0;
    G->tmax = tmax;
    G->is_directed = is_directed;
    G->temporal_edge.resize(tmax + 1);
    G->head_edge.assign(n, nullptr);
    G->degree.assign(n, 0);
    G->head_in_edge.assign(n, nullptr);
    G->in_degree.assign(n, 0);

    for (int t = ts; t <= te; ++t) {
        for (auto it = temporal_edge[t].begin(); it != temporal_edge[t].end(); it++) {
            int u = it->first;
            int v = it->second;
            G->temporal_edge[t].push_back(std::make_pair(u, v));
            G->edge_set.push_back(std::make_pair(std::make_pair(u, v), t));
            G->addEdge(u, v, t);
        }
    }

    return G;
}

TemporalGraph::TemporalGraph(char* graph_file, std::string graph_type, double fraction) {
    int u, v, t;
    std::ifstream fin(graph_file);

    is_directed = graph_type == "Directed";

    while (fin >> u >> v >> t) {
        n = std::max(n, std::max(u, v) + 1);

        while (head_edge.size() < n) {
            head_edge.push_back(nullptr);
            degree.push_back(0);
        }

        while (head_in_edge.size() < n) {
            head_in_edge.push_back(nullptr);
            in_degree.push_back(0);
        }

        tmax = std::max(tmax, t);

        while (temporal_edge.size() < t + 1) {
            temporal_edge.push_back(std::vector<std::pair<int, int>>());
        }

        temporal_edge[t].push_back(std::make_pair(u, v));
        edge_set.push_back(std::make_pair(std::make_pair(u, v), t));
    }

    for (int t = 0; t <= tmax * fraction; t++) {
        for (auto e : temporal_edge[t]) {
            addEdge(e.first, e.second, t);
        }
    }
}

TemporalGraph::~TemporalGraph() {

    std::vector<Edge*>::iterator it;

    for (it = head_edge.begin(); it != head_edge.end(); it++) {
        Edge* temp = *it;
        while (temp) {
            Edge* next = temp->next;
            delete temp;
            temp = next;
        }
    }

    for (it = head_in_edge.begin(); it != head_in_edge.end(); it++) {
        Edge* temp = *it;
        while (temp) {
            Edge* next = temp->next;
            delete temp;
            temp = next;
        }
    }

}