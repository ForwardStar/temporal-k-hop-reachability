#include "temporal_graph.h"

void TemporalGraph::tarjan(int now, int &t) {
    dfsOrder[now] = ++t;
    lowestOrder[now] = t;
    Vis[now] = true;
    Stack.push(now);

    TemporalGraph::Edge* edge = getHeadEdge(now);
    
    while (edge) {
        if (!Vis[edge->to]) {
            tarjan(edge->to, t);
        }
        if (!outOfStack[edge->to]) {
            lowestOrder[now] = std::min(lowestOrder[now], lowestOrder[edge->to]);
        }
        edge = edge->next;
    }

    if (dfsOrder[now] == lowestOrder[now]) {
        std::vector<int> CurrentSCC;
        while (Stack.top() != now) {
            outOfStack[Stack.top()] = true;
            CurrentSCC.push_back(Stack.top());
            Stack.pop();
        }
        outOfStack[Stack.top()] = true;
        CurrentSCC.push_back(Stack.top());
        Stack.pop();
        AllSCC.push_back(CurrentSCC);
    }
}

std::vector<std::vector<int>> TemporalGraph::findSCC() {
    dfsOrder.resize(n);
    lowestOrder.resize(n);
    outOfStack.assign(n, 0);
    Vis.assign(n, 0);
    AllSCC.clear();

    int t = 0;
    for (int u = 0; u < n; u++) {
        if (!Vis[u]) {
            tarjan(u, t);
        }
    }
    
    return AllSCC;
}

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

int TemporalGraph::getDestination(Edge* e) {
    return e->to;
}

int TemporalGraph::getInteractionTime(Edge* e) {
    return e->interaction_time;
}

void TemporalGraph::addEdge(int u, int v, int t, bool repeat) {
    m++;
    degree[u]++;
    if (head_edge[u]) {
        head_edge[u] = new Edge(v, t, head_edge[u]);
    }
    else {
        head_edge[u] = new Edge(v, t, nullptr);
    }
    in_degree[v]++;
    if (head_in_edge[v]) {
        head_in_edge[v] = new Edge(u, t, head_in_edge[v]);
    }
    else {
        head_in_edge[v] = new Edge(u, t, nullptr);
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

TemporalGraph::TemporalGraph(char* graph_file, char* graph_type) {
    int u, v, t;
    std::ifstream fin(graph_file);

    is_directed = std::strcmp(graph_type, "Directed") == 0;

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
        addEdge(u, v, t);
    }
    ++n;
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