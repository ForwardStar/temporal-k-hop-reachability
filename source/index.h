#include "commonfunctions.h"
#include "temporal_graph.h"
#include "heap.h"

class Index {
    private:

    public:
        // vertex_cover: vertices' ids of the vertex_cover;
        std::unordered_set<int> vertex_cover;
        // index[inv_vertex_cover[i]]: all potential (k - 2)-hop vertices of vertex i;
        // index[inv_vertex_cover[i]][j]: j is a potential (k - 2)-hop vertex of vertex i;
        // index[inv_vertex_cover[i]][j]:
        //      it->first: shortest distance (k - 2, k - 1 or k) from i to j;
        //      it->second: corresponding minimal time interval of the shortest distance;
        std::vector<std::unordered_map<int, std::vector<std::set<long long>>>> index;
        std::unordered_map<int, int> inv_vertex_cover;

        int k = 0;

        void shrink_to_fit(int u);

        bool reachable(TemporalGraph* G, int u, int v, int ts, int te, int k_input);

        Index(TemporalGraph* G, int k_input, int t_threshold);

        void solve(TemporalGraph* Graph, char* query_file, char* output_file, int k);
};