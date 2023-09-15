#include "commonfunctions.h"
#include "temporal_graph.h"
#include "heap.h"

class Index {
    private:

    public:
        std::unordered_set<int> vertex_cover;
        std::vector<std::unordered_map<int, std::vector<std::pair<int, int>>>> L;
        std::vector<std::unordered_map<int, std::vector<int>>> cut;
        std::unordered_map<int, int> inv_vertex_cover;

        int k = 0;

        bool reachable(TemporalGraph* G, int u, int v, int ts, int te, int k_input);

        Index(TemporalGraph* G, int k_input, int t_threshold);

        void solve(TemporalGraph* G, char* query_file, char* output_file, int k);
};