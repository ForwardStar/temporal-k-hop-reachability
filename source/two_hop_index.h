#include "commonfunctions.h"
#include "temporal_graph.h"

class TwoHopIndex {
    private:

    public:
        std::unordered_set<int> vertex_cover;
        std::vector<std::unordered_map<int, std::vector<std::pair<int, int>>>> L;
        std::vector<std::unordered_map<int, std::vector<int>>> cut;
        std::unordered_map<int, int> inv_vertex_cover;

        int k = 0;
        std::string index_construct_algorithm;

        int size();

        bool reachable(TemporalGraph* G, int u, int v, int ts, int te, int k_input);

        TwoHopIndex(TemporalGraph* G, int k_input, int t_threshold, std::string algorithm);

        void solve(TemporalGraph* G, char* query_file, char* output_file, int k);
};