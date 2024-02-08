#include "commonfunctions.h"
#include "temporal_graph.h"

class BaselineIndex {
    private:

    public:
        std::unordered_set<int> vertex_cover;
        std::vector<std::unordered_map<int, std::vector<std::pair<int, int>>>> L;
        std::vector<std::unordered_map<int, std::vector<int>>> cut;
        std::unordered_map<int, int> inv_vertex_cover;

        bool is_temporal_path = false;
        int k = 0;
        unsigned long long visited_paths = 0;

        unsigned long long size();

        bool reachable(TemporalGraph* G, int u, int v, int ts, int te, int k_input);

        BaselineIndex(TemporalGraph* G, int k_input, int t_threshold, std::string path_type);

        void solve(TemporalGraph* G, char* query_file, char* output_file);
};