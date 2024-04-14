#include "commonfunctions.h"
#include "temporal_graph.h"

class NaiveIndex {
    private:

    public:
        std::unordered_set<int> vertex_cover;
        // L[ts][te][u][k]: on G_{[t_s, t_e]}, the k-hop neighbours of u.
        std::vector<std::vector<std::vector<std::unordered_map<int, int>>>> L;
        std::unordered_map<int, int> inv_vertex_cover;

        bool is_temporal_path = false;
        int k = 0;
        unsigned long long visited_paths = 0;
        unsigned long long max_number_of_paths = 0;

        unsigned long long size();

        bool reachable(TemporalGraph* G, int u, int v, int ts, int te, int k_input);

        NaiveIndex(TemporalGraph* G, int k_input, int t_threshold, std::string path_type);

        void solve(TemporalGraph* G, char* query_file, char* output_file);
};