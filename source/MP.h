#include "commonfunctions.h"
#include "temporal_graph.h"

class MPIndex {
    private:
        std::vector<int> f;
        std::unordered_set<int> visited_vertices;

    public:
        std::unordered_set<int> vertex_cover;
        std::vector<std::unordered_map<int, std::vector<std::vector<std::pair<int, int>>>>> L;
        std::unordered_map<int, int> inv_vertex_cover;

        int k = 0;
        unsigned long long max_number_of_paths = 0;

        unsigned long long size();

        // Note: k_input >= 2.
        bool reachable(TemporalGraph* G, int u, int v, int ts, int te, int k_input);

        MPIndex(TemporalGraph* G, int k_input);

        void solve(TemporalGraph* G, char* query_file, char* output_file);
};