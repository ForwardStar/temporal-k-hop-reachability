#include "commonfunctions.h"
#include "temporal_graph.h"

class MPIndex {
    private:

    public:
        std::vector<std::unordered_map<int, std::vector<std::vector<std::pair<int, int>>>>> L;

        int k = 0;
        unsigned long long alpha = 0;

        unsigned long long size();

        bool reachable(TemporalGraph* G, int u, int v, int ts, int te, int k_input);

        MPIndex(TemporalGraph* G, int k_input);

        void solve(TemporalGraph* G, char* query_file, char* output_file);
};