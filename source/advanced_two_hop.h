#include "commonfunctions.h"
#include "temporal_graph.h"

class AdvancedTwoHopIndex {
    private:
        void construct_for_a_vertex(TemporalGraph* G, int u, bool reverse, std::vector<std::unordered_map<int, std::vector<std::vector<int>>>> &L, std::string algorithm);
        
    public:
        std::vector<std::unordered_map<int, std::vector<std::vector<int>>>> L_in, L_out;
        
        int k = 0;
        std::string index_construct_algorithm;
        int* order;

        bool reachable(TemporalGraph* G, int u, int v, int ts, int te, int d);

        AdvancedTwoHopIndex(TemporalGraph* G, int k_input, int t_threshold, std::string algorithm);
        ~AdvancedTwoHopIndex();
        
        void solve(TemporalGraph* G, char* query_file, char* output_file, int k);
};