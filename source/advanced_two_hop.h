#include "commonfunctions.h"
#include "temporal_graph.h"

class AdvancedTwoHopIndex {
    private:
        int find_lowest_index(std::vector<std::vector<int>> &L, int v);

        void construct_for_a_vertex(TemporalGraph* G, int u, bool reverse, 
                                    std::vector<std::vector<std::vector<int>>> &L, 
                                    std::vector<std::vector<int>> &next_idx, 
                                    std::vector<std::vector<int>> &cut_idx,
                                    int t_threshold,
                                    std::string algorithm);
        
    public:
        std::vector<std::vector<std::vector<int>>> L_in, L_out;
        std::vector<std::vector<int>> next_in, next_out;
        std::vector<std::vector<int>> cut_in, cut_out;
        
        int k = 0;
        std::string index_construct_algorithm;
        int* order;

        int size();

        bool reachable(TemporalGraph* G, int u, int v, int ts, int te, int d);

        AdvancedTwoHopIndex(TemporalGraph* G, int k_input, int t_threshold, std::string algorithm);
        ~AdvancedTwoHopIndex();
        
        void solve(TemporalGraph* G, char* query_file, char* output_file, int k);
};