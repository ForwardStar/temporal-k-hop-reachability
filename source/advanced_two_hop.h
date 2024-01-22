#include "commonfunctions.h"
#include "temporal_graph.h"

class AdvancedTwoHopIndex {
    private:
        std::unordered_set<int> affected_vertices;
        std::vector<std::vector<std::vector<int>>> inc_index;
        std::vector<std::vector<std::pair<int, int>>> binary_indexed_tree;
        std::queue<std::vector<int>> Q;

        int find_lowest_index(std::vector<std::vector<int>> &L, int v);

        void construct_for_a_vertex(TemporalGraph* G, int u, bool reverse, 
                                    std::vector<std::vector<std::vector<int>>> &L, 
                                    std::vector<std::vector<int>> &next_idx, 
                                    std::vector<std::vector<int>> &cut_idx,
                                    int t_threshold);
        
    public:
        std::vector<std::vector<std::vector<int>>> L_in, L_out;
        std::vector<std::vector<int>> next_in, next_out;
        std::vector<std::vector<int>> cut_in, cut_out;
        
        bool is_temporal_path;
        int k = 0;
        int visited_paths = 0;
        int* order;

        int size();

        bool reachable(TemporalGraph* G, int u, int v, int ts, int te, int d);

        AdvancedTwoHopIndex(TemporalGraph* G, int k_input, int t_threshold, std::string path_type);
        ~AdvancedTwoHopIndex();
        
        void solve(TemporalGraph* G, char* query_file, char* output_file, int k);
};