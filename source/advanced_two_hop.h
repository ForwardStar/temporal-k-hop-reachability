#include "commonfunctions.h"
#include "temporal_graph.h"

class AdvancedTwoHopIndex {
    private:
        std::unordered_set<int> affected_vertices;
        std::vector<std::vector<std::vector<int>>> temp_paths;
        std::vector<std::vector<int>> binary_indexed_tree;
        std::vector<std::vector<int>> temp_binary_indexed_tree;
        std::queue<std::vector<int>> Q;
        std::vector<std::pair<TemporalGraph::Edge*, bool>> deleted_edges;

        int find_index(std::vector<int> &L_in_neighbours, int u);

        void construct_for_a_vertex(TemporalGraph* G, int u, bool reverse, 
                                    std::vector<std::vector<std::vector<std::pair<int, std::vector<std::pair<int, int>>>>>> &L,
                                    std::vector<std::vector<int>> &L_neighbours,
                                    int t_threshold);
        
    public:
        std::vector<std::vector<std::vector<std::pair<int, std::vector<std::pair<int, int>>>>>> L_in, L_out;
        std::vector<std::vector<int>> L_in_neighbours, L_out_neighbours;
        
        bool is_temporal_path;
        int k = 0;
        unsigned long long visited_paths = 0;
        int* order;

        unsigned long long size();

        bool reachable(TemporalGraph* G, int u, int v, int ts, int te, int k);

        AdvancedTwoHopIndex(TemporalGraph* G, int k_input, int t_threshold, std::string path_type);
        ~AdvancedTwoHopIndex();
        
        void solve(TemporalGraph* G, char* query_file, char* output_file);
};