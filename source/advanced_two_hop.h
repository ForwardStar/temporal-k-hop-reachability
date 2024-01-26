#include "commonfunctions.h"
#include "temporal_graph.h"

class AdvancedTwoHopIndex {
    private:
        std::unordered_set<int> affected_vertices;
        std::vector<std::vector<std::vector<int>>> temp_paths;
        std::vector<std::vector<std::pair<int, int>>> binary_indexed_tree;
        std::queue<std::vector<int>> Q;

        void construct_for_a_vertex(TemporalGraph* G, int u, bool reverse, std::vector<std::unordered_map<int, std::vector<std::vector<std::vector<int>>>>> &L, int t_threshold);
        
    public:
        std::vector<std::unordered_map<int, std::vector<std::vector<std::vector<int>>>>> L_in, L_out;
        
        bool is_temporal_path;
        int k = 0;
        int visited_paths = 0;
        int* order;

        int size();

        bool reachable(TemporalGraph* G, int u, int v, int ts, int te, int k);

        AdvancedTwoHopIndex(TemporalGraph* G, int k_input, int t_threshold, std::string path_type);
        ~AdvancedTwoHopIndex();
        
        void solve(TemporalGraph* G, char* query_file, char* output_file);
};