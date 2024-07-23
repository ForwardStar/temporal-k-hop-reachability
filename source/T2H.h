#include "commonfunctions.h"
#include "temporal_graph.h"

class T2HIndex {
    private:
        std::unordered_set<int> affected_vertices;
        std::pair<int, int> null_interval = std::make_pair(-1, -1);

        int find_index(std::vector<int> &L_neighbours, int u);
        std::pair<int, int>& binary_search_ts_Lin(int u, int v, int k, int ts);
        std::pair<int, int>& binary_search_te_Lin(int u, int v, int k, int te);
        std::pair<int, int>& binary_search_ts_Lout(int u, int v, int k, int ts);
        std::pair<int, int>& binary_search_te_Lout(int u, int v, int k, int te);

        void construct_for_a_vertex(TemporalGraph* G, int u, bool reverse);
        
    public:
        std::vector<std::vector<std::vector<std::vector<std::pair<int, int>>>>> L_in, L_out;
        std::vector<std::vector<int>> L_in_neighbours, L_out_neighbours;

        int k = 0;
        unsigned long long visited_paths = 0;
        int* order;

        unsigned long long size();
        double max_number_of_paths();

        bool reachable(int u, int v, int ts, int te, int k);

        T2HIndex(TemporalGraph* G, int k_input);
        ~T2HIndex();
        
        void solve(TemporalGraph* G, char* query_file, char* output_file);
};