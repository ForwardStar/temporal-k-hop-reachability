#include "commonfunctions.h"
#include "temporal_graph.h"

class T2HIndex {
    private:
        std::vector<int> f;
        std::unordered_set<int> visited_vertices;
        std::pair<int, int> null_interval = std::make_pair(-1, -1);
        std::pair<int, int>& binary_search_ts_Lin(int u, int v, int k, int ts);
        std::pair<int, int>& binary_search_te_Lin(int u, int v, int k, int te);
        std::pair<int, int>& binary_search_ts_Lout(int u, int v, int k, int ts);
        std::pair<int, int>& binary_search_te_Lout(int u, int v, int k, int te);

        void construct_for_a_vertex(TemporalGraph* G, int u, bool reverse, std::vector<std::unordered_map<int, std::vector<std::vector<std::pair<int, int>>>>> &L);
        
    public:
        std::vector<std::unordered_map<int, std::vector<std::vector<std::pair<int, int>>>>> L_in, L_out;

        int k = 0;
        unsigned long long visited_paths = 0;
        int* order;

        unsigned long long size();
        unsigned long long max_number_of_paths();

        bool reachable(int u, int v, int ts, int te, int k);

        T2HIndex(TemporalGraph* G, int k_input);
        ~T2HIndex();
        
        void solve(TemporalGraph* G, char* query_file, char* output_file);
};