#include "commonfunctions.h"
#include "temporal_graph.h"
#include "online_search.h"
#include "two_hop_index.h"
#include "advanced_two_hop.h"

bool debug = false;

TemporalGraph* build(char* argv[], std::string graph_type) {
    std::cout << "Building graph..." << std::endl;
    unsigned long long build_graph_start_time = currentTime();
    TemporalGraph* Graph = new TemporalGraph(argv[1], graph_type, 1);
    unsigned long long build_graph_end_time = currentTime();
    std::cout << "Build graph success in " << timeFormatting(difftime(build_graph_end_time, build_graph_start_time)).str() << std::endl;
    std::cout << "n = " << Graph->numOfVertices() << ", m = " << Graph->numOfEdges() << ", tmax = " << Graph->tmax << std::endl;
    return Graph;
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);

    int k, t_threshold;
    std::string sol_type, path_type, graph_type;
    std::cout << "Input the graph type (Undirected/Directed): ";
    std::cin >> graph_type;
    std::cout << "Input k: ";
    std::cin >> k;
    std::cout << "Input maximum size of the query time window: ";
    std::cin >> t_threshold;
    std::cout << "Input the solution to be used (Online/Baseline/Advanced): ";
    std::cin >> sol_type;
    std::cout << "Input the type of paths to be queried (Temporal/Projected): ";
    std::cin >> path_type;

    if (std::strcmp(argv[argc - 1], "Debug") == 0) {
        debug = true;
        argc--;
    }

    unsigned long long start_time = currentTime();

    TemporalGraph* Graph = build(argv, graph_type);
    int vertex_num = Graph->numOfVertices();

    if (sol_type == "Online") {
        for (int i = 2; i < argc - 1; i++) {
            std::cout << "Running online search..." << std::endl;
            unsigned long long online_search_start_time = currentTime();
            online(Graph, argv[i], argv[argc - 1], k, path_type);
            unsigned long long online_search_end_time = currentTime();
            std::cout << "Online search completed in " << timeFormatting(online_search_end_time - online_search_start_time).str() << std::endl;
        }
        delete Graph;
    }

    if (sol_type == "Baseline") {
        for (int i = 2; i < argc - 1; i++) {
            std::cout << "Running index..." << std::endl;
            std::cout << "Constructing the index structure..." << std::endl;
            unsigned long long index_construction_start_time = currentTime();
            TwoHopIndex *index = new TwoHopIndex(Graph, k, t_threshold, path_type);
            unsigned long long index_construction_end_time = currentTime();
            std::cout << "Index construction completed in " << timeFormatting(difftime(index_construction_end_time, index_construction_start_time)).str() << std::endl;
            std::cout << "Number of paths: " << index->size() << std::endl;
            std::cout << "Solving queries..." << std::endl;
            unsigned long long query_start_time = currentTime();
            index->solve(Graph, argv[i], argv[argc - 1], k);
            unsigned long long query_end_time = currentTime();
            std::cout << "Query completed in " << timeFormatting(query_end_time - query_start_time).str() << std::endl;
            std::cout << "Index completed!" << std::endl;
        }
        delete Graph;
    }

    if (sol_type == "Advanced") {
        for (int i = 2; i < argc - 1; i++) {
            std::cout << "Running index..." << std::endl;
            std::cout << "Constructing the index structure..." << std::endl;
            unsigned long long index_construction_start_time = currentTime();
            AdvancedTwoHopIndex *index = new AdvancedTwoHopIndex(Graph, k, t_threshold, path_type);
            unsigned long long index_construction_end_time = currentTime();
            std::cout << "Index construction completed in " << timeFormatting(difftime(index_construction_end_time, index_construction_start_time)).str() << std::endl;
            std::cout << "Number of paths: " << index->size() << std::endl;
            std::cout << "Solving queries..." << std::endl;
            unsigned long long query_start_time = currentTime();
            index->solve(Graph, argv[i], argv[argc - 1], k);
            unsigned long long query_end_time = currentTime();
            std::cout << "Query completed in " << timeFormatting(query_end_time - query_start_time).str() << std::endl;
            std::cout << "Index completed!" << std::endl;
        }
        delete Graph;
    }

    return 0;
}