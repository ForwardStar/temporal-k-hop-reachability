#include "commonfunctions.h"
#include "temporal_graph.h"
#include "online_search.h"
#include "two_hop_index.h"

bool debug = false;

TemporalGraph* build(char* argv[]) {

    std::cout << "Building graph..." << std::endl;
    unsigned long long build_graph_start_time = currentTime();
    TemporalGraph* Graph = new TemporalGraph(argv[1], (char*)"Undirected");
    unsigned long long build_graph_end_time = currentTime();
    std::cout << "Build graph success in " << timeFormatting(difftime(build_graph_end_time, build_graph_start_time)).str() << std::endl;
    std::cout << "n = " << Graph->numOfVertices() << ", m = " << Graph->numOfEdges() << ", tmax = " << Graph->tmax << std::endl;
    return Graph;
    
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);

    int k, t_threshold;
    std::cout << "Input k: ";
    std::cin >> k;
    std::cout << "Input maximum size of the query time window: ";
    std::cin >> t_threshold;

    if (std::strcmp(argv[argc - 1], "Debug") == 0) {
        debug = true;
        argc--;
    }

    std::string algorithm = "BFS-full";
    std::string algorithm_set[] = {
        std::string("PrioritySearch-naive"),
        std::string("PrioritySearch-full"),
        std::string("BFS-naive"),
        std::string("BFS-full"),
        std::string("Experimental")
    };
    if (std::strcmp(argv[argc - 2], "Index") == 0) {
        algorithm = argv[argc - 1];
        argc--;
    }
    bool in_algorithm_set = false;
    for (auto s : algorithm_set) {
        if (s == algorithm) {
            in_algorithm_set = true;
            break;
        }
    }
    if (!in_algorithm_set) {
        algorithm = "BFS-full";
    }

    unsigned long long start_time = currentTime();

    TemporalGraph* Graph = build(argv);
    int vertex_num = Graph->numOfVertices();

    if (std::strcmp(argv[argc - 1], "Online") == 0) {
        for (int i = 2; i < argc - 2; i++) {
            std::cout << "Running online search..." << std::endl;
            unsigned long long online_search_start_time = currentTime();
            online(Graph, argv[i], argv[argc - 2], k);
            unsigned long long online_search_end_time = currentTime();
            std::cout << "Online search completed in " << timeFormatting(online_search_end_time - online_search_start_time).str() << std::endl;
        }
        delete Graph;
    }

    if (std::strcmp(argv[argc - 1], "Index") == 0) {
        for (int i = 2; i < argc - 2; i++) {
            std::cout << "Running index..." << std::endl;
            std::cout << "Constructing the index structure..." << std::endl;
            unsigned long long index_construction_start_time = currentTime();
            TwoHopIndex *index = new TwoHopIndex(Graph, k, t_threshold, algorithm);
            unsigned long long index_construction_end_time = currentTime();
            std::cout << "Index construction completed in " << timeFormatting(difftime(index_construction_end_time, index_construction_start_time)).str() << std::endl;
            std::cout << "Number of intervals: " << index->size() << std::endl;
            std::cout << "Solving queries..." << std::endl;
            unsigned long long query_start_time = currentTime();
            index->solve(Graph, argv[i], argv[argc - 2], k);
            unsigned long long query_end_time = currentTime();
            std::cout << "Query completed in " << timeFormatting(query_end_time - query_start_time).str() << std::endl;
            std::cout << "Index completed!" << std::endl;
        }
        delete Graph;
    }

    return 0;
}