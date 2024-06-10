#include "commonfunctions.h"
#include "temporal_graph.h"
#include "online1.h"
#include "online2.h"
#include "naive.h"
#include "MP.h"
#include "MP_optimized.h"
#include "T2H.h"

bool debug = false;

TemporalGraph* build(char* argv[]) {
    std::cout << "Building graph..." << std::endl;
    unsigned long long build_graph_start_time = currentTime();
    TemporalGraph* Graph = new TemporalGraph(argv[1], 1);
    unsigned long long build_graph_end_time = currentTime();
    std::cout << "Build graph success in " << timeFormatting(difftime(build_graph_end_time, build_graph_start_time)).str() << std::endl;
    std::cout << "n = " << Graph->numOfVertices() << ", m = " << Graph->numOfEdges() << ", tmax = " << Graph->tmax << std::endl;
    return Graph;
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    int k, t_threshold = -1;
    std::string sol_type, path_type;
    path_type = "Temporal";

    if (std::strcmp(argv[argc - 1], "Debug") == 0) {
        debug = true;
        argc--;
    }
    
    if (std::strcmp(argv[argc - 1], "Config") == 0) {
        argc--;
        k = std::atoi(argv[argc - 1]);
        argc--;
        sol_type = argv[argc - 1];
        argc--;
    }
    else {
        std::cout << "Input kmax: ";
        std::cin >> k;
        std::cout << "Input the solution to be used (Online1/Online2/Naive/MP/T2H): ";
        std::cin >> sol_type;
    }

    unsigned long long start_time = currentTime();

    TemporalGraph* Graph = build(argv);
    int vertex_num = Graph->numOfVertices();

    if (sol_type == "Online1") {
        for (int i = 2; i < argc - 1; i++) {
            std::cout << "Running online search..." << std::endl;
            unsigned long long online_search_start_time = currentTime();
            online1(Graph, argv[i], argv[argc - 1]);
            unsigned long long online_search_end_time = currentTime();
            std::cout << "Online search completed in " << timeFormatting(online_search_end_time - online_search_start_time).str() << std::endl;
        }
        delete Graph;
    }

    if (sol_type == "Online2") {
        for (int i = 2; i < argc - 1; i++) {
            std::cout << "Running online search..." << std::endl;
            unsigned long long online_search_start_time = currentTime();
            online2(Graph, argv[i], argv[argc - 1]);
            unsigned long long online_search_end_time = currentTime();
            std::cout << "Online search completed in " << timeFormatting(online_search_end_time - online_search_start_time).str() << std::endl;
        }
        delete Graph;
    }

    if (sol_type == "Naive") {
        std::cout << "Running index..." << std::endl;
        std::cout << "Constructing the index structure..." << std::endl;
        unsigned long long index_construction_start_time = currentTime();
        NaiveIndex *index = new NaiveIndex(Graph, k, t_threshold, path_type);
        unsigned long long index_construction_end_time = currentTime();
        std::cout << "Index construction completed in " << timeFormatting(difftime(index_construction_end_time, index_construction_start_time)).str() << std::endl;
        std::cout << "Number of paths in index: " << index->size() << std::endl;
        for (int i = 2; i < argc - 1; i++) {
            std::cout << "Solving queries..." << std::endl;
            unsigned long long query_start_time = currentTime();
            index->solve(Graph, argv[i], argv[argc - 1]);
            unsigned long long query_end_time = currentTime();
            std::cout << "Query completed in " << timeFormatting(query_end_time - query_start_time).str() << std::endl;
            std::cout << "Index completed!" << std::endl;
        }
        delete Graph;
    }

    if (sol_type == "MP") {
        std::cout << "Running index..." << std::endl;
        std::cout << "Constructing the index structure..." << std::endl;
        unsigned long long index_construction_start_time = currentTime();
        MPIndex *index = new MPIndex(Graph, k);
        unsigned long long index_construction_end_time = currentTime();
        std::cout << "Index construction completed in " << timeFormatting(difftime(index_construction_end_time, index_construction_start_time)).str() << std::endl;
        std::cout << "Number of paths in index: " << index->size() << std::endl;
        std::cout << "Alpha: " << index->alpha << std::endl;
        for (int i = 2; i < argc - 1; i++) {
            std::cout << "Solving queries..." << std::endl;
            unsigned long long query_start_time = currentTime();
            index->solve(Graph, argv[i], argv[argc - 1]);
            unsigned long long query_end_time = currentTime();
            std::cout << "Query completed in " << timeFormatting(query_end_time - query_start_time).str() << std::endl;
            std::cout << "Index completed!" << std::endl;
        }
        delete Graph;
    }

    if (sol_type == "MPO") {
        std::cout << "Running index..." << std::endl;
        std::cout << "Constructing the index structure..." << std::endl;
        unsigned long long index_construction_start_time = currentTime();
        MPIndexO *index = new MPIndexO(Graph, k);
        unsigned long long index_construction_end_time = currentTime();
        std::cout << "Index construction completed in " << timeFormatting(difftime(index_construction_end_time, index_construction_start_time)).str() << std::endl;
        std::cout << "Number of paths in index: " << index->size() << std::endl;
        std::cout << "Alpha: " << index->alpha << std::endl;
        for (int i = 2; i < argc - 1; i++) {
            std::cout << "Solving queries..." << std::endl;
            unsigned long long query_start_time = currentTime();
            index->solve(Graph, argv[i], argv[argc - 1]);
            unsigned long long query_end_time = currentTime();
            std::cout << "Query completed in " << timeFormatting(query_end_time - query_start_time).str() << std::endl;
            std::cout << "Index completed!" << std::endl;
        }
        delete Graph;
    }

    if (sol_type == "T2H") {
        std::cout << "Running index..." << std::endl;
        std::cout << "Constructing the index structure..." << std::endl;
        unsigned long long index_construction_start_time = currentTime();
        T2HIndex *index = new T2HIndex(Graph, k);
        unsigned long long index_construction_end_time = currentTime();
        std::cout << "Index construction completed in " << timeFormatting(difftime(index_construction_end_time, index_construction_start_time)).str() << std::endl;
        std::cout << "Number of paths in index: " << index->size() << std::endl;
        std::cout << "Beta: " << index->max_number_of_paths() << std::endl;
        for (int i = 2; i < argc - 1; i++) {
            std::cout << "Solving queries..." << std::endl;
            unsigned long long query_start_time = currentTime();
            index->solve(Graph, argv[i], argv[argc - 1]);
            unsigned long long query_end_time = currentTime();
            std::cout << "Query completed in " << timeFormatting(query_end_time - query_start_time).str() << std::endl;
            std::cout << "Index completed!" << std::endl;
        }
        delete Graph;
    }

    return 0;
}