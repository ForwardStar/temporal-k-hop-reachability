#include "commonfunctions.h"
#include "temporal_graph.h"
#include "online_search.h"

bool debug = false;

TemporalGraph* build(char* argv[]) {

    std::cout << "Building graph..." << std::endl;
    unsigned long long build_graph_start_time = currentTime();
    TemporalGraph* Graph = new TemporalGraph(argv[1], (char*)"Undirected");
    unsigned long long build_graph_end_time = currentTime();
    std::cout << "Build graph success in " << timeFormatting(difftime(build_graph_end_time, build_graph_start_time)).str() << std::endl;
    std::cout << "n = " << Graph->numOfVertices() << ", m = " << Graph->numOfEdges() << ", tmax = " << Graph->tmax << ", size = " << Graph->size() << " bytes" << std::endl;
    return Graph;
    
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);

    unsigned long long start_time = currentTime();

    TemporalGraph* Graph = build(argv);
    int vertex_num = Graph->numOfVertices();

    if (std::strcmp(argv[argc - 1], "Online") == 0) {
        for (int i = 2; i < argc - 2; i++) {
            std::cout << "Running online search..." << std::endl;
            unsigned long long online_search_start_time = currentTime();
            online(Graph, argv[i], argv[argc - 2]);
            unsigned long long online_search_end_time = currentTime();
            std::cout << "Online search completed in " << timeFormatting(online_search_end_time - online_search_start_time).str() << std::endl;
        }
        delete Graph;
    }

    return 0;
}