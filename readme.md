## Temporal K-hop Reachability
The repository includes four solutions solutions to solve the temporal k-hop reachability (TKR) queries: ``Online``, ``Naive``, ``MP`` and ``T2H``.

How to use it (you need to run the scripts on a Linux platform, with a Python version >= 3.5 and a gcc compiler):

- Run ``graph-gen.py`` and input the name of dataset to generate graph data automatically, which would download datasets from [SNAP](https://snap.stanford.edu/data/index.html) and [KONECT](http://konect.cc/) and process the data into ``graph.txt``:

```sh
python3 graph-gen.py
```

- Run ``query-gen.py``, input No. of queries, km and the length of query intervals to generate query data automatically, which would write queries into ``query.txt``:

```sh
python3 query-gen.py
```

- Run the following command to run the solutions, and input km and the solution to be used:

```sh
sh run.sh
```

The query results are output to the file ``output.txt``.

## Running example
We provide a running example to execute ``Online`` solution on the temporal graph in Figure 1 of the paper and set km = 2.

1. Generate graph data:
```sh
> python3 graph-gen.py
Downloading datasets...
Fetching from http://konect.cc/files/download.tsv.contact.tar.bz2...
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 106k/106k [00:01<00:00, 64.4kB/s]
Fetching from http://konect.cc/files/download.tsv.mit.tar.bz2...
100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 1.07M/1.07M [00:07<00:00, 150kB/s]
Fetching from http://konect.cc/files/download.tsv.facebook-wosn-links.tar.bz2...
100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 4.32M/4.32M [00:10<00:00, 420kB/s]
Fetching from http://konect.cc/files/download.tsv.youtube-u-growth.tar.bz2...
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 41.9M/41.9M [00:25<00:00, 1.75MB/s]
Fetching from http://konect.cc/files/download.tsv.soc-sign-bitcoinotc.tar.bz2...
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 363k/363k [00:04<00:00, 85.9kB/s]
Fetching from https://snap.stanford.edu/data/email-Eu-core-temporal.txt.gz...
100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 1.61M/1.61M [00:02<00:00, 734kB/s]
Fetching from https://snap.stanford.edu/data/CollegeMsg.txt.gz...
100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 337k/337k [00:01<00:00, 279kB/s]
Fetching from http://konect.cc/files/download.tsv.digg-friends.tar.bz2...
100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 11.7M/11.7M [00:13<00:00, 918kB/s]
Fetching from http://konect.cc/files/download.tsv.epinions.tar.bz2...
100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 3.06M/3.06M [00:09<00:00, 324kB/s]
Extracting datasets... done
Datasets:
0. example
1. facebook-wosn-links
2. contact
3. epinions
4. mit
5. email-Eu-core-temporal
6. soc-sign-bitcoinotc
7. youtube-u-growth
8. CollegeMsg
9. digg-friends
Select a graph dataset (0-9): 0
Copying dataset to "graph.txt"... done
Normalizing the graph... done
```

2. Generate random queries:
```sh
> python3 query-gen.py
Enter the number of queries to be generated: 1000
Enter kmax: 2
Enter the length of the query windows (0 < x < 1, or -1 to generate random-length intervals and random k): -1
```

3. Execute ``Online`` solution:
```sh
> sh run.sh
Compiling...
g++ source/commonfunctions.cpp -c -std=c++11 -O3
g++ source/temporal_graph.cpp -c -std=c++11 -O3
g++ source/online_search.cpp -c -std=c++11 -O3
g++ source/naive.cpp -c -std=c++11 -O3
g++ source/baseline.cpp -c -std=c++11 -O3
g++ source/advanced_two_hop.cpp -c -std=c++11 -O3
g++ source/main.cpp -c -std=c++11 -O3
g++ commonfunctions.o temporal_graph.o online_search.o naive.o baseline.o advanced_two_hop.o main.o -o main -O3
Running...
Input kmax: 2
Input the solution to be used (Online/Naive/Baseline/T2H): Online
Building graph...
Build graph success in 145μs (0s)
n = 8, m = 10, tmax = 3
Running online search...
Average: 19μs (0s)
Online search completed in 20043μs (0s)
```