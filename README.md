# disjoint-k-clique

The first line of each file is how to compile the file.

## Usage
### TestCliqueNum.cpp
get the number of k-cliques
```
./TestCliqueNum yourGraphFilePath
```

### TestCliqueDensity.cpp
calculate alpha
```
./TestCliqueDensity yourGraphFilePath
```

### greedyOriginalDegreeDAG.cpp
Algo. 2 in the paper
```
./greedyOriginalDegreeDAG yourGraphFilePath k
```

### greedyAllclique.cpp
Algo. 3 in the paper
```
./greedyAllclique yourGraphFilePath k
```

### thirdGreedyOnceMultiCore.cpp
Algo. 4 in the paper without the pruning strategy
```
./thirdGreedyOnceMultiCore yourGraphFilePath k
```

### thirdGreedyOnceMultiCorePrune.cpp
Algo. 4 in the paper with the pruning strategy
```
./thirdGreedyOnceMultiCorePrune yourGraphFilePath k
```

### bench.cpp
An example of how to conduct update operations. The important operations are listed as follows.
```c++
DynamicGraph dg(yourGraphFilePath, k) // init

dg.Insertion(u, v) // add an edge

dg.Deletion(u, v) // delete an edge
```