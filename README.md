### KNN
This project implements the KNN algorithms for the following paper:
* Yiqi Wang, Long Yuan,  Wenjie Zhang, Xuemin Lin, Zi Chen, Qing Liu, "Simpler is More: Efficient Top-K Nearest Neighbors Search on Large Road Networks", which is submitted to VLDB2024.

### Running Environment

All experiments were conducted on a Linux machine with an Intel Xeon CPU and 384GB of memory, Debian GNU/Linux 12 and the g++ version is 12.2.0. We have implemented all methods using the C++11 standard and turned on the O3 optimization flag.

### Dataset

All the datasets in this paper can be downloaded from [DIMACS](http://www.diag.uniroma1.it/~challenge9/download.shtml)

### Compile

* compile predata.cpp for preprocessing original graph data <br>
  `g++ -std=c++11 -O3 predata.cpp -o pre`
* compile index_sdg.cpp for building bridge neighbor preserved graph for the original <br>
  `g++ -std=c++11 -O3 index_sdg.cpp -o sdg`
* compile query_up.cpp for knn-index construction, knn queries and updating objects <br>
  `g++ -std=c++11 -O3 query_up.cpp -o qu` 

### Preliminary
* There are five txt files in the 'data/test' folder: <br>
  (1) graph.txt; (2) query.txt; (3) objects.txt; (4) deletedObjects.txt; (5) insertedObjects.txt <br>
* A sample data for graph.txt, in this case, contians only 6 nodes and 6 edges, formatted as follows: <br>
  ```
  6  6     (there are total 6 vertices and 6 edges) 
  1  2  8  (an edge between vertex 1 and vertex 2 with the distance of 8)
  1  3  2
  2  4  1
  3  5  2
  4  5  3
  5  6  3
### Test
* Preprocess raw graph data from [DIMACS](http://www.diag.uniroma1.it/~challenge9/download.shtml) <br>
  `./pre oridata targetdata`   
* construct neighbor bridge graph <br>
  `./sdg data/NY.data data/NY.idx`
  
* construct knn-index and query <br>
  `./qu sdg objectset -q queries funcs topk` <br>
  `./qu data/NY.idx data/NY.object -q data/NY.query pkkvc 40`
  
* update for inserting objects into a set of candidate objects <br>
  `./qu data/NY.idx data/NY.object -u data/NY.in pkkvc 40`
  
* update for deleting objects from a set of candidate objects <br>
  `./qu data/NY.idx data/NY.object -u data/NY.de pkkvc 40`

* Arguments
  * dataset: the file path to the dataset
  * xxx: the parameter
  * 

