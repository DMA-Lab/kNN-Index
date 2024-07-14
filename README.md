### KNN
This project implements the KNN algorithms for the following paper:
* Yiqi Wang, Long Yuan,  Wenjie Zhang, Xuemin Lin, Zi Chen, Qing Liu, "Simpler is More: Efficient Top-K Nearest Neighbors Search on Large Road Networks", which is submitted to VLDB2024.

### Running Environment

All experiments were conducted on a Linux machine with an Intel Xeon CPU and 384GB of memory, Debian GNU/Linux 12 and the g++ version is 12.2.0. We have implemented all methods using the C++11 standard and turned on the O3 optimization flag.

### Dataset

All the datasets in this paper can be downloaded from [DIMACS](http://www.diag.uniroma1.it/~challenge9/download.shtml) 

### Preliminary
There are five files in the `data` folder: <br>
 
* `NY.data` stores the whole graph data. A sample data for `NY.data`, in this case, contians only 6 nodes and 6 edges, formatted as follows: <br>
  ```
  6  6     (there are total 6 vertices and 6 edges) 
  1  2  8  (an edge between vertex 1 and vertex 2 with the distance of 8)
  1  3  2
  2  4  1
  3  5  2
  4  5  3
  5  6  3
* `NY.qu` stores all queries. A sample data for `NY.qu` contains 3 queries, formatted as follows: <br>
  ```
  3    (there are total 3 queries)
  1 3  (a query: return top 3 nearest neighbors for query vertex 1)
  2 2
  3 4
* `NY.object` stores all obejcts. A sample data for `NY.object` contians obejcts in the candidate object set, formatted as follows: <br>
  ```
  4  (there are total 4 objects in the candidate object set)
  1  (vertex 1 belongs to the candidate object set)
  3
  5
  6
* `NY.de` stores all objected to be deleted. A sample data for `NY.de` contains objects deleted from the candidate object set, shown as follows: 
  ```
  2  (there are 2 objects, will be deleted from the set)
  3  (vertex 3 will be deleted from the candidate object set)
  5
* `NY.in` stores all objected to be inserted. A sample data for `NY.in` contains objects inserted into the candidate object set, shown as follows:
  ```
  1  (there are 1 objects, will be inserted into the set)
  2  (vertex 2 will be inserted into the candidate object set)
Note: In our experiments, all deleted objects (inserted objects) are valid. All objects to be deleted are part of the original candidate object set. While, all objects to be inserted do not belong to the original candidate object set.

### Compile

* compile predata.cpp for preprocessing the raw graph data <br>
  `g++ -std=c++11 -O3 predata.cpp -o pre`
  
* compile index_bng.cpp for building bridge neighbor preserved graph for the original graph <br>
  `g++ -std=c++11 -O3 index_bng.cpp -o bng`
  
* compile query_up.cpp for knn-index construction, knn queries and updating objects <br>
  `g++ -std=c++11 -O3 query_up.cpp -o qu`
  
### Test
* Preprocess raw graph data from [DIMACS](http://www.diag.uniroma1.it/~challenge9/download.shtml) <br>
  `./pre oridata gendata`  <br>
  * `oridata`: the file path to the raw data <br>
  * `gendata`: the file path to generated data file, which stores the formatted data <br>
  
  eg: `./pre data/NY-road-d.NY.gr data/NY.data` <br>
  
* construct bridge neighbor graph <br>
  `./pre origraph bng`  <br>
  * `origraph`: the file path to a data file, containing the original graph data <br>
  * `bng`: the file path to the data file, storing the bridge neighbor graph <br>
  
  eg: `./bng data/NY.data data/NY.idx` 
 
* construct knn-index, query and update obejcts <br>
  `./qu bng objset -X XF alg topk`
  
  * construct knn-index and query <br>
    `./qu data/NY.idx data/NY.object -q data/NY.query opt 40`
    
  * update for inserting objects into a set of candidate objects <br>
    `./qu data/NY.idx data/NY.object -u data/NY.in opt 40`
    
  * update for deleting objects from a set of candidate objects <br>
    `./qu data/NY.idx data/NY.object -u data/NY.de opt 40`

* Arguments
  * `bng`: the file path to the bridge neighbor graph
  * `objset`: the file path to the object file
  * `-X`: `-q` denotes to query, `-u` denotes to update objects
  * `XF`: the file path to queries document or the path to updates document
  * `alg`: choose a knn-index construction algo, `pri` represents our primary bottom-up computing-sharing algorithm (alg 2 in our paper), `opt` represents our optimized bidirectional construction algorihtm (alg 3 in our paper).
  * `topk`: the paramter for $k$ in the top $k$ nearesr neighbor search 

