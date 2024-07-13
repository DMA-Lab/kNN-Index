### KNN
This project implements the KNN algorithms for the following paper:
* Yiqi Wang, Long Yuan,  Wenjie Zhang, Xuemin Lin, Zi Chen, Qing Liu, "Simpler is More: Efficient Top-K Nearest Neighbors Search on Large Road Networks"

### Dataset

All the datasets in this paper can be downloaded from [DIMACS](http://www.diag.uniroma1.it/~challenge9/download.shtml)

### Compile

* compile index_sdg.cpp for building neighbor bridge graph
  ... g++ -std=c++11 -O3 index_sdg.cpp -o sdg
* compile query_up.cpp for knn queries and updating objects
  ... g++ -std=c++11 -O3 query_up.cpp -o qu

### Preliminary

### Test


#compile index_sdg.cpp for building neighbor bridge graph

g++ -std=c++11 -O3 index_sdg.cpp -o sdg


#construct neighbor bridge graph

./sdg data/NY.data data/NY.idx

#compile query_up.cpp for knn queries and updating objects

g++ -std=c++11 -O3 query_up.cpp -o qu

#query

./qu data/NY.idx data/NY.object -q data/NY.query pkkvc optimal 40

#update for insert

./qu data/NY.idx data/NY.object -u data/NY.in pkkvc optimal 40

#update for delete

./qu data/NY.idx data/NY.object -u data/NY.de pkkvc optimal 40
