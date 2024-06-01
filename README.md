g++ -std=c++11 -O3 index_sdg.cpp -o sdg

# construct neighbor bridge graph

./sdg data/NY.data data/NY.idx



g++ -std=c++11 -O3 query_up.cpp -o qu

#query
./qu data/NY.idx data/NY.object -q data/NY.query pkkvc optimal 40

#update for insert
./qu data/NY.idx data/NY.object -u data/NY.in pkkvc optimal 40

#update for delete
./qu data/NY.idx data/NY.object -u data/NY.de pkkvc optimal 40
