
#ifndef LCTD_FLATTEN_HASH_MAP_H
#define LCTD_FLATTEN_HASH_MAP_H

#include <vector>

namespace koala {
	template<typename kType, typename vType>
	class _iterator {
	
	};
	
	template<typename V>
	class my_openadd_hashmap {
#define LoadfactorInverse 1.5
#define nullele (4294967295U)
	public:
		//构造函数：初始化自定义类
		my_openadd_hashmap() : hash_table(0, std::make_pair(nullele, V())), _size(0) {}
		
		_iterator<unsigned int, V> begin() {
			return _iterator<unsigned int, V>();
		}
		
		_iterator<unsigned int, V> end() {
			return _iterator<unsigned int, V>();
		}
		
		unsigned int iterator() {
			unsigned int p = 0;
			while (p < hash_table.size() && hash_table[p].first == nullele)++p;
			return p;
		}
		
		unsigned int &next(unsigned int &p) {
			++p;
			while (p < hash_table.size() && hash_table[p].first == nullele)++p;
			return p;
		}
		
		bool has_next(unsigned int p) {
			return p != hash_table.size();
		}
		
		inline V &operator[](unsigned int key) {
			int idx = find(key);
			if (idx == -1) {
				idx = insert(key, V());
			}
			return hash_table[static_cast<unsigned int>(idx)].second;
		}
		
		inline std::pair<unsigned int, V> &get_with_idx(unsigned int i) { return hash_table[i]; }
		
		
		inline int insert(unsigned int key, V v) {
			if ((_size + 1) * LoadfactorInverse > hash_table.size()) {
				expand();
			}
			int idx = put(key, v, hash_table);
			if (idx >= 0)_size++;
			return idx;
		}
		
		inline int insert(std::pair<unsigned int, V> kv) {
			return insert(kv.first, kv.second);
		}
		
		/*
		 * if found, return idx of inner data structure
		 * if not found, return -1
		 */
		inline int find(unsigned int key) {
			if (_size == 0)
				return -1;
			unsigned int idx = key % hash_table.size();
			unsigned int tmpi = (idx + hash_table.size() - 1) % hash_table.size();
			while (hash_table[idx].first != nullele && hash_table[idx].first != key && tmpi != idx) {
				++idx;
				idx = idx == hash_table.size() ? 0 : idx;
			}
			if (hash_table[idx].first == key) {
				return static_cast<int>(idx);
			}
			return -1;
		}	
		
		inline void clear() {
			hash_table.clear(), hash_table.shrink_to_fit();
			//shrink_to_fit(): deque、vector 或 string 退回不需要的内存空间
			_size = 0;
		}
		
		inline std::vector<unsigned int> &release_data() {
			return hash_table;
		}
		
		inline unsigned int size() { return _size; }
	
	private:
		static inline int put(unsigned int key, V &v, std::vector<std::pair<unsigned int, V>> &table) {
			unsigned int idx = key % table.size();
			//std::cout<<"ini idx: ("<<idx<<", "<<key<<", "<<v<<")"<<std::endl;
			while (table[idx].first != nullele && table[idx].first != key) {
				++idx;
				idx = idx == table.size() ? 0 : idx;
			}
			//key 已 valid 存在
			//std::cout<<"after while idx: "<<idx<<std::endl;
			if (table[idx].first == nullele) {
				table[idx] = std::make_pair(key, v);
				return static_cast<int>(idx);
				//static_cast<type>(x) 把 idx 的类型强制转化成type类型 
			}
			// 是时候进行插入了

			return -1;
		}
		
		static inline int put(std::pair<unsigned int, V> p, std::vector<std::pair<unsigned int, V>> &table) {
			return put(p.first, p.second, table);
		}
		
		inline void expand() {
			unsigned int _new_size = 0;
			if (hash_table.size() == 0)_new_size = 1;
			else _new_size = hash_table.size() * 2;
			std::vector<std::pair<unsigned int, V>> newtable(_new_size, std::make_pair(nullele, V()));
			for (auto &a:hash_table)if (a.first != nullele)put(a, newtable);
			swap(newtable, hash_table);
		}
		
		std::vector<std::pair<unsigned int, V>> hash_table;
		unsigned int _size;
	};
}
#endif //LCTD_FLATTEN_HASH_MAP_H
#include<cstdio>
#include<cstring>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<set>
#include<map>
#include<queue>
#include<algorithm>
#include<ctime>
#include <sys/time.h>
using namespace std;

clock_t ct;

double GetTime(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

#define INFINITY 999999999
struct Graph{
	int n, m;
	vector<int> V;
	vector< map< int, int > > E;
	vector< vector< pair<int, int> > > Edge;
	vector<int> D;
	int *X, *Y;
	Graph(){
		n = m = 0;
		V.clear();
		E.clear();
		Edge.clear();
	}
	Graph(int tmp, char *file){
		Graph();
		FILE *fin = fopen(file, "r");
		fscanf(fin, "%d", &n);
		fscanf(fin, "%d", &m);
		for (int i = 0; i <= n; i++){
			map< int, int > v;
			v.clear();
			E.push_back(v);
		}
		for (int i = 0; i < m; i++){
			int x, y, z;
			fscanf(fin, "%d %d %d", &x, &y, &z);
			E[x].insert(make_pair(y, z));
		}
		D.clear();
		D.push_back(0);
		for (int i = 1; i <= n; i++)
			D.push_back(E[i].size());
	}
	Graph(char *file){
		Graph();
		FILE *fin = fopen(file, "r");
		fscanf(fin, "%d", &n);
		fscanf(fin, "%d", &m);
		for (int i = 0; i <= n; i++){
			map< int, int > v;
			v.clear();
			E.push_back(v);
		}
		for (int i = 0; i < m; i++){
			int x, y, z = 0;
			fscanf(fin, "%d%d%d", &x, &y, &z);
			if (E[x].find(y) != E[x].end()){
				if (E[x][y] > z){
					E[x][y] = z;
					E[y][x] = z;
				}
			}
			else{
				E[x].insert(make_pair(y, z));
				E[y].insert(make_pair(x, z));
			}
		}
		D.clear();
		D.push_back(0);
		for (int i = 1; i <= n; i++)
			D.push_back(E[i].size());
	}
	void EdgeInitialize(){
		Edge.clear();
		for (int i = 0; i <= n; i++){
			vector< pair<int, int> > Ed;
			Ed.clear();
			for (map<int, int>::iterator it = E[i].begin(); it != E[i].end(); it++){
				Ed.push_back(*it);
			}
			Edge.push_back(Ed);
		}
	}
	bool isEdgeExist(int u, int v){
		if (E[u].find(v) == E[u].end())
			return false;
		else return true;
	}
	void insertEdge(int u, int v, int k){
		if (E[u].find(v) != E[u].end()) return;
		E[u].insert(make_pair(v, k));
		E[v].insert(make_pair(u, k));
		D[u]++;
		D[v]++;
	}
	void deleteEdge(int u, int v){
		if (E[u].find(v) == E[u].end()) return;
		E[u].erase(E[u].find(v));
		E[v].erase(E[v].find(u));
		D[u]--;
		D[v]--;
	}
	void read_coordinate(char *filename){
		printf("read coordinate %s\n", filename);
		X = (int*)malloc(sizeof(int) * (n + 1));
		Y = (int*)malloc(sizeof(int) * (n + 1));
		//	printf("kkk\n");
		FILE *fco = fopen("NY.co", "r");
		//		printf("kkk\n");
		int tmp;
		fscanf(fco, "%d", &tmp);
		//	printf("tmp n: %d %d\n", tmp, n); 
		for (int i = 1; i <= n; i++){
			int v, x, y;
			fscanf(fco, "%d %d %d", &v, &x, &y);
			//		printf("v x y: %d %d %d\n", v, x, y);
			X[v] = x;
			Y[v] = y;
		}
	}
};



	struct PT{
	   	int dis;
		int x;
		PT(){
		}
		PT(int _dis, int _x){
			dis = _dis;
			x = _x;
		}
		bool operator < (const PT _pt) const{
			if (dis == _pt.dis)
				return x > _pt.x;
			return dis > _pt.dis;
		}
	};

	struct PV{
	   	int dis;
		int x;
		int k;
		bool whe_fir;
		PV(){
		}
		PV(int _dis, int _x, int _k, bool _whe_fir){
			dis = _dis;
			x = _x;
			k = _k;
			whe_fir = _whe_fir;
		}
		bool operator < (const PV _pt) const{
			if (dis == _pt.dis)
				return x > _pt.x;
			return dis > _pt.dis;
		}
	};

