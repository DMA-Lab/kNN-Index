#include<cstdio>
#include<cstring>
#include<iostream>
#include<cstdlib>
#include<queue>
#include<sys/time.h>
#include<vector>
#include<xmmintrin.h>
#include<cmath>
#include<set>
#include<algorithm>
#include<fstream>
#include<map>
#include "dynamic_graph.h"
#include "flatten_hash_map.h"
using namespace std;
//#define MAX_K 20
#define INT_MAX 999999999
#define RESERVE_TIME 1

int reset_times = 0;

const int SIZEOFINT = 4;

double cnt_pre_query_time = 0;
char * insert_type, * query_type;


int *height, *pa, *uniqueVertex;
int *belong;
int root, TreeSize;
int **rootToRoot, *rootSite;
int **dis, **pos, **gdis;
int *posSize;
int *chSize;
int ** ch;//, ** kch;
vector<vector<int>> kch;
int *LOG2, *LOGD; 
int rootSize;
int *DFSList, *toDFS;
int ***BS;

int **pv;
int *pvSize;

int **pvd;
int *pvdSize;

int **kpv;
int *kpvSize;

int **kpvd;
int *kpvdSize;

int **kfv;
int *kfvSize;

int **kfvd;
int *kfvdSize;
int *kpa;

vector<map<int,int>> v_d;
//vector<hash_map<int,int>> v_d_hm;
//koala::my_openadd_hashmap<unsigned int> edges;




struct HASH_NODE{
    int a, b, c;
    int next;
    HASH_NODE(){}
    HASH_NODE(int _a, int _b, int _c){
        a = _a;
        b = _b;
        c = _c;
    }
};
class HASH{
public:
    static const int P = 100000007;
    vector<HASH_NODE> nodes;
    vector<int> start;
    HASH(){
        nodes.clear();
        nodes.push_back(HASH_NODE());
        start.resize(P);
        for (int i = 0; i < P; i++)
            start[i] = 0;
    }
    inline int get_id(int a, int b){
        return ((long long )(a) * b) % P;
    }
    inline bool is_exist(int a, int b){
        int t = get_id(a, b);
        int p = start[t];
        while (p != 0){
            if (nodes[p].a == a && nodes[p].b == b)
                return true;
            p = nodes[p].next;
        }
        return false;
    }
    inline int get_value(int a, int b){
        int t = get_id(a, b);
        int p = start[t];
        while (p != 0){
            if (nodes[p].a == a && nodes[p].b == b)
                return nodes[p].c;
            p = nodes[p].next;
        }
        return -1;
    }
    inline bool insert_node(int a, int b, int c){
        if (is_exist(a, b) == true)
            return false;
        int t = get_id(a, b);
        HASH_NODE hn = HASH_NODE(a,b,c);
        hn.next = start[t];
        nodes.push_back(hn);
        start[t] = nodes.size() - 1;
        return true;
    }
    inline bool delete_node(int a, int b){
        int t = get_id(a, b);
        int p = start[t];
        int pre = -1;
        while (p != 0){
            if (nodes[p].a == a && nodes[p].b == b){
                if (pre < 0)
                    start[t] = nodes[p].next;
                else nodes[pre].next = nodes[p].next;
                return true;
            }
            pre = p;
            p = nodes[p].next;
        }
        return false;
    }
};

struct DN{
	int* dis;
	int x;
	DN(){}
	DN(int* _dis, int _x){
		dis = _dis;
		x = _x;
	}
	//根据dis 指向的地址的数值排序
	bool operator < (const DN _pt) const{
		if (*dis == *_pt.dis) return x < _pt.x;
		return *dis < *_pt.dis;
	}
		//~DN()
};

class prio_queue{
	private:
	    deque<DN> data;
	public:
	    void push(DN pt){
			//cout<<"push"<<endl;
			data.push_back(pt);
			//push_heap(data.begin(),data.end()); //O(logn)
		}
		void pop(){
			//cout<<"pop"<<endl;
			//pop_heap(data.begin(), data.end());//delete
			data.pop_front();
		}
		DN top(){
			//cout<<"top"<<endl;
			sort(data.begin(),data.end());//pop_heap(data.begin(), data.end()); //添加+pop_heap
			return data.front();
		}
		bool empty(){
			//cout<<"empty"<<endl;
			return data.empty();
		}
};

	long long queryCnt;	
	//long long aCnt;

	FILE *fin;
	string TT = "";
	void scanIntArray(int *a, int n){
		fread(a, SIZEOFINT, n, fin);
	}
	int* scanIntVector(int *a){
		int _n;
		fread(&_n, SIZEOFINT, 1, fin);
		a = (int*)malloc(sizeof(int) * _n);
		scanIntArray(a, _n);
		return a;
	}

	int n;
	int tree_height = 0;
	void readIndex(char *file){
		//cout<<"com read index"<<endl;
		double _time = GetTime();
		int tree_width = 0, most_sp = 0;
		fin = fopen(file, "rb");
		fread(&n, SIZEOFINT, 1, fin);
		//printf("n: %d\n", n);
		int ts;
		fread(&ts, SIZEOFINT, 1, fin);
		TreeSize = ts;
		height = (int*)malloc(sizeof(int) * (ts + 1));
		pa = (int*)malloc(sizeof(int) * (ts + 1));
		kpa = (int*)malloc(sizeof(int) * (ts + 1));
		uniqueVertex = (int*)malloc(sizeof(int) * (ts + 1));
		int print = 0;
		int pit = 0;
		//if(pit) cout<<"height: ";
		for (int i = 0; i < ts; i++){
			fread(&height[i], SIZEOFINT, 1, fin);
			//if(pit) cout<<height[i]<<" ";
			if (height[i] > tree_height)
				tree_height = height[i];
		}
		//if(pit) cout<<endl;
		if(pit) cout<<"pa: "<<" ";
		for (int i = 0; i < ts; i++){
			fread(&pa[i], SIZEOFINT, 1, fin);
			if(pit) cout<<pa[i]<<" ";
		}
		//pa[0] =-1;
		if (pit) cout<<endl;
		if (pit) cout<<"uV: ";
		for (int i = 0; i < ts; i++){
			fread(&uniqueVertex[i], SIZEOFINT, 1, fin);
			if (pit) cout<<uniqueVertex[i]<<" ";
		}
		if (pit) cout<<endl;
		belong = (int*)malloc(sizeof(int) * (n + 1));
	  	fread(belong, SIZEOFINT, n + 1, fin);
		fread(&root, SIZEOFINT, 1, fin);

		if (pit) cout<<"belong: ";
		for (int i = 1; i < ts+1; i++){
			if (pit) cout<<belong[i]<<" ";
		}
		if (pit) cout<<endl;

		if (pit) cout<<"root: "<<root<<endl;

		posSize = (int*)malloc(sizeof(int) * (n + 1));
		pos = (int**)malloc(sizeof(int*) * (TreeSize));

		pvSize = (int*)malloc(sizeof(int) * (n + 1));
		pv = (int**)malloc(sizeof(int*) * (TreeSize));

		pvdSize = (int*)malloc(sizeof(int) * (n + 1));
		pvd = (int**)malloc(sizeof(int*) * (TreeSize));

		kpvSize = (int*)malloc(sizeof(int) * (n + 1));
		kpv = (int**)malloc(sizeof(int*) * (TreeSize));

		kpvdSize = (int*)malloc(sizeof(int) * (n + 1));
		kpvd = (int**)malloc(sizeof(int*) * (TreeSize));

		kfvSize = (int*)malloc(sizeof(int) * (n + 1));
		kfv = (int**)malloc(sizeof(int*) * (TreeSize));

		kfvdSize = (int*)malloc(sizeof(int) * (n + 1));
		kfvd = (int**)malloc(sizeof(int*) * (TreeSize));

		dis = (int**)malloc(sizeof(int*) * (TreeSize));
		chSize = (int*)malloc(sizeof(int) * (TreeSize));
		ch = (int**)malloc(sizeof(int*) * (TreeSize));

		kch.resize(TreeSize);

		v_d.clear();
        map<int,int> vemp;
	    vemp.clear();
	    for (int i = 0; i <= n; i++){
		    v_d.push_back(vemp);
        }
		for (int i = 0; i < TreeSize; i++){
			fread(&chSize[i], SIZEOFINT, 1, fin);
			ch[i] = (int*)malloc(sizeof(int) * chSize[i]);

			for (int j = 0; j < chSize[i]; j++){
				int x;
				fread(&x, SIZEOFINT, 1, fin);
				ch[i][j] = x;
			}
		}
		for (int i = 0; i < TreeSize; i++){
			// read vert[]
			int x =0;
			fread(&x, SIZEOFINT, 1, fin);
			//cout<<"x: "<<x<<" uniqueVertex[x]: "<<uniqueVertex[x]<<endl;
			fread(&posSize[x], SIZEOFINT, 1, fin);
			if (posSize[x] > tree_width) tree_width = posSize[x];
		
			//posSize[x]
			pos[x] = (int*)malloc(sizeof(int) * (posSize[x] + 1));
			fread(pos[x], SIZEOFINT, posSize[x], fin);
			//cout<<"(pos, belong): ";
			//for(int j =0;j<posSize[x];j++) cout<<" ("<<pos[x][j]<<", "<<belong[pos[x][j]]<<") ";
			//cout<<endl;
			//cout<<"uni[x]: "<< uniqueVertex[x]<<" pa: "<<uniqueVertex[pa[x]]<<" belong[pa]: "<< pa[x]<<endl;
			
			int _n;
			fread(&_n, SIZEOFINT, 1, fin);
			dis[x] = (int*)malloc(sizeof(int) * _n);
			fread(dis[x], SIZEOFINT, _n, fin);

			fread(&pvSize[x], SIZEOFINT, 1, fin);
			pv[x] = (int*)malloc(sizeof(int) * pvSize[x]);
			fread(pv[x], SIZEOFINT, pvSize[x], fin);

			fread(&pvdSize[x], SIZEOFINT, 1, fin);
			pvd[x] = (int*)malloc(sizeof(int) * pvdSize[x]);
			fread(pvd[x], SIZEOFINT, pvdSize[x], fin);

			fread(&kpvSize[x], SIZEOFINT, 1, fin);
			kpv[x] = (int*)malloc(sizeof(int) * kpvSize[x]);
			fread(kpv[x], SIZEOFINT, kpvSize[x], fin);
         
			fread(&kpvdSize[x], SIZEOFINT, 1, fin);
			kpvd[x] = (int*)malloc(sizeof(int) * kpvdSize[x]);
			fread(kpvd[x], SIZEOFINT, kpvdSize[x], fin);

			fread(&kfvSize[x], SIZEOFINT, 1, fin);
			kfv[x] = (int*)malloc(sizeof(int) * kfvSize[x]);
			fread(kfv[x], SIZEOFINT, kfvSize[x], fin);

			fread(&kfvdSize[x], SIZEOFINT, 1, fin);
			kfvd[x] = (int*)malloc(sizeof(int) * kfvdSize[x]);
			fread(kfvd[x], SIZEOFINT, kfvdSize[x], fin);
            
			kpa[x] = -1;
			for (int i=0;i<kfvdSize[x];i++){
				if(belong[kfv[x][i]]> kpa[x]) kpa[x] = belong[kfv[x][i]];
				kch[kpa[x]].push_back(x);
			}

			for(int j=0; j<_n; j++){
				v_d[x].insert(make_pair(pos[x][j],dis[x][j]));
			}

			if(print) cout<<"x: "<< uniqueVertex[x] <<" pv[x][j]: ";
			if(print) for(int j=0; j<pvSize[x]; j++) cout<<" ("<<pv[x][j]<<","<<pvd[x][j]<<") ";
			if(print) cout<<endl;
			
			if(print) cout<<"x: "<< uniqueVertex[x] <<"kpv[x][j]: ";
			if(print) for(int j=0; j<kpvSize[x]; j++) cout<<" ("<<kpv[x][j]<<","<<kpvd[x][j]<<") ";
			if(print) cout<<endl;

			if(print) cout<<"x: "<< uniqueVertex[x] <<"fv[x][j]: ";
			if(print) for(int j=0; j<posSize[x]; j++) cout<<" ("<<pos[x][j]<<","<<dis[x][j]<<") ";
			if(print) cout<<endl;

			if(print) cout<<"x: "<< uniqueVertex[x] <<"kfv[x][j]: ";
			if(print) for(int j=0; j<kfvSize[x]; j++) cout<<" ("<<kfv[x][j]<<","<<kfvd[x][j]<<") ";
			if(print) cout<<endl;
		}
		//for(int i = 0; i < TreeSize; i++) cout<<" ("<<uniqueVertex[i]<<", "<<uniqueVertex[kpa[i]]<<") ";
		//cout<<endl;
		
		fclose(fin);
		//printf("Load Index Time : %lf sec\n", (GetTime() - _time));
		//printf("tree height: %d\n", tree_height);
		//printf("tree width: %d\n", tree_width);
		//printf("most search space: %d\n", most_sp);
	}

	void readIndex_up(char *file){
		//cout<<"com read index"<<endl;
		double _time = GetTime();
		int tree_width = 0, most_sp = 0;
		fin = fopen(file, "rb");
		fread(&n, SIZEOFINT, 1, fin);
		//printf("n: %d\n", n);
		int ts;
		fread(&ts, SIZEOFINT, 1, fin);
		TreeSize = ts;
		height = (int*)malloc(sizeof(int) * (ts + 1));
		pa = (int*)malloc(sizeof(int) * (ts + 1));
		uniqueVertex = (int*)malloc(sizeof(int) * (ts + 1));
		int print = 0;
		if(print) cout<<"height: ";
		for (int i = 0; i < ts; i++){
			fread(&height[i], SIZEOFINT, 1, fin);
			if(print) cout<<height[i]<<" ";
			if (height[i] > tree_height)
				tree_height = height[i];
		}
		if(print) cout<<endl;
		for (int i = 0; i < ts; i++){
			fread(&pa[i], SIZEOFINT, 1, fin);
		}
		for (int i = 0; i < ts; i++){
			fread(&uniqueVertex[i], SIZEOFINT, 1, fin);
		}
		belong = (int*)malloc(sizeof(int) * (n + 1));
	  	fread(belong, SIZEOFINT, n + 1, fin);
		fread(&root, SIZEOFINT, 1, fin);
		posSize = (int*)malloc(sizeof(int) * (n + 1));
		pos = (int**)malloc(sizeof(int*) * (TreeSize));

		pvSize = (int*)malloc(sizeof(int) * (n + 1));
		pv = (int**)malloc(sizeof(int*) * (TreeSize));

		pvdSize = (int*)malloc(sizeof(int) * (n + 1));
		pvd = (int**)malloc(sizeof(int*) * (TreeSize));

		kpvSize = (int*)malloc(sizeof(int) * (n + 1));
		kpv = (int**)malloc(sizeof(int*) * (TreeSize));

		kpvdSize = (int*)malloc(sizeof(int) * (n + 1));
		kpvd = (int**)malloc(sizeof(int*) * (TreeSize));

		kfvSize = (int*)malloc(sizeof(int) * (n + 1));
		kfv = (int**)malloc(sizeof(int*) * (TreeSize));

		kfvdSize = (int*)malloc(sizeof(int) * (n + 1));
		kfvd = (int**)malloc(sizeof(int*) * (TreeSize));

		dis = (int**)malloc(sizeof(int*) * (TreeSize));
		gdis = (int**)malloc(sizeof(int*) * (TreeSize));
		chSize = (int*)malloc(sizeof(int) * (TreeSize));
		ch = (int**)malloc(sizeof(int*) * (TreeSize));

		v_d.clear();
        map<int,int> vemp;
	    vemp.clear();
	    for (int i = 0; i <= n; i++){
		    v_d.push_back(vemp);
        }
		for (int i = 0; i < TreeSize; i++){
			fread(&chSize[i], SIZEOFINT, 1, fin);
			ch[i] = (int*)malloc(sizeof(int) * chSize[i]);
			for (int j = 0; j < chSize[i]; j++){
				int x;
				fread(&x, SIZEOFINT, 1, fin);
				ch[i][j] = x;
			}
		}

		for (int i = 0; i < TreeSize; i++){
			// read vert[]
			int x =0;
			fread(&x, SIZEOFINT, 1, fin);
			//cout<<"x: "<<x<<" uniqueVertex[x]: "<<uniqueVertex[x]<<endl;
			fread(&posSize[x], SIZEOFINT, 1, fin);
			if (posSize[x] > tree_width) tree_width = posSize[x];
			// (posSize[x] + 1)
			//posSize[x]
			pos[x] = (int*)malloc(sizeof(int) * (posSize[x] + 1));
			fread(pos[x], SIZEOFINT, posSize[x], fin);
			//read VL[]
			int _n;
			fread(&_n, SIZEOFINT, 1, fin);
			dis[x] = (int*)malloc(sizeof(int) * _n);
			fread(dis[x], SIZEOFINT, _n, fin);

			fread(&_n, SIZEOFINT, 1, fin);
			gdis[x] = (int*)malloc(sizeof(int) * _n);
			fread(gdis[x], SIZEOFINT, _n, fin);

			fread(&pvSize[x], SIZEOFINT, 1, fin);
			pv[x] = (int*)malloc(sizeof(int) * pvSize[x]);
			fread(pv[x], SIZEOFINT, pvSize[x], fin);

			fread(&pvdSize[x], SIZEOFINT, 1, fin);
			pvd[x] = (int*)malloc(sizeof(int) * pvdSize[x]);
			fread(pvd[x], SIZEOFINT, pvdSize[x], fin);

			fread(&kpvSize[x], SIZEOFINT, 1, fin);
			kpv[x] = (int*)malloc(sizeof(int) * kpvSize[x]);
			fread(kpv[x], SIZEOFINT, kpvSize[x], fin);
         
			fread(&kpvdSize[x], SIZEOFINT, 1, fin);
			kpvd[x] = (int*)malloc(sizeof(int) * kpvdSize[x]);
			fread(kpvd[x], SIZEOFINT, kpvdSize[x], fin);

			fread(&kfvSize[x], SIZEOFINT, 1, fin);
			kfv[x] = (int*)malloc(sizeof(int) * kfvSize[x]);
			fread(kfv[x], SIZEOFINT, kfvSize[x], fin);

			fread(&kfvdSize[x], SIZEOFINT, 1, fin);
			kfvd[x] = (int*)malloc(sizeof(int) * kfvdSize[x]);
			fread(kfvd[x], SIZEOFINT, kfvdSize[x], fin);

			//read v_d 
			//v_d[i] = (map<int,int>*)malloc(2*sizeof(int) * (_n));
			for(int j=0; j<_n; j++){
				v_d[x].insert(make_pair(pos[x][j],gdis[x][j]));
				//v_d_hm[x].insert(make_pair(pos[x][j],gdis[x][j]));
			}

			if(print) cout<<"x: "<< uniqueVertex[x] <<" pv[x][j]: ";
			if(print) for(int j=0; j<pvSize[x]; j++) cout<<" ("<<pv[x][j]<<","<<pvd[x][j]<<") ";
			if(print) cout<<endl;
			
			if(print) cout<<"x: "<< uniqueVertex[x] <<"kpv[x][j]: ";
			if(print) for(int j=0; j<kpvSize[x]; j++) cout<<" ("<<kpv[x][j]<<","<<kpvd[x][j]<<") ";
			if(print) cout<<endl;

			if(print) cout<<"x: "<< uniqueVertex[x] <<"fv[x][j]: ";
			if(print) for(int j=0; j<posSize[x]; j++) cout<<" ("<<pos[x][j]<<","<<dis[x][j]<<") ";
			if(print) cout<<endl;

			if(print) cout<<"x: "<< uniqueVertex[x] <<"kfv[x][j]: ";
			if(print) for(int j=0; j<kfvSize[x]; j++) cout<<" ("<<kfv[x][j]<<","<<kfvd[x][j]<<") ";
			if(print) cout<<endl;
		}
		fclose(fin);
		//printf("Load Index Time : %lf sec\n", (GetTime() - _time));
		//printf("tree height: %d\n", tree_height);
		//printf("tree width: %d\n", tree_width);
		//printf("most search space: %d\n", most_sp);
	}

	void readIndex_ori(char *file){
		//cout<<"com read index"<<endl;
		double _time = GetTime();
		int tree_height = 0, tree_width = 0, most_sp = 0;
		fin = fopen(file, "rb");
		fread(&n, SIZEOFINT, 1, fin);
		//printf("n: %d\n", n);
		int ts;
		fread(&ts, SIZEOFINT, 1, fin);
		TreeSize = ts;
		height = (int*)malloc(sizeof(int) * (ts + 1));
		pa = (int*)malloc(sizeof(int) * (ts + 1));
		uniqueVertex = (int*)malloc(sizeof(int) * (ts + 1));
		for (int i = 0; i < ts; i++){
			fread(&height[i], SIZEOFINT, 1, fin);
			if (height[i] > tree_height)
				tree_height = height[i];
		}
		for (int i = 0; i < ts; i++){
			fread(&pa[i], SIZEOFINT, 1, fin);
		}
		for (int i = 0; i < ts; i++){
			fread(&uniqueVertex[i], SIZEOFINT, 1, fin);
		}
		//cout<<"111"<<endl;
		belong = (int*)malloc(sizeof(int) * (n + 1));
	  	fread(belong, SIZEOFINT, n + 1, fin);
		//cout<<"111"<<endl;
		fread(&root, SIZEOFINT, 1, fin);
		cout << "root: " << root << endl;
		posSize = (int*)malloc(sizeof(int) * (n + 1));
		pos = (int**)malloc(sizeof(int*) * (TreeSize));
		dis = (int**)malloc(sizeof(int*) * (TreeSize));
		chSize = (int*)malloc(sizeof(int) * (TreeSize));
		ch = (int**)malloc(sizeof(int*) * (TreeSize));
		//初始化！！
		//v_d = (map<int,int>*)malloc(2*sizeof(int) * (TreeSize));
		v_d.clear();
        map<int,int> vemp;
	    vemp.clear();
	    for (int i = 0; i <= n; i++){
		    v_d.push_back(vemp);
        }
		for (int i = 0; i < TreeSize; i++){
			fread(&chSize[i], SIZEOFINT, 1, fin);
			ch[i] = (int*)malloc(sizeof(int) * chSize[i]);
			for (int j = 0; j < chSize[i]; j++){
				int x;
				fread(&x, SIZEOFINT, 1, fin);
				ch[i][j] = x;
			}
		}

		for (int i = 0; i < TreeSize; i++){
			// read vert[]
			int x =0;
			fread(&x, SIZEOFINT, 1, fin);
			//cout<<"x: "<<x<<" uniqueVertex[x]: "<<uniqueVertex[x]<<endl;
			fread(&posSize[x], SIZEOFINT, 1, fin);
			if (posSize[x] > tree_width)
				tree_width = posSize[x];
			pos[x] = (int*)malloc(sizeof(int) * (posSize[x] + 1));
			fread(pos[x], SIZEOFINT, posSize[x], fin);
			//read VL[]
			int _n;
			fread(&_n, SIZEOFINT, 1, fin);
			dis[x] = (int*)malloc(sizeof(int) * _n);
			fread(dis[x], SIZEOFINT, _n, fin);
			//read v_d 
			//v_d[i] = (map<int,int>*)malloc(2*sizeof(int) * (_n));
			for(int j=0; j<_n; j++){
				v_d[x].insert(make_pair(pos[x][j],dis[x][j]));
				//cout<<"pos[x][j]: "<<pos[x][j]<<endl;
				//cout<<"dis[x][j]: "<<dis[x][j]<<endl;
			}
		}
		fclose(fin);
		
		printf("Load Index Time : %lf sec\n", (GetTime() - _time));
		//printf("tree height: %d\n", tree_height);
		//printf("tree width: %d\n", tree_width);
		//printf("most search space: %d\n", most_sp);
		
	}
	


int *is_current_object, *current_distance, current_stamp, *current_state;

int *curup_dist, curup_stamp, *curup_affect, *curup_fvSize;


int *ug_affect, ug_stamp, *uVToUg_ID; 

int djik_stamp, *djik_state, *distance_computed;



class kNN{
public:	
int M_K;
//FILE *fout;
kNN(int x){
	M_K = x;
}
//双向链表
struct list_node{
	int previous, next, key, dist;
	list_node(){
		previous = -1;
		next = -1;
		key = -1;         // vertex
		dist = -1;        // distance
	}
};

struct object_saveing_structure{
	vector<list_node> a;
	vector<int> trash_can;
	int current, size_num; 
	object_saveing_structure(){
		a.clear();
		list_node _a;
		_a.key = -1;
		_a.dist = -1;
		_a.previous = 0;
		_a.next = 0;
		a.push_back(_a);
		trash_can.clear();
		current = 0;
		size_num = 0;
	}
};

    HASH hash;
    vector<double> times_period[10];
    int period = 8;
	vector<object_saveing_structure> OSS;
 	vector<object_saveing_structure> OSS_global;
    vector<int> object_number;
	

	//initailize knn structure
	void create_kNN_index(){

		//initialize OSS
		OSS.clear();
		for (int i = 0; i < TreeSize; i++){
			object_saveing_structure oss;
			OSS.push_back(oss);
		}

		//initialize OSS_global
		OSS_global.clear();
		for (int i = 0; i < TreeSize; i++){
			object_saveing_structure oss;
			OSS_global.push_back(oss);
		}
        
	}
	//对knn 进行删除object的操作
	void delete_element(int p, int x){
		OSS[p].size_num--;
		int pre = OSS[p].a[x].previous;
		int ne = OSS[p].a[x].next;
		OSS[p].a[ne].previous = pre;
		OSS[p].a[pre].next = ne;
	}
	void delete_element_global(int p, int x){
		OSS_global[p].size_num--;
		int pre = OSS_global[p].a[x].previous;
		int ne = OSS_global[p].a[x].next;
		OSS_global[p].a[ne].previous = pre;
		OSS_global[p].a[pre].next = ne;
	}
	//OSS OSS_GLOBAL push
	static bool object_compare(int a, int b){
		if (current_distance[a] < current_distance[b])
			return true;
		else return false;
	}

	void OSS_push_back(int p, int key, int dist){
		//	printf("p, key, dist: %d %d %d\n", p, key, dist);
		list_node _a;
		_a.previous = OSS[p].a[0].previous;
		_a.next = 0;
		_a.key = key;
		_a.dist = dist;
		OSS[p].a.push_back(_a);
		OSS[p].a[_a.previous].next = OSS[p].a.size() - 1;
		OSS[p].a[_a.next].previous = OSS[p].a.size() - 1;
		OSS[p].size_num++;
	}
	void OSS_push_front(int p, int key, int dist){
		list_node _a;
		_a.previous = 0;
		_a.next = OSS[p].a[0].next;
		_a.key = key;
		_a.dist = dist;
		OSS[p].a.push_back(_a);
		OSS[p].a[_a.previous].next = OSS[p].a.size() - 1;
		OSS[p].a[_a.next].previous = OSS[p].a.size() - 1;
		OSS[p].size_num++;
	}
	void OSS_global_push_back(int p, int key, int dist){
		//	printf("p, key, dist: %d %d %d\n", p, key, dist);
		list_node _a;
		_a.previous = OSS_global[p].a[0].previous;
		_a.next = 0;
		_a.key = key;
		_a.dist = dist;
		OSS_global[p].a.push_back(_a);
		OSS_global[p].a[_a.previous].next = OSS_global[p].a.size() - 1;
		OSS_global[p].a[_a.next].previous = OSS_global[p].a.size() - 1;
		OSS_global[p].size_num++;
	}
	void OSS_global_push_front(int p, int key, int dist){
		list_node _a;
		_a.previous = 0;
		_a.next = OSS_global[p].a[0].next;
		_a.key = key;
		_a.dist = dist;
		OSS_global[p].a.push_back(_a);
		OSS_global[p].a[_a.previous].next = OSS_global[p].a.size() - 1;
		OSS_global[p].a[_a.next].previous = OSS_global[p].a.size() - 1;
		OSS_global[p].size_num++;
	}

	static bool rank_compare(int a, int b){
		if (belong[a] < belong[b])
			return true;
		else return false;
	}

	//knn construction: many methods
	//naive_up
	void join_sbt_up_kpv_pro(int p, int x, int* pv, int* pvd, int pvsize){
		current_stamp++;
		current_distance[x] = 0;
		vector<int> b; b.clear();
		for (int i = 0; i < pvsize; i++){
			int org_t = pv[i];
			int t;
			int dis_orgt_x = pvd[i];
			int pdis = 0;
			int a_i = belong[pv[i]];
			for(int j = 1; j < (OSS[a_i].size_num+1); j++ ){
				up_cand_time++;
				t = OSS[a_i].a[j].key;
				pdis = dis_orgt_x + OSS[a_i].a[j].dist;
				if(current_state[t] != current_stamp){
					up_cand_size++;
					current_state[t] = current_stamp;
					current_distance[t] = pdis;
					b.push_back(t);
				}else if(current_distance[t] > pdis){
					current_distance[t] = pdis;
				}
			}
		}
		sort(b.begin(), b.end(), object_compare);
		
		if(is_current_object[x] == 1){
			OSS_push_front(p, x, 0);
			for (int i = 0; (i < b.size()) && (i < (M_K-1)); i++){
				OSS_push_back(p, b[i], current_distance[b[i]]);
        	}
		}else{
			for (int i = 0; (i < b.size()) && (i < M_K); i++){
				OSS_push_back(p, b[i], current_distance[b[i]]);
        	}
		}
	}
    void dfs_up_pro(){
		for(int p=n-1;p>=0;p--){
			join_sbt_up_kpv_pro(p, uniqueVertex[p], kpv[p], kpvd[p], kpvSize[p]);
		}
	}
	//naive_down
	Graph ug_con(int p){
		int x = uniqueVertex[p];
		Graph ug; 
		ug_stamp++;
		vector<pair<int,int>> v_dis; v_dis.clear(); //用来初始化 ug的E的
		ug.V.emplace_back(x); ug.Edge.emplace_back(v_dis); uVToUg_ID[x] = ug.n; ug.n++;
		queue<int> cug; cug.push(x);
		//kw*kh
		while(!cug.empty()){
			int p = belong[cug.front()], x = cug.front(); cug.pop();
			for(int i=0;i<kfvdSize[p];i++){
				int y = kfv[p][i]; 
				if (ug_affect[y] != ug_stamp){
					//cout<<"y: "<<y<<endl;
					ug.V.emplace_back(y); ug.Edge.emplace_back(v_dis); uVToUg_ID[y] = ug.n; ug.n++;
					cug.push(y);
					ug_affect[y] = ug_stamp;
				} 
				//边是一定要加的！！
				ug.Edge[uVToUg_ID[x]].emplace_back(uVToUg_ID[y], kfvd[p][i]);
				ug.Edge[uVToUg_ID[y]].emplace_back(uVToUg_ID[x], kfvd[p][i]);
			}
		}
		ug.D.clear();
		for (int i = 0; i < ug.n; i++)
			ug.D.push_back(ug.Edge[i].size());
		//cout<<"ug.n: "<<ug.n<< " max_ug_size: "<<max_ug_size<<" sum_ug_size: "<<sum_ug_size<<endl;
		if(max_ug_size < ug.n) max_ug_size = ug.n;
		sum_ug_size += ug.n;

		return ug;
	}
	vector<pair<int,int>> dijk_ug(Graph ug, int p){
		vector<pair<int,int>> v_dis; v_dis.clear();
		deque<int> cand; cand.clear();
		djik_stamp++;
		int x = uniqueVertex[p]; cand.push_back(x);
		current_distance[x] =0;
		while(!cand.empty()){
			x = cand.front(); cand.pop_front();
			distance_computed[x]=djik_stamp;
			v_dis.emplace_back(x, current_distance[x]);
			for(int i=0; i < ug.D[uVToUg_ID[x]]; i++){
				int cur_v = ug.V[ug.Edge[uVToUg_ID[x]][i].first];
				if(distance_computed[cur_v] != djik_stamp){
					if(djik_state[cur_v]==djik_stamp){
						int cur_dis = current_distance[x] + ug.Edge[uVToUg_ID[x]][i].second;
						if(cur_dis < current_distance[cur_v]) current_distance[cur_v] = cur_dis;
					}else{
						cand.push_back(cur_v);
						current_distance[cur_v] = current_distance[x] + ug.Edge[uVToUg_ID[x]][i].second;
						djik_state[cur_v] = djik_stamp;
					}
				}
			}
			sort(cand.begin(),cand.end(),object_compare);
		}
		int pin = 0;
		if(pin){
			cout<<uniqueVertex[p]<<" UG: ";
			for(int i=0;i<ug.n;i++){
				cout<<"( "<<v_dis[i].first<<", "<<v_dis[i].second<<") ";
			}
			cout<<endl;
		}
		return v_dis;
	}
	void ug_knn_con(int p){
		Graph ug = ug_con(p);
		vector<pair<int,int>> v_dist = dijk_ug(ug,p);
		current_stamp++;
		vector<int> b; b.clear();

		for(int i=0;i<ug.n;i++){
			int v = belong[v_dist[i].first], dist = v_dist[i].second;
			for(int k = 1; k <= OSS[v].size_num; k++){
				down_cand_time++;
				int t = OSS[v].a[k].key;
				int up_dist = dist + OSS[v].a[k].dist;
				if(current_state[t] != current_stamp){
					down_cand_size++;
					current_state[t] = current_stamp;
					current_distance[t] = up_dist;
					b.push_back(t);
				}else if(current_distance[t] > up_dist){
					current_distance[t] = up_dist;
				}
			}
		}
		sort(b.begin(), b.end(), object_compare);

        for (int i = 0; (i < b.size()) && (i < M_K); i++){
			OSS_global_push_back(p, b[i], current_distance[b[i]]);
        }
	}
	void knn_con_all(int n){
		for(int i=0;i<n;i++) ug_knn_con(i);
	}

	void test_h(int n){
		for(int i=0;i<n;i++) ug_con(i);
	}

    int up_cand_size, up_cand_time;
	long down_cand_size, down_cand_time;

	int max_ug_size, sum_ug_size;
	int *query_mark, query_mark_stamp;
	//max_dist + 取k次 
	void join_sbt_up_kpv_mdk(int p, int x, int* pv, int* pvd, int pvsize){
		current_stamp++;
		current_distance[x] = 0;
		vector<int> b; b.clear();
		int a_i;
		for (int i = 0; i < pvsize; i++){
			int org_t = pv[i];
			int t;
			int dis_orgt_x = pvd[i];
			int pdis = 0;
			a_i = belong[pv[i]];
			for(int j= OSS[a_i].a[0].next; j !=0; j = OSS[a_i].a[j].next){
				up_cand_time++;
				//cout<<"j: "<<j<<endl;
				t = OSS[a_i].a[j].key;
				pdis = dis_orgt_x + OSS[a_i].a[j].dist;
				if(current_state[t] != current_stamp){
					up_cand_size++;
					current_state[t] = current_stamp;
					current_distance[t] = pdis;
					b.push_back(t);
				}else if(current_distance[t] > pdis){
					current_distance[t] = pdis;
				}
				//cout<<"key: "<<t<<" dist: "<<current_distance[t]<<endl;
				//j = OSS[a_i].a[j].next;
			}
		}
		sort(b.begin(), b.end(), object_compare);
		
		if (is_current_object[x] == 1){
			OSS_push_back(p, x, 0);
			for (int i = 0; (i < b.size()) && (i < (M_K-1)); i++)
				OSS_push_back(p, b[i], current_distance[b[i]]);
		}else{
			for (int i = 0; (i < b.size()) && (i < M_K); i++)
				OSS_push_back(p, b[i], current_distance[b[i]]);
			
		}
	}
    void dfs_up_kvcmdk(){
		for(int p=n-1;p>=0;p--){
			join_sbt_up_kpv_mdk(p, uniqueVertex[p], kpv[p], kpvd[p], kpvSize[p]);
		}
	}
	void join_sbt_down_kfv_mdk(int p, int x, int* fv, int* fvd, int fvSize){
		query_mark_stamp++;
		int t, pdis;
		current_distance[x] = 0;

		int max_dis = INT_MAX;
		int a_i;
		//max_dist 计算
		for(int i = 0; i < fvSize-1; i++){
			a_i = belong[fv[i]];
			current_distance[fv[i]] = fvd[i];
			OSS_global[a_i].current = OSS_global[a_i].a[0].next;
			if(OSS_global[a_i].size_num != 0){
				down_cand_time++;
				down_cand_size++;
				int v = OSS_global[a_i].a[0].previous;
				int tmp = fvd[i] + OSS_global[a_i].a[v].dist; 
				if(tmp < max_dis) max_dis = tmp;
			}
		}
		//for(int j= OSS[a_i].a[0].next; j !=0; j = OSS[a_i].a[j].next)

		OSS[p].current = OSS[p].a[0].next;
		//取M_K次
		for (int j = 0; j < M_K; j++){
			int k = -1, dist_k = INT_MAX, i=0, _dist=0;
			for (i; i < fvSize-1; i++){
				int q = belong[fv[i]];
				while (OSS_global[q].current != 0 && query_mark[OSS_global[q].a[OSS_global[q].current].key] == query_mark_stamp){
					down_cand_time++;
					OSS_global[q].current = OSS_global[q].a[OSS_global[q].current].next;
				}
				if (OSS_global[q].current != 0){
					down_cand_size++;
					_dist = fvd[i] + OSS_global[q].a[OSS_global[q].current].dist;
					if (_dist > max_dis) OSS_global[q].current =0;
					if (k < 0 || (dist_k > _dist)){
						k = q; 
						dist_k = _dist;
					}
				}
			}
			if(i==fvSize-1){
				while (OSS[p].current != 0 && query_mark[OSS[p].a[OSS[p].current].key] == query_mark_stamp){
					down_cand_time++;
					OSS[p].current = OSS[p].a[OSS[p].current].next;
				}
				if (OSS[p].current != 0){
					down_cand_size++;
					_dist = OSS[p].a[OSS[p].current].dist;
					if (_dist > max_dis) OSS[p].current =0;
					if (k < 0 || (dist_k > _dist)){
						k = p; 
						dist_k = _dist;
					}
				} 
			}
			
			if (k < 0) break;
			int y = 0;
			if (k==p){ 
				y = OSS[k].a[OSS[k].current].key;
			}else{
				y = OSS_global[k].a[OSS_global[k].current].key;
			}
			
			OSS_global_push_back(p, y, dist_k);
			query_mark[y] = query_mark_stamp;
		}
		
	}
	void dfs_down_kvcmdk(){
		for(int p=0;p<n;p++){ 
			join_sbt_down_kfv_mdk(p, uniqueVertex[p], kfv[p], kfvd[p], kfvSize[p]);
		}
	}

	//knn index con
	void initialize_object_ori(){
		for (int i = 0; i < n; i++)
			for (int j = OSS[i].a[0].next; j != 0; j = OSS[i].a[j].next )
				delete_element(i, j);
		for (int i = 0; i < n; i++)
			for (int j = OSS_global[i].a[0].next; j != 0; j = OSS_global[i].a[j].next )
				delete_element_global(i, j);
		
		/*
		up_cand_size = 0;
		up_cand_time = 0;
		down_cand_size = 0;
		down_cand_time = 0;
		*/

		//max_ug_size=0, sum_ug_size=0;

		if(strcmp(insert_type, "pkkvc") == 0){
			//kvc naive: MAX_dis 做了一个bound + 取k次
			dfs_up_kvcmdk();
			dfs_down_kvcmdk();
		}else if(strcmp(insert_type, "naive") == 0){
			//naive kvc
			dfs_up_pro();
			knn_con_all(n);
		}else if(strcmp(insert_type, "h") == 0){
			//cout<<"hhhhh"<<endl;
			test_h(n);
		}

		
       
		//printf("max_ug_size: %.3lf \n", max_ug_size);	
		//printf("sum_ug_size: %.3lf \n", sum_ug_size);	
		/*
		cout<<"TREE SIZE: "<<TreeSize<<endl;
		float avg_up_cand_size = ((float)up_cand_size)/TreeSize;
		float avg_up_cand_time = ((float)up_cand_time)/TreeSize;
		float avg_down_cand_size = ((float)down_cand_size)/TreeSize;
		float avg_down_cand_time = ((float)down_cand_time)/TreeSize;

		//cout<<"max_ug_size: "<< max_ug_size<<endl;
		//cout<<"sum_ug_size: "<< sum_ug_size<<endl;
		//cout<<"ave_ug_size: "<< (sum_ug_size/TreeSize)<<endl;

		printf("avg_up_cand_size: %.3lf \n", avg_up_cand_size);	
		printf("avg_up_cand_time: %.3lf \n", avg_up_cand_time );	 
		printf("avg_down_cand_size: %.3lf \n", avg_down_cand_size);	
		printf("avg_down_cand_time: %.3lf \n", avg_down_cand_time);	
		*/

        /* vector<int> a;
		if (strcmp(insert_type, "sort") == 0) dfs_sort(root, a);
		else{
			dfs_up(root);
			dfs_down(root);
		}*/
	}

    //initialize obj, pointers and vectors
	void object_setting_ori(int n){
		is_current_object = (int*)malloc(sizeof(int) * (n + 1));

		current_distance = (int*)malloc(sizeof(int) * (n + 1));
		current_state = (int*)malloc(sizeof(int) * (n + 1));
		current_stamp = 0;

		/*
		djik_state = (int*)malloc(sizeof(int) * (n + 1));
		distance_computed = (int*)malloc(sizeof(int) * (n + 1));
		ug_stamp = 0;

		ug_affect = (int*)malloc(sizeof(int) * (n + 1));
		uVToUg_ID = (int*)malloc(sizeof(int) * (n + 1));
		*/
		query_mark = (int*)malloc(sizeof(int) * (n + 1));
		query_mark_stamp = 0;

		for (int i = 0; i <= n; i++){
			current_state[i] = 0;
			is_current_object[i] = 0;
			query_mark[i] = 0;
			/*
			ug_affect[i] = 0;
			uVToUg_ID[i] = 0;
			djik_state[i] = 0;
			distance_computed[i] = 0;*/
		}
	}

	void update_setting(int n){
		//compute_object_number();
		curup_dist = (int*)malloc(sizeof(int) * (n + 1));
		curup_affect = (int*)malloc(sizeof(int) * (n + 1));
		curup_stamp = 0;
		query_mark_stamp = 0;
		for (int i = 0; i <= n; i++){
			//current_state[i] = 0;
			query_mark[i] = 0;
			curup_dist[i] = 99999999;
			curup_affect[i] = 0;
		}
	}

	void process_insert(int p, int x, int disxy){
		int i = 0;
		while (OSS_global[p].a[i].previous != 0 && OSS_global[p].a[OSS_global[p].a[i].previous].dist > disxy)
            i = OSS_global[p].a[i].previous;
		
        list_node _a;
        _a.next = i;
        _a.previous = OSS_global[p].a[i].previous;
        _a.key = x;
        _a.dist = disxy;
		
        //双向链表更新一下
		OSS_global[p].a.push_back(_a);
		OSS_global[p].a[_a.previous].next = OSS_global[p].a.size() - 1;
		OSS_global[p].a[_a.next].previous = OSS_global[p].a.size() - 1;
		OSS_global[p].size_num++;
		if (OSS_global[p].size_num > M_K * RESERVE_TIME)
			delete_element_global(p, OSS_global[p].a[0].previous);

		//cout<<uniqueVertex[p]<<" G: ";
		//for (int i = OSS_global[p].a[0].next; i != 0; i = OSS_global[p].a[i].next){
		//	cout<<" ("<<OSS_global[p].a[i].key<<", "<<OSS_global[p].a[i].dist<<") ";
		//}
		//cout<<endl;
		
	}

	bool check_insert(int cur, int x, int disxy){
		//if (OSS_global[cur].size_num >= M_K){
		if (OSS_global[cur].a[OSS_global[cur].a[0].previous].dist <= disxy){
			return false; 
		}
		//}
		return true;
	}

	void insert_ori(int x){
		if(is_current_object[x] == 1) return; 
		is_current_object[x] = 1;
		
		int p = belong[x];
        curup_stamp++;
		curup_dist[x] = 0; current_state[x] = curup_stamp;
		
		queue<int> cug; curup_affect[x] = curup_stamp;

		int mark_size=0;
		vector<int> mark; mark.clear(); mark.emplace_back(x); mark_size++;
		
		for(int i=0;i<kfvdSize[p];i++){
			curup_dist[kfv[p][i]] = kfvd[p][i]; current_state[kfv[p][i]] = curup_stamp;
			if(check_insert(belong[kfv[p][i]],x,curup_dist[kfv[p][i]])) {
				cug.push(kfv[p][i]);
				curup_affect[kfv[p][i]] = curup_stamp;
				mark.emplace_back(kfv[p][i]); mark_size++;
			}
		}

		for(int i=0;i<kpvdSize[p];i++){
			curup_dist[kpv[p][i]] = kpvd[p][i]; current_state[kpv[p][i]] = curup_stamp; 
			if(check_insert(belong[kpv[p][i]],x,curup_dist[kpv[p][i]])) {
				cug.push(kpv[p][i]);
				curup_affect[kpv[p][i]] = curup_stamp;
				mark.emplace_back(kpv[p][i]); mark_size++;
			}
		}

		while(!cug.empty()){
			p = belong[cug.front()]; cug.pop();
			
			for(int i=0;i<kfvdSize[p];i++){
				
				if (current_state[kfv[p][i]] == curup_stamp) {
					int cur_dis = kfvd[p][i] + curup_dist[uniqueVertex[p]];
					if (curup_dist[kfv[p][i]] > cur_dis){
						curup_dist[kfv[p][i]] = cur_dis;
						if(curup_affect[kfv[p][i]] == curup_stamp) continue;
						if(check_insert(belong[kfv[p][i]],x,curup_dist[kfv[p][i]])){ curup_affect[kfv[p][i]] = curup_stamp;cug.push(kfv[p][i]); mark.emplace_back(kfv[p][i]); mark_size++;} 
					} 
				}else{
					curup_dist[kfv[p][i]] = kfvd[p][i] + curup_dist[uniqueVertex[p]];
					current_state[kfv[p][i]] = curup_stamp;
					if(check_insert(belong[kfv[p][i]],x,curup_dist[kfv[p][i]])){
						cug.push(kfv[p][i]);
						curup_affect[kfv[p][i]] = curup_stamp;
						mark.emplace_back(kfv[p][i]); mark_size++;
					} 
				}
			}
			for(int i=0;i<kpvdSize[p];i++){
				if (current_state[kpv[p][i]] == curup_stamp) {
					int cur_dis = kpvd[p][i] + curup_dist[uniqueVertex[p]];
					if (curup_dist[kpv[p][i]] > cur_dis) {
						curup_dist[kpv[p][i]] = cur_dis;
						if(curup_affect[kpv[p][i]] == curup_stamp) continue;
						if(check_insert(belong[kpv[p][i]],x,curup_dist[kpv[p][i]])){ curup_affect[kpv[p][i]] = curup_stamp;cug.push(kpv[p][i]); mark.emplace_back(kpv[p][i]); mark_size++;}  
					}
				}else{
					curup_dist[kpv[p][i]] = kpvd[p][i] + curup_dist[uniqueVertex[p]];
					current_state[kpv[p][i]] = curup_stamp;
					if(check_insert(belong[kpv[p][i]],x,curup_dist[kpv[p][i]])) {
						curup_affect[kpv[p][i]] = curup_stamp;
						cug.push(kpv[p][i]);
						mark.emplace_back(kpv[p][i]); mark_size++;
					}
				}
			}
		}
		
		for(int i=0;i<mark_size;i++){
			process_insert(belong[mark[i]],x,curup_dist[mark[i]]);
		}
	}

	void insert(int x){
		if(is_current_object[x] == 1) return; 
		is_current_object[x] = 1;

        curup_stamp++;
		curup_dist[x] = 0; current_state[x] = curup_stamp;
		
		queue<int> cug; curup_affect[x] = curup_stamp; cug.push(x);

		int mark_size=0;
		vector<int> mark; mark.clear(); mark.emplace_back(x); mark_size++;

		int p,cur_dis;
		while(!cug.empty()){
			p = belong[cug.front()]; cug.pop();
			for(int i=0;i<kfvdSize[p];i++){
				if (current_state[kfv[p][i]] == curup_stamp) {
					cur_dis = kfvd[p][i] + curup_dist[uniqueVertex[p]];
					if (curup_dist[kfv[p][i]] > cur_dis){
						curup_dist[kfv[p][i]] = cur_dis;
						if(curup_affect[kfv[p][i]] == curup_stamp) continue;
						if(check_insert(belong[kfv[p][i]],x,curup_dist[kfv[p][i]])){ 
							curup_affect[kfv[p][i]] = curup_stamp;
							cug.push(kfv[p][i]); mark.emplace_back(kfv[p][i]); mark_size++;
						} 
					} 
				}else{
					curup_dist[kfv[p][i]] = kfvd[p][i] + curup_dist[uniqueVertex[p]];
					current_state[kfv[p][i]] = curup_stamp;
					if(check_insert(belong[kfv[p][i]],x,curup_dist[kfv[p][i]])){
						curup_affect[kfv[p][i]] = curup_stamp;
						cug.push(kfv[p][i]); mark.emplace_back(kfv[p][i]); mark_size++;
					} 
				}
			}
			for(int i=0;i<kpvdSize[p];i++){
				if (current_state[kpv[p][i]] == curup_stamp) {
					int cur_dis = kpvd[p][i] + curup_dist[uniqueVertex[p]];
					if (curup_dist[kpv[p][i]] > cur_dis) {
						curup_dist[kpv[p][i]] = cur_dis;
						if(curup_affect[kpv[p][i]] == curup_stamp) continue;
						if(check_insert(belong[kpv[p][i]],x,curup_dist[kpv[p][i]])){ 
							curup_affect[kpv[p][i]] = curup_stamp;
							cug.push(kpv[p][i]); mark.emplace_back(kpv[p][i]); mark_size++;}  
					}
				}else{
					curup_dist[kpv[p][i]] = kpvd[p][i] + curup_dist[uniqueVertex[p]];
					current_state[kpv[p][i]] = curup_stamp;
					if(check_insert(belong[kpv[p][i]],x,curup_dist[kpv[p][i]])) {
						curup_affect[kpv[p][i]] = curup_stamp;
						cug.push(kpv[p][i]);
						mark.emplace_back(kpv[p][i]); mark_size++;
					}
				}
			}
		}
		
		for(int i=0;i<mark_size;i++){
			process_insert(belong[mark[i]],x,curup_dist[mark[i]]);
		}
	}

	void process_delete(int p, int up_stamp){
		query_mark_stamp++;

		int max_dis = INT_MAX;
		int a_i;

		for(int i=0;i<kfvdSize[p];i++){
			a_i = belong[kfv[p][i]];
			OSS_global[a_i].current = OSS_global[a_i].a[0].next;
			if(OSS_global[a_i].size_num != 0){
				int tmp = kfvd[p][i] + OSS_global[a_i].a[OSS_global[a_i].a[0].previous].dist; 
				if(tmp < max_dis) max_dis = tmp;
			}
		}

		for(int i=0;i<kpvdSize[p];i++){
			a_i = belong[kpv[p][i]];
			OSS_global[a_i].current = OSS_global[a_i].a[0].next;
			if(curup_affect[kpv[p][i]] == up_stamp) continue;
			if(OSS_global[a_i].size_num != 0){
				int tmp = kpvd[p][i] + OSS_global[a_i].a[ OSS_global[a_i].a[0].previous].dist; 
				if(tmp < max_dis) max_dis = tmp;
			}
		}
		//取M_K次
		int q, y;
		int topk = M_K;
		if(is_current_object[uniqueVertex[p]] == 1) {
			topk--;
			query_mark[uniqueVertex[p]] = query_mark_stamp;
		}
		for (int j = 0; j < M_K; j++){
			int k = -1, dist_k = INT_MAX, _dist=0;
			for (int i=0; i < kfvdSize[p]; i++){
				q = belong[kfv[p][i]];
				while (OSS_global[q].current != 0 && query_mark[OSS_global[q].a[OSS_global[q].current].key] == query_mark_stamp){
					OSS_global[q].current = OSS_global[q].a[OSS_global[q].current].next;
				}
				if (OSS_global[q].current != 0){
					_dist = kfvd[p][i] + OSS_global[q].a[OSS_global[q].current].dist;
					if (_dist > max_dis) OSS_global[q].current =0;
					if (k < 0 || (dist_k > _dist)){
						k = q; 
						dist_k = _dist;
					}
				}
			}
			for (int i=0; i < kpvdSize[p]; i++){
				q = belong[kpv[p][i]];
				while (OSS_global[q].current != 0 && query_mark[OSS_global[q].a[OSS_global[q].current].key] == query_mark_stamp){
					OSS_global[q].current = OSS_global[q].a[OSS_global[q].current].next;
				}
				if (OSS_global[q].current != 0){
					_dist = kpvd[p][i] + OSS_global[q].a[OSS_global[q].current].dist;
					if (_dist > max_dis) OSS_global[q].current =0;
					if (k < 0 || (dist_k > _dist)){
						k = q; 
						dist_k = _dist;
					}
				}
			}
			
			if (k < 0) break;
			
			y = OSS_global[k].a[OSS_global[k].current].key;
			
			query_mark[y] = query_mark_stamp;
			if(j==(topk-1)) {
				OSS_global_push_back(p, y, dist_k); break;
			}
		}
	}

	bool check_delete(int p, int x, int disxy){
		if (OSS_global[p].a[OSS_global[p].a[0].previous].dist < disxy){
            return false;
        }
		int i = 0;
		while (OSS_global[p].a[i].next != 0 && OSS_global[p].a[OSS_global[p].a[i].next].dist <= disxy){
			i = OSS_global[p].a[i].next;
			if (OSS_global[p].a[i].key == x) {
				delete_element_global(p, i);
				return true;
			} 
		}
		return false;
	}

	void obj_delete(int x){
		if(is_current_object[x] == 0) return; 
		is_current_object[x] = 0;

		curup_stamp++;
		curup_dist[x] = 0; current_state[x] = curup_stamp;

		delete_element_global(belong[x], OSS_global[belong[x]].a[0].next);

		queue<int> cug; curup_affect[x] = curup_stamp; cug.push(x);

		int mark_size=0;
		vector<int> mark; mark.clear(); mark.emplace_back(x); mark_size++;

		int p;
		while(!cug.empty()){
			p = belong[cug.front()]; cug.pop();
			for(int i=0;i<kfvdSize[p];i++){
				if (current_state[kfv[p][i]] == curup_stamp) {
					int cur_dis = kfvd[p][i] + curup_dist[uniqueVertex[p]];
					if (curup_dist[kfv[p][i]] > cur_dis){
						curup_dist[kfv[p][i]] = cur_dis;
						if(curup_affect[kfv[p][i]] == curup_stamp) continue;
						if(check_delete(belong[kfv[p][i]],x,curup_dist[kfv[p][i]])){ curup_affect[kfv[p][i]] = curup_stamp; cug.push(kfv[p][i]); mark.emplace_back(kfv[p][i]); mark_size++;} 
					} 
				}else{
					curup_dist[kfv[p][i]] = kfvd[p][i] + curup_dist[uniqueVertex[p]];
					current_state[kfv[p][i]] = curup_stamp;
					if(check_delete(belong[kfv[p][i]],x,curup_dist[kfv[p][i]])){
						cug.push(kfv[p][i]);
						curup_affect[kfv[p][i]] = curup_stamp;
						mark.emplace_back(kfv[p][i]); mark_size++;
					} 
				}
			}
			for(int i=0;i<kpvdSize[p];i++){
				if (current_state[kpv[p][i]] == curup_stamp) {
					int cur_dis = kpvd[p][i] + curup_dist[uniqueVertex[p]];
					if (curup_dist[kpv[p][i]] > cur_dis) {
						curup_dist[kpv[p][i]] = cur_dis;
						if(curup_affect[kpv[p][i]] == curup_stamp) continue;
						if(check_delete(belong[kpv[p][i]],x,curup_dist[kpv[p][i]])){ 
						curup_affect[kpv[p][i]] = curup_stamp;
						cug.push(kpv[p][i]); mark.emplace_back(kpv[p][i]); mark_size++;}  
					}
				}else{
					curup_dist[kpv[p][i]] = kpvd[p][i] + curup_dist[uniqueVertex[p]];
					current_state[kpv[p][i]] = curup_stamp;
					if(check_delete(belong[kpv[p][i]],x,curup_dist[kpv[p][i]])) {
						curup_affect[kpv[p][i]] = curup_stamp;
						cug.push(kpv[p][i]);
						mark.emplace_back(kpv[p][i]); mark_size++;
					}
				}
			}
		}
		

		sort(mark.begin(), mark.end(), rank_compare);
		for(int i=0;i<mark_size;i++){
			process_delete(belong[mark[i]],curup_stamp);
		}

	}

    void traversal(int p){
        int x = uniqueVertex[p];
        if (is_current_object[x] == 1)
            object_number[p]++;
		for (int i = 0; i < chSize[p]; i++){
			object_number[p] += object_number[ch[p][i]];
    	}
		
		//cout<<x<<" object_number[p]: "<< object_number[p]<<endl;
    }
	
    void compute_object_number(){
		//初始化 object_number
        object_number.resize(n + 1);
        for (int i = 0; i <= n; i++)
            object_number[i] = 0;
		for(int i=n-1;i>=0;i--){
			traversal(i);
		}
    }
	
	vector<pair<int,int>> query(int x, int top_k){
		vector<pair<int,int>> result; result.clear();

		int p = belong[x];
		//cout<<x<<" G: ";
		for (int i = OSS_global[p].a[0].next; i != 0; i = OSS_global[p].a[i].next){
			result.push_back(make_pair(OSS_global[p].a[i].key, OSS_global[p].a[i].dist));
			//cout<<" ("<<OSS_global[p].a[i].key<<", "<<OSS_global[p].a[i].dist<<") ";
		}
		//cout<<endl;
		return result;
	}

	void print(int root){
		int cnt = 0;
		for (int i = OSS[root].a[0].next; i != 0; i = OSS[root].a[i].next){
			printf("(%d, %d)", OSS[root].a[i].key, OSS[root].a[i].dist);
			cnt++;
		}
		printf("\n");

		if (cnt != OSS[root].size_num){
			cout << "cnt: " << cnt << "     " << "OSS[root].size_num: " << OSS[root].size_num << endl; 
			while (1);
		}
	}
	bool double_objects(int p){
		for (int i = OSS[p].a[0].next; i != 0; i = OSS[p].a[i].next){
			//a[i] 当前的 和 previous 如果相等-- 是一个vertex，return false
			//a[i] 不重复
			if (OSS[p].a[i].key == OSS[p].a[OSS[p].a[i].previous].key)
				return false;
			//保证a[i] dist 前面的比后面的小
			if ((OSS[p].a[i].next != 0) && (OSS[p].a[OSS[p].a[i].next].dist < OSS[p].a[i].dist))
				return false;
		}
		return true;
	}
	bool double_objects_ori(int p){
		for (int i = OSS_global[p].a[0].next; i != 0; i = OSS_global[p].a[i].next){
			//a[i] 当前的 和 previous 如果相等-- 是一个vertex，return false
			//a[i] 不重复
			if (OSS_global[p].a[i].key == OSS_global[p].a[OSS_global[p].a[i].previous].key)
				return false;
			//保证a[i] dist 前面的比后面的小
			if ((OSS_global[p].a[i].next != 0) && (OSS_global[p].a[OSS_global[p].a[i].next].dist < OSS_global[p].a[i].dist))
				return false;
		}
		return true;
	}
	bool check_everyone(){
		for (int i = 0; i < n; i++){
			if (double_objects(i) == false) return false;
			if (double_objects_ori(i) == false) return false;
		}
			
		return true;
	}
	
	void stop(){
		while (1);
	}
	void cnt_oss(){
		int sum_oss = 0;
		int max_oss = 0;
		int oss_size = 0;

		for (int i=0;i<TreeSize;i++){
			oss_size = OSS_global[i].size_num;
			if(oss_size > max_oss) max_oss = oss_size;
			sum_oss += oss_size;
		}
		cout<<"sum_oss: "<<sum_oss<<endl;
		cout<<"max_oss: "<<max_oss<<endl;
		cout<<"avg_oss: "<<(sum_oss/TreeSize)<<endl;
	}
};

double get_mean_double(vector<double>& times) {
    double mean = 0.0;
    for (double val : times) {
        mean += val;
    }
	//cout<<mean<<endl;
	//printf("query all time: %.6lf s \n", mean);	
	//cout<<"times.size(): "<<times.size()<<endl;
    return mean / times.size();
}
double get_var_double(vector<double>& times) {
    double mean = get_mean_double(times);
    double var = 0.0;
    for (double val : times) {
        var += (val - mean) * (val - mean);
    }
    return var / times.size();
}


	//-- query NY-t.index NY-t.obj NY-t.query sort
	//./querychscomp example.index1 example.obj1 -q example.query2 neighbor optimal vcgindex topk
	//                 1              2          3   4              5        6       7        8
	//out_file 改三处: main // kNN // vct
int main(int argc, char *argv[]){
	srand((int)(time(0)));
    cout << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << " " << argv[6] << " " << argv[7] << endl;
	readIndex(argv[1]);
	//cout<<"n: "<<n<<endl;
	//cout<<"finish readindex"<<endl;

    //----------------------prepare-------------
	
    insert_type = argv[5];
    query_type = argv[6];
    
	//-------------test------------------
	
	int topk = atoi(argv[7]); //atoi(argv[8]);
    cout<<"topk: "<<topk<<endl;
	kNN knn(topk);
	//initialize knn structure
	knn.create_kNN_index();

	FILE *fobj = fopen(argv[2], "r");
	int number_object;
	fscanf(fobj, "%d", &number_object);
	double start_time = GetTime();
	//initialize pointer, vector, objects
	knn.object_setting_ori(n);
	for (int i = 0; i < number_object; i++){
		int x;
		fscanf(fobj, "%d", &x);
		is_current_object[x] = 1;
	}
	
	int cnt_object = 0;
	for (int i = 1; i <= n; i++) if(is_current_object[i] == 1) cnt_object ++;
	
	//construct knn
	knn.initialize_object_ori();

	double end_time = GetTime();	
    //printf("object initialization time: %.6lf ms\n", (end_time - start_time) * 1e3);
	printf("object initialization time: %.6lf s\n", (end_time - start_time) );

	//----------------------query------------------------------
	
    vector<double> times;
    times.clear();
    for (int i = 0; i < knn.period; i++)
        knn.times_period[i].clear();
	
	if (argv[3][1] == 'q'){
		FILE *fquery = fopen(argv[4], "r");
		int q_n;
		// query number
		fscanf(fquery, "%d", &q_n);
		//printf("n: %d\n", n);
		
		knn.query_mark_stamp = 0;
		
		for (int i = 0; i <= n; i++)
			knn.query_mark[i] = 0;
		
		start_time = GetTime();
	    vector<double> time_array;
	    time_array.clear();
        
		double cnt_delay_time = 0;
		//tmp_dis = (int*)malloc(sizeof(int) * (n + 1));
		vector<pair<int,int>> res;
		//double _start_time = GetTime();
		for (int i = 0; i < q_n; i++){
			int x, k;
			fscanf(fquery, "%d %d", &x, &k);
            double _start_time = GetTime();
            start_time = GetTime();
            res = knn.query(x, k);
            times.push_back(GetTime() - start_time);
			
		}
		//end_time = GetTime();	
	    //double ave_time = (end_time - _start_time) / q_n;
        printf("Average query time: %.6lf us\n", get_mean_double(times) * 1e6);
		//cout<<"ave_time : "<< (ave_time * 1e6) <<endl;
        printf("Var query time: %.6lf us\n", get_var_double(times) * 1e6);
    	
	}
	else if (argv[3][1] == 'u'){
		knn.update_setting(n);
		FILE *fupdate = fopen(argv[4], "r");
		int u_n;
		fscanf(fupdate, "%d", &u_n);
		cout<<"u_n: "<<u_n<<endl;
		char st[20];
		//u_n = 100;
		double sum_time = 0.0;	
		for (int i = 0; i < u_n; i++){
			fscanf(fupdate, "%s", st);
			int x;
			fscanf(fupdate, "%d", &x);
            start_time = GetTime();
			if (st[0] == 'i'){
				//knn.insert(x); 
				knn.insert_ori(x); 
			}else{
				//cout<<"com dele"<<endl;
				knn.obj_delete(x);
			}
			sum_time += (GetTime() - start_time);
		}
		double ave_time = sum_time / u_n;
		cout<<"Average Update Time(us) : "<< (ave_time * 1e6) <<endl;
		
	}
    else{
		//update 重头算update之后的index
		knn.update_setting(n);
		FILE *fupdate = fopen(argv[4], "r");
		int u_n;
		fscanf(fupdate, "%d", &u_n);
		char st[20];

		start_time = GetTime();
		for (int i = 0; i < u_n; i++){
			fscanf(fupdate, "%s", st);
			int x;
			fscanf(fupdate, "%d", &x);
			if (st[0] == 'i'){
				is_current_object[x] = 1;
			}
			else {
				is_current_object[x] = 0;
			}
		}
		
		knn.initialize_object_ori();
	
		//	knn.print(root);
		if (knn.check_everyone()){
			printf("right\n");
		}else printf("wrong\n");
	}

	//fclose(knn.fout);

}
