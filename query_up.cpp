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
#define INT_MAX 999999999
#define RESERVE_TIME 1

const int SIZEOFINT = 4;

char * insert_type;

int *uniqueVertex; // order -> original id
int *belong;

int **kpv;
int *kpvSize;

int **kpvd;
int *kpvdSize;

int **kfv;
int *kfvSize;

int **kfvd;
int *kfvdSize;

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

FILE *fin;

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
void readIndex(char *file){

	fin = fopen(file, "rb");
	fread(&n, SIZEOFINT, 1, fin);

	kpv = (int**)malloc(sizeof(int*) * (n));
	kpvSize = (int*)malloc(sizeof(int) * (n));
	kpvd = (int**)malloc(sizeof(int*) * (n));
	kpvdSize = (int*)malloc(sizeof(int) * (n));
	kfv = (int**)malloc(sizeof(int*) * (n));
	kfvSize = (int*)malloc(sizeof(int) * (n));
	kfvd = (int**)malloc(sizeof(int*) * (n));
	kfvdSize = (int*)malloc(sizeof(int) * (n));

	int a = 0;
	fread(&a, SIZEOFINT, 1, fin);
	uniqueVertex = (int*)malloc(sizeof(int) * (a));
	for (int i = 0; i < a; i++){
		fread(&uniqueVertex[i], SIZEOFINT, 1, fin);
	}

	fread(&a, SIZEOFINT, 1, fin);
	belong = (int*)malloc(sizeof(int) * (a));
	for (int i = 0; i < a; i++){
		fread(&belong[i], SIZEOFINT, 1, fin);
	}

	for(int x=0;x<n;++x){
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
	}
	fclose(fin);
}
int *is_current_object, *current_distance, current_stamp, *current_state;

int *curup_dist, curup_stamp, *curup_affect, *curup_fvSize;

int *ug_affect, ug_stamp, *uVToUg_ID; 

int djik_stamp, *djik_state, *distance_computed;

class kNN{
public:	
int M_K;

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

    
    vector<double> times_period[10];
    int period = 8;
	vector<object_saveing_structure> OSS;
 	vector<object_saveing_structure> OSS_global;
    vector<int> object_number;
	
	//initailize knn structure
	void create_kNN_index(){
		//initialize OSS
		OSS.clear();
		for (int i = 0; i < n; i++){
			object_saveing_structure oss;
			OSS.push_back(oss);
		}
		//initialize OSS_global
		OSS_global.clear();
		for (int i = 0; i < n; i++){
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
		++current_stamp;
		current_distance[x] = 0;
		vector<int> b; b.clear();
		for (int i = 0; i < pvsize; i++){
			int org_t = pv[i];
			int t;
			int dis_orgt_x = pvd[i];
			int pdis = 0;
			int a_i = belong[pv[i]];
			for(int j = 1; j < (OSS[a_i].size_num+1); j++ ){
				
				t = OSS[a_i].a[j].key;
				pdis = dis_orgt_x + OSS[a_i].a[j].dist;
				if(current_state[t] != current_stamp){
					
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
		for(int p=0;p<n;++p){
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
				int t = OSS[v].a[k].key;
				int up_dist = dist + OSS[v].a[k].dist;
				if(current_state[t] != current_stamp){
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
				
				//cout<<"j: "<<j<<endl;
				t = OSS[a_i].a[j].key;
				pdis = dis_orgt_x + OSS[a_i].a[j].dist;
				if(current_state[t] != current_stamp){
					
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
			for (int i = 0; (i < b.size()) && (i < (M_K-1)); i++){
				OSS_push_back(p, b[i], current_distance[b[i]]);
			}
		}else{
			for (int i = 0; (i < b.size()) && (i < M_K); i++)
				OSS_push_back(p, b[i], current_distance[b[i]]);
			
		}
	}
    void dfs_up_kvcmdk(){
		for(int p=0;p<n;++p){
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
				int v = OSS_global[a_i].a[0].previous;
				int tmp = fvd[i] + OSS_global[a_i].a[v].dist; 
				if(tmp < max_dis) max_dis = tmp;
			}
		}

		OSS[p].current = OSS[p].a[0].next;
		//取M_K次
		for (int j = 0; j < M_K; j++){
			int k = -1, dist_k = INT_MAX, i=0, _dist=0;
			for (i; i < fvSize-1; i++){
				int q = belong[fv[i]];
				while (OSS_global[q].current != 0 && query_mark[OSS_global[q].a[OSS_global[q].current].key] == query_mark_stamp){
					OSS_global[q].current = OSS_global[q].a[OSS_global[q].current].next;
				}
				if (OSS_global[q].current != 0){
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
					OSS[p].current = OSS[p].a[OSS[p].current].next;
				}
				if (OSS[p].current != 0){
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
		for(int p=n-1;p>=0;--p){
			//cout<<"p: "<<p<<" x: "<<uniqueVertex[p]<<endl;
			join_sbt_down_kfv_mdk(p, uniqueVertex[p], kfv[p], kfvd[p], kfvSize[p]);
		}
	}


	//initialize obj, pointers and vectors
	void object_setting(int n){
		current_distance = (int*)malloc(sizeof(int) * (n + 1));
		current_state = (int*)malloc(sizeof(int) * (n + 1));
		current_stamp = 0;
		query_mark = (int*)malloc(sizeof(int) * (n + 1));
		query_mark_stamp = 0;

		for (int i = 0; i <= n; i++){
			current_state[i] = 0;
			query_mark[i] = 0;
		}
	}

	void object_setting_naive(int n){
		current_distance = (int*)malloc(sizeof(int) * (n + 1));
		current_state = (int*)malloc(sizeof(int) * (n + 1));
		current_stamp = 0;
		
		djik_state = (int*)malloc(sizeof(int) * (n + 1));
		distance_computed = (int*)malloc(sizeof(int) * (n + 1));
		ug_stamp = 0;

		ug_affect = (int*)malloc(sizeof(int) * (n + 1));
		uVToUg_ID = (int*)malloc(sizeof(int) * (n + 1));
		
		query_mark = (int*)malloc(sizeof(int) * (n + 1));
		query_mark_stamp = 0;

		for (int i = 0; i <= n; i++){
			current_state[i] = 0;
			query_mark[i] = 0;
			
			ug_affect[i] = 0;
			uVToUg_ID[i] = 0;
			djik_state[i] = 0;
			distance_computed[i] = 0;
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
		
		//max_ug_size=0, sum_ug_size=0;

		if(strcmp(insert_type, "opt") == 0){
			object_setting(n);
			//kvc naive: MAX_dis 做了一个bound + 取k次
			dfs_up_kvcmdk();
			dfs_down_kvcmdk();
		}else if(strcmp(insert_type, "pri") == 0){
			//naive kvc
			object_setting_naive(n);
			dfs_up_pro();
			knn_con_all(n);
		}

		//printf("max_ug_size: %.3lf \n", max_ug_size);	
		//printf("sum_ug_size: %.3lf \n", sum_ug_size);	
	}

    
	void update_setting(int n){
		
		curup_dist = (int*)malloc(sizeof(int) * (n + 1));
		curup_affect = (int*)malloc(sizeof(int) * (n + 1));
		curup_stamp = 0;
		query_mark_stamp = 0;
		for (int i = 0; i <= n; i++){
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
		if (OSS_global[cur].a[OSS_global[cur].a[0].previous].dist <= disxy){
			return false; 
		}
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

	vector<pair<int,int>> query(int x, int top_k){
		vector<pair<int,int>> result; result.clear();
		int p = belong[x];
		//cout<<x<<" G: ";
		for (int i = OSS_global[p].a[0].next; i != 0; i = OSS_global[p].a[i].next){
			result.push_back(make_pair(OSS_global[p].a[i].key, OSS_global[p].a[i].dist));
			//cout<<"("<<OSS_global[p].a[i].key<<", "<<OSS_global[p].a[i].dist<<") ";
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

		for (int i=0;i<n;i++){
			oss_size = OSS_global[i].size_num;
			if(oss_size > max_oss) max_oss = oss_size;
			sum_oss += oss_size;
		}
		cout<<"sum_oss: "<<sum_oss<<endl;
		cout<<"max_oss: "<<max_oss<<endl;
		cout<<"avg_oss: "<<(sum_oss/n)<<endl;
	}
};

double get_mean_double(vector<double>& times) {
    double mean = 0.0;
    for (double val : times) {
        mean += val;
    }
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
//./querychscomp example.index1 example.obj1 -q example.query2 con_method topk
//                 1              2          3   4              5          6   

int main(int argc, char *argv[]){
	srand((int)(time(0)));
    cout << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << " " << argv[6] << endl;
	readIndex(argv[1]);
	cout<<"n: "<<n<<endl;
	cout<<"finish readindex"<<endl;

    //----------------------prepare-------------
	
    insert_type = argv[5];
    
	//-------------test------------------
	
	int topk = atoi(argv[6]);
    cout<<"topk: "<<topk<<endl;
	kNN knn(topk);
	//initialize knn structure
	knn.create_kNN_index();

	FILE *fobj = fopen(argv[2], "r");
	int number_object;
	fscanf(fobj, "%d", &number_object);
	double start_time = GetTime();
	//initialize pointer, vector, objects
	//knn.object_setting_ori(n);

	is_current_object = (int*)malloc(sizeof(int) * (n + 1));
	for(int i=0;i<=n;++i) is_current_object[i] = 0;

	for (int i = 0; i < number_object; i++){
		int x;
		fscanf(fobj, "%d", &x);
		is_current_object[x] = 1;
	}
	
	int cnt_object = 0;
	for (int i = 1; i <= n; i++) if(is_current_object[i] == 1) ++cnt_object;
	
	//construct knn
	knn.initialize_object_ori();

	double end_time = GetTime();	
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
		printf("q_n: %d\n", q_n);
		
		knn.query_mark_stamp = 0;
		
		for (int i = 0; i <= n; i++) knn.query_mark[i] = 0;
		
		start_time = GetTime();
	    vector<double> time_array;
	    time_array.clear();
		double cnt_delay_time = 0;
		vector<pair<int,int>> res;
		for (int i = 0; i < q_n; i++){
			int x, k;
			fscanf(fquery, "%d %d", &x, &k);
            double _start_time = GetTime();
            start_time = GetTime();
            res = knn.query(x, k);
            times.push_back(GetTime() - start_time);
		}
        printf("Average query time: %.6lf us\n", get_mean_double(times) * 1e6);
        printf("Var query time: %.6lf us\n", get_var_double(times) * 1e6);
	}
	else if (argv[3][1] == 'u'){
		knn.update_setting(n);
		FILE *fupdate = fopen(argv[4], "r");
		int u_n;
		fscanf(fupdate, "%d", &u_n);
		char st[20];
		double sum_time = 0.0;	
		for (int i = 0; i < u_n; i++){
			fscanf(fupdate, "%s", st);
			int x;
			fscanf(fupdate, "%d", &x);
            start_time = GetTime();
			if (st[0] == 'i'){
				knn.insert_ori(x); 
			}else{
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
	
		if (knn.check_everyone()){
			printf("right\n");
		}else printf("wrong\n");
	}

	//fclose(knn.fout);

}
