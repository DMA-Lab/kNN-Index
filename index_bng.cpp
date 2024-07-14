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
#include "flatten_hash_map.h"
using namespace std;

clock_t ct;
int cnt, tree_width = 0;
vector<vector<pair<int,int>>> Edge;
// initialize graph data
struct Graph{
	int n, m;
	vector<int> V;
	vector<map<int,int>> E;
	vector<int> D;
	Graph(){
		n = m = 0;
		V.clear();
		E.clear();
		Edge.clear();
	}

	Graph(char *file){
		//	cout << "file:" << file << endl;
		Graph();
		FILE *fin = fopen(file, "r");
		fscanf(fin, "%d", &n);
		fscanf(fin, "%d", &m);
		printf("n m: %d %d\n", n, m);

		Edge.resize(n+1);
		int x, y, z = 0;
		for (unsigned long i = 0; i < m; ++i) {
			fscanf(fin, "%d%d%d", &x, &y, &z);
			Edge[x].emplace_back(y, z);
			Edge[y].emplace_back(x, z);
		}
	}

	void EdgeInitialize(){
		Edge.clear();
		for (int i = 0; i <= n; i++){
			vector<pair<int, int>> Ed;
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
};

struct Node{
	vector<int> VL, KVL, pv, pvd, kfvd, kpv, kpvd;
	vector<int> vert, kfv;
	int uniqueVertex;
	//, pos, pos2, dis
	//new code
	map<int,int> v_d;

	vector<int> ch;
	koala::my_openadd_hashmap<unsigned int> edges; // contains: vert & VL
	int height;
	int pa;
	Node(){
		vert.clear();
		VL.clear();
		KVL.clear();
		kfv.clear();
		kfvd.clear();
		kpv.clear();
		kpvd.clear();
		ch.clear();
		pv.clear();
		pvd.clear();
		v_d.clear();
		edges.clear();
		pa = -1;
		uniqueVertex = -1;
		height = 0;
	}
};

class Decomposition{
public: 
	FILE *fout, *fin;
	Graph G, H;
	int maxSize;
	int *curup_affect, curup_stamp;
	Decomposition(){
	}
	vector<int> ord;
	int heightMax;

	vector<int> vertexOrder, root_vertexes;
	vector<Node> tnodes;
	int vorder = 0;
	//int treewidth = 0;
	
	void reduction(){
		vector<vector<int>> degree2nodeQ; //degree->a list of vertex of this degree
		vector<int> ve;
		ve.clear();
		degree2nodeQ.push_back(ve);
		vector<pair<int,int>> vPosition(Edge.size()); //(degree,idx)
		int degree2nodeQ_size =0;
		for (int v = 1; v < Edge.size(); ++v) {
			int degree = Edge[v].size();
			if (degree >= degree2nodeQ.size()) {
				degree2nodeQ.resize(degree + 1);
			}
			vPosition[v] = make_pair(degree, degree2nodeQ[degree].size());
			degree2nodeQ[degree].push_back(v);
		}
		vector<koala::my_openadd_hashmap<unsigned int>> shortcuts(Edge.size());
		for (int s = 1; s < Edge.size(); ++s) {
			for (unsigned int i = 0; i < Edge[s].size(); ++i) {
				shortcuts[s].insert(Edge[s][i].first, Edge[s][i].second);
				//if(s==5) cout<<"("<<Edge[s][i].first<<","<< Edge[s][i].second<<")";
			}
		}
		//cout<<endl;

		vertexOrder.resize(Edge.size(), -1);
		ord.clear();
		int mindegree = 0;
		tnodes.resize(Edge.size());
		while (true) {
			int cnt =0;
			while (degree2nodeQ[mindegree].size() == 0) {
				cnt += degree2nodeQ[mindegree].size();
				mindegree++;
				if (mindegree == degree2nodeQ.size() && cnt==0 ) break;
			}
			if (mindegree == degree2nodeQ.size() && cnt==0) break;

			int vid = degree2nodeQ[mindegree].back();
			degree2nodeQ[mindegree].pop_back();
			
			ord.push_back(vid);
			vertexOrder[vid] = vorder++;
			//cout<<"vid: "<<vid<<" vertexOrder[vid]: "<<vertexOrder[vid]<<endl;
			
			koala::my_openadd_hashmap<unsigned int> &v = shortcuts[vid];
			vector<unsigned int> valid_neighbor_index;
			int cnt_width = 0;
			
			for (unsigned int i = v.iterator(); v.has_next(i); v.next(i)) {
				if (vertexOrder[v.get_with_idx(i).first] == -1) {
					valid_neighbor_index.push_back(i);
					int vec = v.get_with_idx(i).first, di = v.get_with_idx(i).second;
					tnodes[vid].edges.insert(vec, di);
					tnodes[vid].vert.push_back(vec);
					tnodes[vid].VL.push_back(di);
					//tnodes[vid].KVL.push_back(v.get_with_idx(i).second);
					//cout<<"("<<v.get_with_idx(i).first<<","<<v.get_with_idx(i).second<<") ";
					cnt_width++;
				}else{
					tnodes[vid].pv.push_back(v.get_with_idx(i).first);
					tnodes[vid].pvd.push_back(v.get_with_idx(i).second);
				}
			}//cout<<endl;
			//if (cnt_width>treewidth) treewidth = cnt_width;
			
			vector<int> neighbor_degree_increase_cnt(valid_neighbor_index.size(), -1);
			
			for (unsigned int ii = 0; ii < valid_neighbor_index.size(); ++ii) {
				for (unsigned int jj = ii + 1; jj < valid_neighbor_index.size(); ++jj) {
					unsigned int i = valid_neighbor_index[ii], j = valid_neighbor_index[jj];
					if (shortcuts[v.get_with_idx(i).first].find(v.get_with_idx(j).first) != -1) {//exist, update with min value
						shortcuts[v.get_with_idx(j).first][v.get_with_idx(i).first]
								= shortcuts[v.get_with_idx(i).first][v.get_with_idx(j).first]
								= min(shortcuts[v.get_with_idx(i).first][v.get_with_idx(j).first],
								      v.get_with_idx(i).second + v.get_with_idx(j).second);
					} else {//doesn't exist, add shortcut
						shortcuts[v.get_with_idx(j).first][v.get_with_idx(i).first] = shortcuts[v.get_with_idx(
								i).first][v.get_with_idx(j).first] =
								v.get_with_idx(i).second + v.get_with_idx(j).second;
						neighbor_degree_increase_cnt[ii]++;
						neighbor_degree_increase_cnt[jj]++;
					}
				}
			}
			
			for (unsigned int i = 0; i < valid_neighbor_index.size(); ++i) {
				if (neighbor_degree_increase_cnt[i] != 0) {
					unsigned int &vid = v.get_with_idx(valid_neighbor_index[i]).first;
					pair<int, int> &p = vPosition[vid];
					//swap and delete 
					if (degree2nodeQ[p.first][p.second] != degree2nodeQ[p.first].back()){
						degree2nodeQ[p.first][p.second] = degree2nodeQ[p.first].back();
						vPosition[degree2nodeQ[p.first].back()].second = p.second;
					}
					degree2nodeQ[p.first].pop_back();
					//place in a new position
					p.first += neighbor_degree_increase_cnt[i];
					if (p.first >= degree2nodeQ.size()) degree2nodeQ.resize(p.first + 1);
					p.second = degree2nodeQ[p.first].size();
					degree2nodeQ[p.first].push_back(vid);
					if (p.first < mindegree) mindegree = p.first;
				}
			}
			//cout<<"(edges, vert): ";
			//for (unsigned int i = 0; i < tnodes[vid].vert.size(); ++i) {
			//	tnodes[vid].edges.insert(tnodes[vid].vert[i],tnodes[vid].VL[i]);
			//	cout<<"("<<tnodes[vid].edges.get_with_idx(i).first<<", "<<tnodes[vid].vert[i]<<")";
			//}cout<<endl;
			/*
			cout<<"(2edges, vert): ";
			auto& e = tnodes[vid].edges;
			for (unsigned int i = e.iterator(); e.has_next(i); e.next(i)) {
			//for (unsigned int i = 0; i < tnodes[vid].vert.size(); i++) {
				cout<<"("<<i<<", "<<tnodes[vid].edges.get_with_idx(i).first<<", "<<tnodes[vid].vert[i]<<")";
			}cout<<endl;cout<<endl;*/
			/*
			if (tnodes[vid].vert.size() == 0) {
				root_vertexes.push_back(vid);
			}*/
		}

		/*
		int root_cnt =0;
		for (auto &l1:degree2nodeQ) {//collect remaining shortcuts, as root
			root_cnt += l1.size();
			for (int &vid:l1) {
				root_vertexes.push_back(vid);
				koala::my_openadd_hashmap<unsigned int> &v = shortcuts[vid];
			  	for (unsigned int i = v.iterator(); v.has_next(i); v.next(i)) {
					if (vertexOrder[v.get_with_idx(i).first]==-1) {
					  	tnodes[vid].edges.insert(v.get_with_idx(i));
					}
			    }
			}
		}
		//for(int i=0;i<root_vertexes.size();i++) cout<<root_vertexes[i]<<endl;
		cout<<"root_vertexes.size(): "<<root_vertexes.size()<<endl;
		cout<<"there are still "<< root_cnt <<" root vertexes not reduced yet "<< endl;*/
	}

	vector<Node> Tree;
	int root = 0;
	
	//只是在找parents 和 children 
	void pa_ch(){
		int cpin =0;
		int order = 0;
		int nearest =-1;
		int trlen = ord.size()-1;
		tnodes[ord[trlen]].height = 1;
		tnodes[ord[trlen]].pa = -1;
		root = ord[trlen];
		trlen--;

		for(; trlen >= 0; trlen--){
			int vi = ord[trlen];
			order = H.n;
			nearest = -1;
			int neighbor = 0;
			auto& edges = tnodes[vi].edges;
			for (unsigned int i = edges.iterator(); edges.has_next(i); edges.next(i)) {
				neighbor = edges.get_with_idx(i).first;
				if (vertexOrder[neighbor] < order){
					nearest = neighbor;
					order = vertexOrder[neighbor];
				}
			}
			int pa = nearest;
			//nod.pa = pa;
			if(cpin) cout<<"vi: "<<vi<<" pa: "<<pa<<endl;
			tnodes[vi].pa = pa;
			//ch = pv?
			tnodes[pa].ch.push_back(vi);
			//height 什么用 ？？除非用于统计L
			tnodes[vi].height = tnodes[pa].height + 1;
			if (tnodes[vi].height > heightMax) heightMax = tnodes[vi].height;
		}
	}

	inline unsigned int dist(int s, int t) {
		if (s == t) return 0;
		else if (vertexOrder[s] < vertexOrder[t])
			swap(s, t);
		return tnodes[t].edges[s];
	}

	void kvc(){	
		//int cpin = 0;
		int trlen = ord.size()-1;
		//tnodes[ord[trlen]].height = 1;
		//heightMax = 0;
		int g_n = trlen+1;
		int order = g_n;
		
		vector<bool> c_kvc;
		c_kvc.resize(g_n, true);
		int nearest=0, neighbor=0; 
		int fv =0;
		trlen--;
		//按照 order 在取vertex 
		for (; trlen >= 0; trlen--){
			int x = ord[trlen];
			//if(cpin) cout<<"x: "<<x<<endl;
			order = g_n;
			auto& edges = tnodes[x].edges;
			int ki =0, kj=0;
			for (unsigned int i = edges.iterator(); edges.has_next(i); edges.next(i)) {
				auto &neib1 = edges.get_with_idx(i);
				for (unsigned int j = i; edges.has_next(j); edges.next(j)) {
					auto &neib2 = edges.get_with_idx(j);
					unsigned int d = dist(neib1.first, neib2.first);
					if (neib1.second > (neib2.second + d)) {
						neib1.second = neib2.second + d;
						c_kvc[i] = false;
					}
					else if (neib2.second > (neib1.second + d)) {
						neib2.second = neib1.second + d;
						c_kvc[j] = false;
					}
				}
				neighbor = neib1.first;
				if (vertexOrder[neighbor] < order){
					nearest = neighbor;
					order = vertexOrder[neighbor];
				}
				if(c_kvc[i]){
				 	fv = neib1.first;
					tnodes[x].kfv.push_back(fv);
					tnodes[x].kfvd.push_back(neib1.second);
					tnodes[fv].kpv.push_back(x);
					tnodes[fv].kpvd.push_back(neib1.second);
				}else{
					c_kvc[i] = true;
				}
			}
		/*	
		if(cpin) {
			cout<<"nbr: ";
			for (unsigned int i = edges.iterator(); edges.has_next(i); edges.next(i)) {
			//for(int i = 0; i<tnode_size; i++){
				cout<<"("<<i<<", "<<edges.get_with_idx(i).first<<") ";
			}cout<<endl;
			
			cout<<"vert: ";
			for(int i = 0; i<tnodes[x].vert.size(); i++){
				cout<<tnodes[x].vert[i]<<" ";
			}cout<<endl;

			cout<<"dis: ";
			for(int i = 0; i<tnodes[x].VL.size(); i++){
				cout<<tnodes[x].VL[i]<<" ";
			}cout<<endl;

			cout<<"gdis: ";
			for (unsigned int i = edges.iterator(); edges.has_next(i); edges.next(i)) {
			//for(int i = 0; i<tnode_size; i++){
				cout<<edges.get_with_idx(i).second<<" ";
			}cout<<endl;
			
			cout<<"kfv: ";
			for(int i = 0; i<tnodes[x].kfv.size(); i++){
				cout<<tnodes[x].kfv[i]<<" ";
			}cout<<endl;
		}*/
			//int pa = nearest;
			//tnodes[x].pa = pa;
			//tnodes[pa].ch.push_back(x);
			//tnodes[x].height = tnodes[pa].height + 1;
			//if (tnodes[x].height > heightMax) heightMax = tnodes[x].height;
		}
	}

	void printIntVector(vector<int> &a){
		if (a.size() == 0){
			int x = 0;
			fwrite(&x, SIZEOFINT, 1, fout);
			return;
		}
		int x = a.size();
		fwrite(&x, SIZEOFINT, 1, fout);
		for (int i = 0; i < a.size(); i++){
			fwrite(&a[i], SIZEOFINT, 1, fout);
		}
	}
	
	void saveIndex(){

		int pin = 0;

		/*int cnt_fv = 0;
		int cnt_pv = 0;
		int cnt_kfv = 0;
		int cnt_kpv = 0;*/

		fwrite(&root, SIZEOFINT, 1, fout);
		Tree[root].vert.push_back(Tree[root].uniqueVertex);
		printIntVector(Tree[root].vert);
		Tree[root].vert.pop_back();
		printIntVector(Tree[root].KVL);
		printIntVector(Tree[root].pv);
		printIntVector(Tree[root].pvd);
		printIntVector(Tree[root].kpv);
		printIntVector(Tree[root].kpvd);			
		Tree[root].kfv.push_back(Tree[root].uniqueVertex);
		printIntVector(Tree[root].kfv);
		Tree[root].kfv.pop_back();
		printIntVector(Tree[root].kfvd);
		if(pin){
			cout<<"current v: "<<Tree[root].uniqueVertex<<endl;
			for(int i = 0; i<Tree[root].pv.size(); i++) cout<<"Tree[root].pv[i]: "<<Tree[root].pv[i]<<endl;
			for(int i = 0; i<Tree[root].kpv.size(); i++) cout<<"Tree[root].kpv[i]: "<<Tree[root].kpv[i]<<endl;
		}
		cnt = 0;

		/*cnt_fv = Tree[root].vert.size();
		cnt_pv = Tree[root].pv.size();
		cnt_kfv = Tree[root].kfv.size();
		cnt_kpv = Tree[root].kpv.size();*/

		//new code
		queue<int> Q;
		while (!Q.empty()) Q.pop();

		for (int i = 0; i < Tree[root].ch.size(); i++) Q.push(Tree[root].ch[i]);

		while(!Q.empty()){

			auto v = Q.front();
			Q.pop();

			/*cnt_fv += Tree[v].vert.size();
			cnt_pv += Tree[v].pv.size();
			cnt_kfv += Tree[v].kfv.size();
			cnt_kpv += Tree[v].kpv.size();*/

			for(auto a: Tree[v].ch){
				Q.push(a);
			}
            
			fwrite(&v, SIZEOFINT, 1, fout);

			/*if(pin){
				cout<<"current v: "<<Tree[v].uniqueVertex<<endl;
				for(int i = 0; i<Tree[v].pv.size(); i++) cout<<"Tree[v].pv[i]: "<<Tree[v].pv[i]<<endl;
				for(int i = 0; i<Tree[v].kpv.size(); i++) cout<<"Tree[v].kpv[i]: "<<Tree[v].kpv[i]<<endl;
				for(int i = 0; i<Tree[v].vert.size(); i++) cout<<"Tree[v].vert[i]: "<<Tree[v].vert[i]<<endl;
				for(int i = 0; i<Tree[v].kfv.size(); i++) cout<<"Tree[v].kfv[i]: "<<Tree[v].kfv[i]<<endl;
			}*/
			
		    Tree[v].vert.push_back(Tree[v].uniqueVertex);
		    printIntVector(Tree[v].vert);
		    Tree[v].vert.pop_back();
		    printIntVector(Tree[v].KVL);
			printIntVector(Tree[v].pv);
			printIntVector(Tree[v].pvd);
			printIntVector(Tree[v].kpv);
			printIntVector(Tree[v].kpvd);
			Tree[v].kfv.push_back(Tree[v].uniqueVertex);
			printIntVector(Tree[v].kfv);
			Tree[v].kfv.pop_back();
			printIntVector(Tree[v].kfvd);
			//printMapVector(Tree[v].v_d);
		}
		/*
		int pin_cnt=0;
		if(pin_cnt){
			cout<<"cnt_fv: "<<cnt_fv<<endl;
			cout<<"cnt_pv: "<<cnt_pv<<endl;
			cout<<"cnt_kfv: "<<cnt_kfv<<endl;
			cout<<"cnt_kpv: "<<cnt_kpv<<endl;
		}*/
	}
    int ug_size(int p){
		curup_stamp++;
		int ug_size = 0; 
		queue<int> cug; cug.push(p); ug_size++;
        int x = p;
		//kw*kh
		while(!cug.empty()){
			int p = cug.front(); cug.pop();
			for(int i=0;i<tnodes[p].kfvd.size();i++){
				if (curup_affect[tnodes[p].kfv[i]] != curup_stamp) {
					ug_size++; 
					cug.push(tnodes[p].kfv[i]);
					curup_affect[tnodes[p].kfv[i]] = curup_stamp;
				}
			}
		}
		
		return ug_size;
	}
	int max_ug_size(){
		int trlen = ord.size();
		curup_stamp = 0;
		curup_affect = (int*)malloc(sizeof(int) * (trlen+1));
		int max_ug_size = 0;
		for(int i = 0;i<trlen;i++){
			int us = ug_size(ord[i]);
			if(us > max_ug_size) max_ug_size = us;
		}
		return max_ug_size;
	}
	void cntSize(){
		int tree_size = ord.size();
		int kw = 0, w = 0;
		for( int i = 0; i < tree_size; ++i ) {

			int cur_w = tnodes[ord[i]].edges.size();
			if(cur_w > w) w = cur_w;
			
			int cur_kw = tnodes[ord[i]].kfv.size();
			if(cur_kw > kw) kw = cur_kw;
			
		}
		cout<<"w: "<<(w+1)<<endl;
		cout<<"kw: "<<(kw+1)<<endl;
		int kh = max_ug_size();
		cout<<"kh: "<<kh<<endl;
	}
	
	static const int SIZEOFINT = 4;
	void printIntArray(int * a, int n){
		fwrite(a, SIZEOFINT, n, fout);
	}
	
	void print(){
		//G.n
		fwrite(&G.n, SIZEOFINT, 1, fout); 
		printf("G.n %d\n", G.n);
		int x = Tree.size();
		fwrite(&x, SIZEOFINT, 1, fout);
		for (int i = 0; i < Tree.size(); i++){
			fwrite(&Tree[i].height, SIZEOFINT, 1, fout);
		}
		for (int i = 0; i < Tree.size(); i++){
			fwrite(&Tree[i].pa, SIZEOFINT, 1, fout);
		}
		for (int i = 0; i < Tree.size(); i++){
			fwrite(&Tree[i].uniqueVertex, SIZEOFINT, 1, fout);
		}
		//belong
		//printIntArray(belong, H.n + 1);
        int pt=0;
		//rootDistance
		fwrite(&root, SIZEOFINT, 1, fout);
		for (int i = 0; i < Tree.size(); i++){
			if(pt) cout<<"ch[]"<<Tree[i].uniqueVertex<<": ";
			int t = Tree[i].ch.size();
			fwrite(&t, SIZEOFINT, 1, fout);
			for (int j = 0; j < t; j++){
				fwrite(&Tree[i].ch[j], SIZEOFINT, 1, fout);
				if(pt) cout<<Tree[i].ch[j]<<" ";
			}
			if(pt) cout<<endl;
		}	
	}
};

// floyd deal with vertices which is deleted
int main(int argc, char *argv[])
{
	srand(time(0));
	int operation = 1;
	if (operation == 1){ // index 
		string filest;
		char *file, *fileout;
		//int i;
		file = argv[1];
		cout << "file: " << file << endl;
		fileout = argv[2];
		Decomposition td;
		td.fout = fopen(fileout, "wb");
		// read graph data
		td.G = Graph(file);
		
		td.H = td.G;
		cout<<"start reduction"<<endl;
		clock_t start = clock();
		td.reduction();
		double redu_time =  (double)(clock() - start)/ CLOCKS_PER_SEC;
		cout<<"finish reduce: "<< redu_time <<endl;

		//clock_t pc_begin = clock();
		//td.pa_ch();
		//cout << "pa_ch time: " << (double)(clock() - pc_begin) / CLOCKS_PER_SEC << endl; 
		
		clock_t kvc_begin = clock();
		td.kvc();
		double kvc_time =  (double)(clock() - kvc_begin)/ CLOCKS_PER_SEC;
		cout << "KVC: " << kvc_time << endl; 
		cout<<"all ch time: "<< (kvc_time + redu_time) <<endl;
		//cout << "MakeIndex+reduce time: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		//td.cntSize();

		//td.print();
        //cout << "Tree height: " << td.heightMax << endl;

		fclose(stdout);
	}

}
