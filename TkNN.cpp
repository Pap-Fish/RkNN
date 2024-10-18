#include <cstdio>
#include <cstring>
#include <iostream>
#include <cstdlib>
#include <queue>
#include <sys/time.h>
#include <vector>
#include <xmmintrin.h>
#include <cmath>
#include <set>
#include <algorithm>
#include <fstream>
#include "graph.h"
#include "TkNN.h"
using namespace std;
#define MAX_K 100
#define INT_MAX 999999999
#define RESERVE_TIME 1

int reset_times = 0;

const int infinity = 999999999;
const int SIZEOFINT = 4;

double cnt_pre_query_time = 0;

int *toRMQ, *height, *pa, *uniqueVertex, **RMQIndex;
int *belong;
int root, TreeSize;
int **rootToRoot, *rootSite;
int **dis, **pos, **pos2;
int *posSize, *pos2Size;
int *chSize;
int **ch;
int *LOG2, *LOGD;
int rootSize;
int *DFSList, *toDFS;
int ***BS;
int *tmp_dis, *higher, *cloest_higher;

int **SPN;
// vector<map<int,int>> SPN;
int LCAQuery(int _p, int _q)
{
	int p = toRMQ[_p], q = toRMQ[_q];
	//	printf("p q : %d %d\n", p, q);
	if (p > q)
	{
		int x = p;
		p = q;
		q = x;
	}
	int len = q - p + 1;
	//	printf("len: %d\n", len);
	int i = LOGD[len], k = LOG2[len];
	//	printf("i, k: %d %d\n", i, k);
	q = q - i + 1;
	if (height[RMQIndex[k][p]] < height[RMQIndex[k][q]])
		return RMQIndex[k][p];
	else
		return RMQIndex[k][q];
}

long long queryCnt;
// long long aCnt;

int distanceQuery(int p, int q)
{
	if (p == q)
		return 0;
	int x = belong[p], y = belong[q];

	if (toRMQ[x] > toRMQ[y])
	{
		int k = x;
		x = y;
		y = k;
	}
	//	return dis[y][pos[x][posSize[x] - 1]];
	//	if (p == 9817 && q == 9816)
	//		printf("HHHHHHH x, y, height[x], height[y]: %d %d %d %d\n", x, y, height[x], height[y]);
	return dis[y][height[x] - 1];
}


int distanceQueryFull(int p, int q){
	if (p == q) return 0;
	int x = belong[p], y = belong[q];

	int lca = LCAQuery(x, y);

	if (lca == x || lca == y){
		queryCnt++;
		if (lca == y){
			int v = y;
			y = x;
			x = v;
			v = p;
			p = q;
			q = v;
		}
		return dis[y][height[x]-1];
	}
	else {
		int res = infinity;
		int *dx = dis[x], *dy = dis[y],*p2 = pos2[lca];
		_mm_prefetch(dx, _MM_HINT_T0);
		_mm_prefetch(dy, _MM_HINT_T0);
		_mm_prefetch(p2, _MM_HINT_T0);
		int ps = pos2Size[lca];
		for (int i = 0; i < ps; i++){
			queryCnt ++;
			int tmp = dx[p2[i]] + dy[p2[i]];
			if (res > tmp)
				res = tmp;
		}
		return res;
	}

}
// extern int *Degree;
// extern int **Neighbor, **Weight;
// int SPNQuery(int x, int target){
// 	if(x == target)
// 		return x;
// 	int p = belong[x], q = belong[target];
// 	int lca = LCAQuery(p, q);
// 	if(lca == q){
// 		// if(lca == p){
// 		// 	int v = p;
// 		// 	p = q;
// 		// 	q = v;
// 		// }
// 		return SPN[p][height[q] - 1];
// 	}
// 	else if(lca == p){
// 		int dist = distanceQueryFull(x,target);
// 		for(int i=0; i<Degree[x]; i++){
// 			int nb = Neighbor[x][i];
// 			int we = Weight[x][i];
// 			int _d = distanceQueryFull(nb,target);
// 			if(_d+we == dist){
// 				return nb;
// 			}
// 		}
// 	}
// 	else{
// 		int sp_dist = infinity;
// 		int *dx = dis[p], *dy = dis[q], *p2 = pos2[lca];
// 		_mm_prefetch(dx, _MM_HINT_T0);
// 		_mm_prefetch(dy, _MM_HINT_T0);
// 		_mm_prefetch(p2, _MM_HINT_T0);
// 		int ps = pos2Size[lca];
// 		int anc_pos = 0;
// 		for(int i=0; i<ps; i++){
// 			int tmp = dx[p2[i]] + dy[p2[i]];
// 			if(sp_dist > tmp){
// 				sp_dist = tmp;
// 				anc_pos = p2[i];
// 			}
// 		}
// 		// int t_anc;
// 		// for(int i=0; i<posSize[lca]; i++){
// 		// 	int anc = pos[lca][i];
// 		// 	if(height[belong[anc]]-1 == anc_pos){
// 		// 		t_anc = anc;
// 		// 		break;
// 		// 	}
// 		// }
// 		return SPN[p][anc_pos];
// 	}
// }

// int *Degree;
// int **Neighbor, **Weight;
// void readGraph(char *filename)
// {
// 	FILE *file = fopen(filename, "r");
// 	int n, m;
// 	fscanf(file, "%d %d", &n, &m);
// 	Degree = (int *)malloc(sizeof(int) * (n + 1));
// 	vector<vector<pair<int, int>>> nb;
// 	vector<pair<int, int>> v;
// 	v.clear();
// 	for (int i = 0; i <= n; i++)
// 	{
// 		//	Degree[i] = 0;
// 		nb.push_back(v);
// 	}
// 	//	cout << n << " " << m << endl;
// 	for (int i = 0; i < m; i++)
// 	{
// 		int x, y, z;
// 		fscanf(file, "%d %d %d", &x, &y, &z);
// 		//		Degree[x]++;
// 		//		cout << x << " " << y << " " << z << endl;
// 		nb[x].push_back(make_pair(y, z));
// 	}
// 	Neighbor = (int **)malloc(sizeof(int *) * (n + 1));
// 	Weight = (int **)malloc(sizeof(int *) * (n + 1));
// 	for (int i = 1; i <= n; i++)
// 	{
// 		Degree[i] = nb[i].size();
// 		Neighbor[i] = (int *)malloc(sizeof(int) * nb[i].size());
// 		Weight[i] = (int *)malloc(sizeof(int) * nb[i].size());
// 		for (int j = 0; j < nb[i].size(); j++)
// 		{
// 			Neighbor[i][j] = nb[i][j].first;
// 			Weight[i][j] = nb[i][j].second;
// 		}
// 	}
// }
// inline int shortestPathQuery(int p, int q)
// {
// 	int res = 0;
// 	while (p != q)
// 	{
// 		res++;
// 		int pq = distanceQuery(p, q);
// 		for (int i = 0; i < Degree[p]; i++)
// 		{
// 			int x = Neighbor[p][i];
// 			//	int y = Weight[p][i];
// 			int xq = distanceQuery(x, q);
// 			if (xq + Weight[p][i] == pq)
// 			{
// 				p = x;
// 				break;
// 			}
// 		}
// 	}
// 	return res;
// }

int test(char *file)
{
	//	cout << "test: " << file << " BEGIN" << endl;
	FILE *fin = fopen(file, "r");
	int x, y, dis, res = 0;
	vector<int> X, Y, DIS;
	X.clear();
	Y.clear();
	DIS.clear();
	while (fscanf(fin, "%d %d %d", &x, &y, &dis) != EOF)
	{
		if (x <= 0 || y <= 0)
			break;
		X.push_back(x);
		Y.push_back(y);
		DIS.push_back(dis);
		res++;
	}
	cout << X.size() << endl;

	//
	queryCnt = 0;
	//	aCnt = 0;
	//
	double t = 0;
	int *ANS;
	ANS = (int *)malloc(sizeof(int) * (X.size() + 1));
	double start_time = GetTime();
	for (int i = 0; i < X.size(); i++)
	{
		ANS[i] = distanceQuery(X[i], Y[i]);
		//	if (ANS[i] != DIS[i]){
		//		cout << X[i] << " " << Y[i] << " " << DIS[i] << " " << ANS[i] << endl;
		//		while(1);
		//	}
	}
	double end_time = GetTime();
	for (int i = 0; i < X.size(); i++)
	{
		t += ANS[i];
	}
	//	cout << end_time - start_time << endl;
	//	cout << res << endl;
	cout << "Check Count: " << double(queryCnt) / res << endl;
	//	cout << "aCnt: " << double(aCnt) / res << endl;
	printf("Distance Query Time : %lf usec\n", (end_time - start_time) * 1e6 / res);
	printf("average distance: %.6lf\n", t / res);
	/*
	start_time = GetTime();
	long long step = 0;
	for (int i = 0; i < X.size(); i++){
		step += shortestPathQuery(X[i], Y[i]);
	}
	end_time = GetTime();
	cout << step;
	printf("Shortest Path Query Time : %lf usec\n", (end_time - start_time) * 1e6 / res);
	printf("average step: %.6lf\n", double(step) / double(res));
	*/
	return res;
}

FILE *fin;
string TT = "";
void scanIntArray(int *a, int n)
{
	fread(a, SIZEOFINT, n, fin);
}
int *scanIntVector(int *a)
{
	int _n;
	fread(&_n, SIZEOFINT, 1, fin);
	a = (int *)malloc(sizeof(int) * _n);
	scanIntArray(a, _n);
	return a;
}

int n;
int *EulerSeq;
void readIndex(char *file)
{
	double _time = GetTime();
	int tree_height = 0, tree_width = 0, most_sp = 0;
	fin = fopen(file, "rb");
	fread(&n, SIZEOFINT, 1, fin);
	printf("n: %d\n", n);
	int ts;
	fread(&ts, SIZEOFINT, 1, fin);
	printf("ts: %d\n", ts);
	TreeSize = ts;
	height = (int *)malloc(sizeof(int) * (ts + 1));
	pa = (int *)malloc(sizeof(int) * (ts + 1));
	uniqueVertex = (int *)malloc(sizeof(int) * (ts + 1));
	for (int i = 0; i < ts; i++)
	{
		fread(&height[i], SIZEOFINT, 1, fin);
	}
	for (int i = 0; i < ts; i++)
	{
		fread(&pa[i], SIZEOFINT, 1, fin);
	}
	for (int i = 0; i < ts; i++)
	{
		fread(&uniqueVertex[i], SIZEOFINT, 1, fin);
	}
	belong = (int *)malloc(sizeof(int) * (n + 1));
	fread(belong, SIZEOFINT, n + 1, fin);
	toRMQ = (int *)malloc(sizeof(int) * (n + 1));
	fread(toRMQ, SIZEOFINT, n + 1, fin);
	int ris;
	fread(&ris, SIZEOFINT, 1, fin);
	//	printf("ris: %d\n", ris);
	fread(&ts, SIZEOFINT, 1, fin);
	//	printf("ts: %d\n", ts);
	EulerSeq = (int *)malloc(sizeof(int) * (ts + 1));
	RMQIndex = (int **)malloc(sizeof(int *) * (ris + 1));
	for (int i = 0; i < ris; i++)
	{
		RMQIndex[i] = scanIntVector(RMQIndex[i]);
	}
	fread(&root, SIZEOFINT, 1, fin);
	//	cout << "root: " << root << endl;

	posSize = (int *)malloc(sizeof(int) * (n + 1));
	pos2Size = (int *)malloc(sizeof(int) * (n + 1));
	pos = (int **)malloc(sizeof(int *) * (TreeSize));
	pos2 = (int **)malloc(sizeof(int *) * (TreeSize));
	dis = (int **)malloc(sizeof(int *) * (TreeSize));
	chSize = (int *)malloc(sizeof(int) * (TreeSize));
	ch = (int **)malloc(sizeof(int *) * (TreeSize));
	// SPN = (int **)malloc(sizeof(int *) * (TreeSize));

	for (int i = 0; i < TreeSize; i++)
	{
		fread(&chSize[i], SIZEOFINT, 1, fin);
		ch[i] = (int *)malloc(sizeof(int) * chSize[i]);
		for (int j = 0; j < chSize[i]; j++)
		{
			int x;
			fread(&x, SIZEOFINT, 1, fin);
			ch[i][j] = x;
		}
	}
	for (int i = 0; i < TreeSize; i++)
	{
		int x;
		fread(&x, SIZEOFINT, 1, fin);
		fread(&posSize[x], SIZEOFINT, 1, fin);
		pos[x] = (int *)malloc(sizeof(int) * (posSize[x] + 1));
		fread(pos[x], SIZEOFINT, posSize[x], fin);
		if (posSize[x] > tree_width)
			tree_width = posSize[x];
		int _n;
		fread(&_n, SIZEOFINT, 1, fin);
		dis[x] = (int *)malloc(sizeof(int) * _n);
		fread(dis[x], SIZEOFINT, _n, fin);
		if (_n > tree_height)
			tree_height = _n;
		// int _s;
		// fread(&_s, SIZEOFINT, 1, fin);
		// SPN[x] = (int *)malloc(sizeof(int) * _s);
		// fread(SPN[x], SIZEOFINT, _s, fin);
	}
	//	printf("dis read finished!\n");
	for (int i = 0; i < TreeSize; i++)
	{
		int x;
		fread(&x, SIZEOFINT, 1, fin);
		fread(&pos2Size[x], SIZEOFINT, 1, fin);
		pos2[x] = (int *)malloc(sizeof(int) * (pos2Size[x] + 1));
		fread(pos2[x], SIZEOFINT, pos2Size[x], fin);
		if (pos2Size[x] > most_sp)
			most_sp = pos2Size[x];
	}

	fclose(fin);
	higher = (int *)malloc(sizeof(int) * (n + 1));
	cloest_higher = (int *)malloc(sizeof(int) * (n + 1));
	for (int i = 0; i < n; i++)
	{
		higher[i] = n;
		cloest_higher[i] = INT_MAX;
		for (int j = 0; j < posSize[i]; j++)
		{
			if (height[belong[pos[i][j]]] < higher[i])
				higher[i] = height[belong[pos[i][j]]];
			int _d = distanceQuery(uniqueVertex[i], pos[i][j]);
			if (cloest_higher[i] > _d)
				cloest_higher[i] = _d;
		}
	}
	LOG2 = (int *)malloc(sizeof(int) * (n * 2 + 10));
	LOGD = (int *)malloc(sizeof(int) * (n * 2 + 10));
	int k = 0, j = 1;
	for (int i = 0; i < n * 2 + 10; i++)
	{
		if (i > j * 2)
		{
			j *= 2;
			k++;
		}
		LOG2[i] = k;
		LOGD[i] = j;
	}
	// printf("Load Index Time : %lf sec\n", (GetTime() - _time));
	// printf("tree height: %d\n", tree_height);
	// printf("tree width: %d\n", tree_width);
	// printf("most search space: %d\n", most_sp);
}
int cnt;
void getDFSListDFS(int p)
{
	toDFS[p] = cnt;
	DFSList[cnt++] = p;
	for (int i = 0; i < chSize[p]; i++)
	{
		getDFSListDFS(ch[p][i]);
	}
	BS[p] = (int **)malloc(sizeof(int *) * chSize[p]);
	for (int i = 0; i < chSize[p]; i++)
	{
		BS[p][i] = (int *)malloc(sizeof(int) * chSize[p]);
		for (int j = 0; j < chSize[p]; j++)
		{
			if (posSize[ch[p][i]] < posSize[ch[p][j]])
				BS[p][i][j] = ch[p][i];
			else
				BS[p][i][j] = ch[p][j];
		}
	}
}
void getDFSList()
{
	DFSList = (int *)malloc(sizeof(int) * (TreeSize + 1));
	toDFS = (int *)malloc(sizeof(int) * (TreeSize + 1));
	BS = (int ***)malloc(sizeof(int **) * (TreeSize + 1));
	cnt = 0;
	getDFSListDFS(root);
}

bool anc_compare(int p, int q)
{
	if (tmp_dis[p] < tmp_dis[q])
		return true;
	return false;
}

void stop()
{
	while (1)
		;
}

int *is_current_object,*is_current_user,*current_distance, *group_height, *current_state, current_stamp;

static bool object_compare(int a, int b)
{
	if (current_distance[a] < current_distance[b])
		return true;
	else
		return false;
}

void TkNN::initialize_knn(char *index_file, char *obj_file)
{
	srand((int)(time(0)));
	readIndex(index_file);

	create_kNN_index();
	FILE *fobj = fopen(obj_file, "r");
	int num_obj, num_user;
	fscanf(fobj, "%d %d\n", &num_obj, &num_user);
	double start_time = GetTime();
	object_setting(n);

	for (int i = 0; i <= n; i++){
		is_current_object[i] = 0;
		is_current_user[i] = 0;
	}
	for (int i = 0; i < num_obj + num_user; i++)
	{
		int x;
		char c;
		fscanf(fobj, "%c %d\n", &c, &x);
		if(c == 'q')
			is_current_object[x] = 1;
		else if(c == 'u')
			is_current_user[x] = 1;
	}
	int cnt_object = 0, cnt_user = 0;
    for (int i = 1; i <= n; i++){
        if (is_current_object[i] == 1)
            cnt_object++;
        if (is_current_user[i] == 1)
            cnt_user++;
    }
	printf("cnt_object: %d, cnt_user: %d\n", cnt_object, cnt_user);
	double bef_init = GetTime();
	initialize_object();

	printf("knn.object_number[root]: %d\n", object_number[root]);
	double end_time = GetTime();
	printf("index initialization time: %.6lf ms\n", (end_time - start_time) * 1e3);

	//----------------------query------------------------------
	if (check_everyone())
	{
		printf("right\n");
	}
	else
		printf("wrong\n");
	for (int i = 0; i < period; i++)
		times_period[i].clear();
	tmp_dis = (int *)malloc(sizeof(int) * (n + 1));
	query_mark = (int *)malloc(sizeof(int) * (n + 1));
	query_mark_stamp = 0;
	for (int i = 0; i <= n; i++)
		query_mark[i] = 0;
}

void TkNN::create_kNN_index()
{
	// /	object_set.clear();
	vector<pair<int, int>> _v;
	_v.clear();
	for (int i = 0; i <= TreeSize; i++)
		path_from_root.push_back(_v);
	OSS.clear();
	for (int i = 0; i < TreeSize; i++)
	{
		object_saveing_structure oss;
		OSS.push_back(oss);
	}
	printf("TreeSize: %d\n", TreeSize);
}
void TkNN::delete_element(int p, int x)
{
	OSS[p].size_num--;
	int pre = OSS[p].a[x].previous;
	int ne = OSS[p].a[x].next;
	OSS[p].a[ne].previous = pre;
	OSS[p].a[pre].next = ne;
}
inline void TkNN::insert(int p, int x)
{
	//	printf("p, x: %d %d\n", p, x);
	//  if (x == 30137)
	//       cout << "(--" << p << ", " << x << "--)";

	int y = uniqueVertex[p];
	//	if (p == 59532)
	//		printf("p, y, x: %d %d %d\n", p, y, x);
	//	printf("y: %d\n", y);

	//    if (p == 1476 && x == 30137)
	//        cout << "heihei" << endl;
	int disxy = distanceQuery(y, x);
	if (OSS[p].size_num >= int(MAX_K))
		if (OSS[p].a[OSS[p].a[0].previous].dist <= disxy)
			return;
	//	printf("y, x, distanceQuery(y, x): %d %d %d\n", y, x, disxy);
	//	if (y == 9817)
	//		printf("y, disxy: %d %d\n", y, disxy);
	int i = 0;
	//	printf("OSS.size(): %d\n", OSS.size());
	//	printf("OSS[p].a.size(): %d\n", OSS[p].a.size());
	//	printf("(*a)[0].key (*a)[0].dist, (*a)[0].next, (*a)[0].previous: %d %d %d %d\n", OSS[p].a[0].key, OSS[p].a[0].dist, OSS[p].a[0].next, OSS[p].a[0].previous);
	while (OSS[p].a[i].previous != 0 && OSS[p].a[OSS[p].a[i].previous].dist > disxy)
		i = OSS[p].a[i].previous;
	list_node _a;
	_a.next = i;
	_a.previous = OSS[p].a[i].previous;
	_a.key = x;
	_a.dist = disxy;

	OSS[p].a.push_back(_a);
	//
	hash.insert_node(p, x, OSS[p].a.size() - 1);
	//
	OSS[p].a[_a.previous].next = OSS[p].a.size() - 1;
	OSS[p].a[_a.next].previous = OSS[p].a.size() - 1;
	OSS[p].size_num++;
	if (OSS[p].size_num > MAX_K * RESERVE_TIME)
		delete_element(p, OSS[p].a[0].previous);
}
inline void TkNN::insert(int x)
{
	//	printf("insert x: %d\n", x);
	//	if (object_set.find(x) != object_set.end())
	//		return;
	//	object_set.insert(x);
	if (is_current_object[x] == 1)
		return;
	is_current_object[x] = 1;

	int p = belong[x];
	//	vector<pair<int, int> > list;
	//	list.clear();
	while (p >= 0)
	{
		//			printf("insert p, x: %d %d\n", p, x);
		object_number[p]++;
		insert(p, x);
		//		list.push_back(make_pair(p, OSS[p].a.size() - 1));
		p = pa[p];
	}
	//	printf("list.size(): %d\n", list.size());
	//	for (int i = list.size() - 1; i >= 0; i--)
	//		path_from_root[x].push_back(list[i]);
	//	printf("x, path_from_root[x].size(): %d %d\n", x, path_from_root[x].size());
}
void TkNN::get_all_object(int p, vector<int> &a)
{
	int x = uniqueVertex[p];
	if (is_current_object[x])
	{
		a.push_back(x);
	}
	for (int i = 0; i < chSize[p]; i++)
	{
		get_all_object(ch[p][i], a);
	}
}

void TkNN::OSS_push_back(int p, int key, int dist)
{
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

	//	printf("p(%d):", p);
	//	for (int i = OSS[p].a[0].next; i != 0; i = OSS[p].a[i].next)
	//		printf("(%d, %d)", OSS[p].a[i].key, OSS[p].a[i].dist);
	//	printf("\n");
}
void TkNN::OSS_push_front(int p, int key, int dist)
{
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
void TkNN::get_subtree(int p, int key, vector<int> &a)
{
	//	printf("p, key, a.size(): %d %d %d\n", p, key, a.size());
	bool is_containing = false;
	int x = uniqueVertex[p];
	if (x == key)
		is_containing = true;
	else
	{
		//		printf("posSize[p]: %d\n", posSize[p]);
		for (int i = 0; i < posSize[p]; i++)
		{
			//			printf("pos[p][i]: %d\n", pos[p][i]);
			if (pos[p][i] == key)
			{
				is_containing = true;
				break;
			}
		}
	}
	//	printf("is_containing: %d\n", is_containing);
	if (is_containing)
	{
		a.push_back(p);

		//		printf("a.size(), chSize[p]: %d %d\n", a.size(), chSize[p]);
		for (int i = 0; i < chSize[p]; i++)
			get_subtree(ch[p][i], key, a);
	}
}
void TkNN::dfs_sort(int p, vector<int> &a)
{
	int _MAX = int(MAX_K * RESERVE_TIME);
	vector<int> b;
	a.clear();
	int x = uniqueVertex[p];
	if (is_current_object[x] == 1)
		a.push_back(x);
	for (int i = 0; i < chSize[p]; i++)
	{
		dfs_sort(ch[p][i], b);
		a.insert(a.end(), b.begin(), b.end());
	}
	int current_k;
	for (int i = 0; i < a.size(); i++)
		current_distance[a[i]] = distanceQuery(x, a[i]);
	if (_MAX < a.size())
	{
		nth_element(a.begin(), a.begin() + _MAX, a.end(), object_compare);
		current_k = _MAX;
	}
	else
	{
		current_k = a.size();
	}
	//	printf("a.size(), current_k: %d %d\n", a.size(), current_k);
	// sort(a.begin(), a.begin() + current_k, object_compare);
	sort(a.begin(), a.end(), object_compare);

	//	if (OSS[p].size_num > MAX_K)
	//		delete_element(p, OSS[p].a[0].previous);
	if (p == root)
	{
		//		for (int i = 0; i < a.size(); i++)
		//			printf("%d ", a[i]);
		//		printf("\n");
	}
	for (int i = 0; i < current_k; i++)
	{
		OSS_push_back(p, a[i], current_distance[a[i]]);
	}
	if (!double_objects(p))
	{
		print(p);
		stop();
	}
	//	printf("p, a.size()):%d %d\n", p, a.size());
}

void TkNN::join_subtree(int p, int x, vector<int> &a)
{

	priority_queue<PT> Q;
	while (!Q.empty())
		Q.pop();
	vector<int> iter;
	iter.clear();
	iter.push_back(0);
	current_stamp++;
	for (int i = 1; i < a.size(); i++)
	{
		iter.push_back(OSS[a[i]].a[0].next);
		int t;
		while (iter[i] != 0)
		{
			t = OSS[a[i]].a[iter[i]].key;
			//	cout << "a[i]: " << a[i] << ", iter[i]: " << iter[i] << endl;
			//	cout << "t: " << t << " "  << current_state[t] << " " << current_stamp << " " << current_distance[t] << " " << distanceQuery(t, x) << endl;;
			if (current_state[t] != current_stamp)
			{
				current_state[t] = current_stamp;
				current_distance[t] = distanceQuery(t, x);
				break;
			}
			else
				iter[i] = OSS[a[i]].a[iter[i]].next;
		}
		if (iter[i] != 0)
		{
			Q.push(PT(current_distance[t], i));
			//		cout << "push: " << p << ", " << a[i] << ", " << current_distance[t] << ", " << OSS[a[i]].a[iter[i]].dist << endl;q4
		}
	}
	vector<int> b;
	b.clear();
	for (int i = 0; i < int(RESERVE_TIME * MAX_K); i++)
	{
		if (!Q.empty())
		{
			PT pt = Q.top();
			Q.pop();
			//	OSS_push_back(p, OSS[a[pt.x]].a[iter[pt.x]].key, pt.dis);
			b.push_back(OSS[a[pt.x]].a[iter[pt.x]].key);
			//		cout << "pop: " << p << ", " << OSS[a[pt.x]].a[iter[pt.x]].key << ", " << pt.dis << endl;
			int t;
			while (iter[pt.x] != 0)
			{
				t = OSS[a[pt.x]].a[iter[pt.x]].key;
				//			cout << "a[pt.x]: " << a[pt.x] << ", iter[pt.x]: " << iter[pt.x] << endl;
				//			cout << current_state[t] << " " << current_stamp << " " << current_distance[t] << " " << distanceQuery(t, x) << endl;
				if (current_state[t] != current_stamp)
				{
					current_state[t] = current_stamp;
					current_distance[t] = distanceQuery(t, x);
					break;
				}
				else
					iter[pt.x] = OSS[a[pt.x]].a[iter[pt.x]].next;
			}
			if (iter[pt.x] != 0)
			{
				Q.push(PT(current_distance[t], pt.x));
				//		cout << "push: " << p << ", " << a[pt.x] << ", " << current_distance[t]<< ", " << OSS[a[pt.x]].a[iter[pt.x]].dist << endl;

				//		if (OSS[a[pt.x]].a[iter[pt.x]].key == 258384){
				//		cout << a[pt.x] << " " << uniqueVertex[a[pt.x]] << endl;
				//		}
			}
			//	else cout << pt.x << " " << a[pt.x] << " end " << endl;
		}
		else
			break;
	}
	//	cout << endl;
	sort(b.begin(), b.end(), object_compare);
	for (int i = 0; i < b.size(); i++)
	{
		OSS_push_back(p, b[i], current_distance[b[i]]);
		hash.insert_node(p, b[i], OSS[p].a.size() - 1);
	}
}
void TkNN::dfs_neighbor(int p)
{
	//    printf("dfs_neighbor: %d\n", p);
	vector<int> a;
	a.clear();
	int x = uniqueVertex[p];
	if (is_current_object[x] == 1)
	{
		//           printf("%d is object\n", x);
		OSS_push_front(p, x, 0);
	}
	for (int i = 0; i < chSize[p]; i++)
	{
		dfs_neighbor(ch[p][i]);
	}
	get_subtree(p, x, a);
	//	cout <<"sub" << uniqueVertex[p] << ": ";
	//	for (int i = 0; i < a.size(); i++)
	//		cout << uniqueVertex[a[i]] << ", ";
	//	cout << endl;
	join_subtree(p, x, a);

	if (!double_objects(p))
	{
		print(p);
		stop();
	}
}
void TkNN::traversal(int p)
{
	//    printf("p uniqueVertex[p] %d %d: ", p, uniqueVertex[p]);
	//        for (int i = 0; i < posSize[p]; i++)
	//        printf(" %d", pos[p][i]);
	//    printf("\n");
	int x = uniqueVertex[p];
	if (is_current_object[x] == 1)
		object_number[p]++;
	if (is_current_user[x] == 1)
		user_number[p]++;
	
	for (int i = 0; i < chSize[p]; i++)
	{
		traversal(ch[p][i]);
		object_number[p] += object_number[ch[p][i]];
		user_number[p] += user_number[ch[p][i]];
	}
}
void TkNN::compute_object_number()
{
	object_number.resize(n + 1);
	user_number.resize(n + 1);
	for (int i = 0; i <= n; i++){
		object_number[i] = 0;
		user_number[i] = 0;
	}
	traversal(root);
}
void TkNN::initialize_object()
{
	compute_object_number();
	printf("start initialize:\n");
	//	traversal(root);
	// type 1 sort
	vector<int> a;
	for (int i = 0; i < n; i++)
		for (int j = OSS[i].a[0].next; j != 0; j = OSS[i].a[j].next)
			delete_element(i, j);
	//
	//	dfs_merge(root, a);
	// if (strcmp(insert_type, "sort") == 0)
	// 	dfs_sort(root, a);
	// else
	dfs_neighbor(root);
	cout << "initialize object finished" << endl;
}

vector<pair<int, int>> TkNN::query(int x, int top_k, int limit_dis)
{
	query_mark_stamp++;
	int cnt1 = 0, cnt2 = 0, cnt3 = 0;
	vector<pair<int, int>> result;
	result.clear();
	int p = belong[x];
	vector<int> anc;
	anc.clear();
	int MAX_DIS = INT_MAX;
	//	return result;
	//
	vector<int> a;
	a.clear();
	a.push_back(p);
	for (int i = 0; i < posSize[p]; i++)
		a.push_back(belong[pos[p][i]]);
	//	printf("a.size(): %d\n", a.size());
	//	printf("MAX_DIS: %d\n", MAX_DIS);
	for (int i = 0; i < a.size(); i++)
	{
		int q = a[i];
		if (OSS[q].size_num <= top_k)
			continue;
		tmp_dis[q] = distanceQuery(uniqueVertex[q], x);
		if (tmp_dis[q] >= MAX_DIS)
			continue;
		if (OSS[q].a[0].next == 0)
			continue;
		int j = OSS[q].a[0].next;
		int _cnt = 1;
		while (_cnt < top_k)
		{
			_cnt++;
			j = OSS[q].a[j].next;
			cnt2++;
		}
		if (MAX_DIS > tmp_dis[q] + OSS[q].a[j].dist)
			MAX_DIS = tmp_dis[q] + OSS[q].a[j].dist;
	}

	int max_height = height[p];
	while (p >= 0 && height[p] >= max_height)
	{
		cnt3++;
		int t = OSS[p].a[0].next;
		tmp_dis[p] = distanceQuery(uniqueVertex[p], x);
		if (t > 0 && MAX_DIS >= tmp_dis[p] + OSS[p].a[t].dist)
		{
			anc.push_back(p);
			OSS[p].current = OSS[p].a[0].next;
			if (OSS[p].size_num >= top_k)
			{
				int v = OSS[p].a[0].previous;
				int tmp = tmp_dis[p] + OSS[p].a[v].dist;
				if (tmp < MAX_DIS)
					MAX_DIS = tmp;
			}
		}
		if (tmp_dis[p] + cloest_higher[p] < MAX_DIS)
			if (max_height > higher[p])
			{
				max_height = higher[p];
			}
		p = pa[p];
	}
	sort(anc.begin(), anc.end(), anc_compare);
	//	printf("anc.size(): %d\n", anc.size());
	p = belong[x];
	int _cnt = 0;
	for (int i = 0; i < top_k; i++)
	{
		//	printf("i result.size(): %d %d\n", i, result.size());
		//--------------find minimum-------------
		int k = -1, dist_k = -1;
		for (int j = 0; j < anc.size(); j++)
		{
			_cnt++;
			int q = anc[j];
			//	printf("q: %d\n", q);
			while (OSS[q].current != 0 && query_mark[OSS[q].a[OSS[q].current].key] == query_mark_stamp)
				OSS[q].current = OSS[q].a[OSS[q].current].next;
			if (OSS[q].current == 0)
				continue;
			int _dist = tmp_dis[q] + OSS[q].a[OSS[q].current].dist;
			if (k < 0 || (dist_k > _dist))
			{
				k = q;
				dist_k = _dist;
			}
		}
		if (k < 0 || dist_k > limit_dis)
			break;
		result.push_back(make_pair(OSS[k].a[OSS[k].current].key, dist_k));

		//--------------delete and update ------------------

		int y = OSS[k].a[OSS[k].current].key;
		query_mark[y] = query_mark_stamp;

	}
	//   cout << "_cnt: " << _cnt << endl;

	return result;
}

extern int farest;
extern int far_v;
bool TkNN::Verify(int x, int target, int top_k){
	query_mark_stamp++;
	int cnt1 = 0, cnt2 = 0, cnt3 = 0;
	bool is_rknn = false;
	int p = belong[x];
	vector<int> anc;
	anc.clear();
	int MAX_DIS = INT_MAX;

	int x2target = distanceQueryFull(x, target);
	if (x2target > farest){
		farest = x2target;
		far_v = x;
	}
	//	return result;
	//
	vector<int> a;
	a.clear();
	a.push_back(p);
	for (int i = 0; i < posSize[p]; i++)
		a.push_back(belong[pos[p][i]]);
	//	printf("a.size(): %d\n", a.size());
	//	printf("MAX_DIS: %d\n", MAX_DIS);
	for (int i = 0; i < a.size(); i++)
	{
		int q = a[i];
		if (OSS[q].size_num <= top_k)
			continue;
		tmp_dis[q] = distanceQuery(uniqueVertex[q], x);
		if (tmp_dis[q] >= MAX_DIS)
			continue;
		if (OSS[q].a[0].next == 0)
			continue;
		int j = OSS[q].a[0].next;
		int _cnt = 1;
		while (_cnt < top_k)
		{
			_cnt++;
			j = OSS[q].a[j].next;
			cnt2++;
		}
		if (MAX_DIS > tmp_dis[q] + OSS[q].a[j].dist)
			MAX_DIS = tmp_dis[q] + OSS[q].a[j].dist;
	}

	int max_height = height[p];
	while (p >= 0 && height[p] >= max_height)
	{
		cnt3++;
		int t = OSS[p].a[0].next;
		tmp_dis[p] = distanceQuery(uniqueVertex[p], x);
		if (t > 0 && MAX_DIS >= tmp_dis[p] + OSS[p].a[t].dist)
		{
			anc.push_back(p);
			OSS[p].current = OSS[p].a[0].next;
			if (OSS[p].size_num >= top_k)
			{
				int v = OSS[p].a[0].previous;
				int tmp = tmp_dis[p] + OSS[p].a[v].dist;
				if (tmp < MAX_DIS)
					MAX_DIS = tmp;
			}
		}
		if (tmp_dis[p] + cloest_higher[p] < MAX_DIS)
			if (max_height > higher[p])
			{
				max_height = higher[p];
			}
		p = pa[p];
	}
	sort(anc.begin(), anc.end(), anc_compare);
	//	printf("anc.size(): %d\n", anc.size());
	p = belong[x];
	int _cnt = 0;
	int find_k_dist;
	for (int i = 0; i < top_k; i++)
	{
		//	printf("i result.size(): %d %d\n", i, result.size());
		//--------------find minimum-------------
		int k = -1, dist_k = -1;
		for (int j = 0; j < anc.size(); j++)
		{
			_cnt++;
			int q = anc[j];
			//	printf("q: %d\n", q);
			while (OSS[q].current != 0 && query_mark[OSS[q].a[OSS[q].current].key] == query_mark_stamp)
				OSS[q].current = OSS[q].a[OSS[q].current].next;
			if (OSS[q].current == 0)
				continue;
			int _dist = tmp_dis[q] + OSS[q].a[OSS[q].current].dist;
			if (k < 0 || (dist_k > _dist))
			{
				k = q;
				dist_k = _dist;
			}
		}
		if (k < 0)
			break;
		find_k_dist = dist_k;
		int y = OSS[k].a[OSS[k].current].key;
		if(y == target){
			is_rknn = true;
			break;
		}
		query_mark[y] = query_mark_stamp;

	}
	if(!is_rknn){
		if(x2target <= find_k_dist)
			is_rknn = true;
	}
	return is_rknn;
}

vector<pair<int, int>> TkNN::query_naive(int x, int top_k)
{
	time_save.clear();
	double pre_time = GetTime();
	query_mark_stamp++;
	int cnt1 = 0, cnt2 = 0, cnt3 = 0;
	vector<pair<int, int>> result;
	result.clear();
	int p = belong[x];

	vector<int> anc;
	anc.clear();
	int MAX_DIS = INT_MAX;
	vector<int> a;
	a.clear();
	a.push_back(p);
	for (int i = 0; i < posSize[p]; i++)
		a.push_back(belong[pos[p][i]]);
	for (int i = 0; i < a.size(); i++)
	{
		int q = a[i];
		if (OSS[q].size_num <= top_k)
			continue;
		tmp_dis[q] = distanceQuery(uniqueVertex[q], x);
		if (tmp_dis[q] >= MAX_DIS)
			continue;
		if (OSS[q].a[0].next == 0)
			continue;
		int j = OSS[q].a[0].next;
		int _cnt = 1;
		while (_cnt < top_k)
		{
			_cnt++;
			j = OSS[q].a[j].next;
			cnt2++;
		}
		if (MAX_DIS > tmp_dis[q] + OSS[q].a[j].dist)
			MAX_DIS = tmp_dis[q] + OSS[q].a[j].dist;
	}

	int max_height = height[p];
	while (p >= 0 && height[p] >= max_height)
	{
		cnt3++;
		int t = OSS[p].a[0].next;
		tmp_dis[p] = distanceQuery(uniqueVertex[p], x);
		if (t > 0 && MAX_DIS >= tmp_dis[p] + OSS[p].a[t].dist)
		{
			anc.push_back(p);
			OSS[p].current = OSS[p].a[0].next;
			if (OSS[p].size_num >= top_k)
			{
				int v = OSS[p].a[0].previous;
				int tmp = tmp_dis[p] + OSS[p].a[v].dist;
				if (tmp < MAX_DIS)
					MAX_DIS = tmp;
			}
		}
		if (tmp_dis[p] + cloest_higher[p] < MAX_DIS)
			if (max_height > higher[p])
			{
				max_height = higher[p];
			}
		p = pa[p];
	}
	sort(anc.begin(), anc.end(), anc_compare);
	p = belong[x];
	cnt_pre_query_time += GetTime() - pre_time;

	//--------------find minimum-------------
	//    int k = -1, dist_k = -1;

	a.clear();
	for (int j = 0; j < anc.size(); j++)
	{
		cnt1++;
		int q = anc[j];
		//    printf("q: %d\n", q);
		int i = 0;
		OSS[q].current = OSS[q].a[0].next;
		while (OSS[q].current != 0 && i < top_k)
		{
			i++;
			//   if (OSS[q].current == 0)
			//     break;
			int t_t = tmp_dis[q] + OSS[q].a[OSS[q].current].dist;
			// distanceQuery(q, OSS[q].a[OSS[q].current].key);
			//  cout << t_t << " " << tmp_dis[q] << " " <<  OSS[q].a[OSS[q].current].dist << endl;
			//     cout << q << " " << OSS[q].a[OSS[q].current].key << endl;
			//  cout << OSS[q].a[OSS[q].current].dist << " " << endl;
			int v = OSS[q].a[OSS[q].current].key;
			//   if (v == -1){
			//        cout << q << " " << OSS[q].current << " " << OSS[q].a[OSS[q].current].key << endl;
			//       while (1);
			//   }
			//    cout << "q: " << q << "  v: " << v << endl;
			if (query_mark[OSS[q].a[OSS[q].current].key] == query_mark_stamp)
			{
				if (t_t < tmp_dis[v])
					tmp_dis[v] = t_t;
			}
			else
			{
				query_mark[v] = query_mark_stamp;
				a.push_back(v);
				tmp_dis[v] = t_t;
			}
			OSS[q].current = OSS[q].a[OSS[q].current].next;
		}
	}
	/*
	for (int i = 0; i < a.size(); i++){
		tmp_dis[a[i]] = distanceQuery(x, a[i]);
		cout << x << " " << a[i] << " " << distanceQuery(a[i], x) << endl;
	}
	  */
	//    cout << "a.size(): " << a.size() << endl;
	//  for (int i = 0; i < a.size(); i++)
	//     cout << a[i] << " " << tmp_dis[a[i]] << endl;
	sort(a.begin(), a.end(), anc_compare);
	for (int i = 0; i < a.size(); i++)
	{
		if (i >= top_k)
			break;
		result.push_back(make_pair(a[i], tmp_dis[a[i]]));
	}

	double current_time = GetTime();
	time_save.push_back(current_time - pre_time);

	return result;
}
vector<pair<int, int>> TkNN::query_delay(int x, int top_k)
{
	//     printf("query_delay\n");
	/*
	for (int i = OSS[1476].a[0].next; i != 0; i = OSS[1476].a[i].next){
		cout << "(" << OSS[1476].a[i].key << ", " << OSS[1476].a[i].dist << "), ";
	}
	 */
	//   cout << "step1"<< endl;
	//   cout << endl;
	time_save.clear();
	double pre_time = GetTime();
	double start_time = pre_time;
	vector<pair<int, int>> result;
	result.clear();
	int p = belong[x];

	//    printf("x, k, p: %d %d %d\n", x, k, p);
	vector<int> anc;
	anc.clear();
	//    printf("OSS[7745].a.size(), OSS[7745].a[0].next: %d %d\n", OSS[7745].a.size(), OSS[7745].a[0].next);
	/*
	int MAX_DIS = INFINITY;
	while (p >= 0){
		int t = OSS[p].a[0].next;
		tmp_dis[p] = distanceQuery(uniqueVertex[p], x);
		//    printf("p :%d\n", p);
	  //  if (t > 0 &&  MAX_DIS > tmp_dis[p] + OSS[p].a[t].dist){
			anc.push_back(p);
			OSS[p].current = OSS[p].a[0].next;
			if (OSS[p].a.size() >= (top_k + 1)){
				int v = OSS[p].a[0].previous;
				int tmp = tmp_dis[p] + OSS[p].a[v].dist;
				if (tmp < MAX_DIS)
				MAX_DIS = tmp;
			}
	 //   }
		p = pa[p];
	}
	 */
	int MAX_DIS = INT_MAX;
	vector<int> a;
	a.clear();
	a.push_back(p);

	//   cout << "step2"<< endl;
	for (int i = 0; i < posSize[p]; i++)
		a.push_back(belong[pos[p][i]]);
	for (int i = 0; i < a.size(); i++)
	{
		int q = a[i];
		if (OSS[q].size_num <= top_k)
			continue;
		tmp_dis[q] = distanceQuery(uniqueVertex[q], x);
		if (tmp_dis[q] >= MAX_DIS)
			continue;
		if (OSS[q].a[0].next == 0)
			continue;
		int j = OSS[q].a[0].next;
		int _cnt = 1;
		while (j != 0 && _cnt < top_k)
		{
			j = OSS[q].a[j].next;
		}
		if (MAX_DIS > tmp_dis[q] + OSS[q].a[j].dist)
			MAX_DIS = tmp_dis[q] + OSS[q].a[j].dist;
	}

	//   cout << "step3"<< endl;
	int max_height = height[p];
	while (p >= 0 && height[p] >= max_height)
	{
		int t = OSS[p].a[0].next;
		tmp_dis[p] = distanceQuery(uniqueVertex[p], x);
		if (t > 0 && MAX_DIS >= tmp_dis[p] + OSS[p].a[t].dist)
		{
			anc.push_back(p);
			OSS[p].current = OSS[p].a[0].next;
			if (OSS[p].size_num >= top_k)
			{
				int v = OSS[p].a[0].previous;
				int tmp = tmp_dis[p] + OSS[p].a[v].dist;
				if (tmp < MAX_DIS)
					MAX_DIS = tmp;
			}
		}
		if (tmp_dis[p] + cloest_higher[p] < MAX_DIS)
			if (max_height > higher[p])
			{
				max_height = higher[p];
			}
		p = pa[p];
	}
	//     printf("anc.size(): %d\n", anc.size());
	//    printf("OSS[7745].current: %d\n", OSS[7745].current);
	sort(anc.begin(), anc.end(), anc_compare);

	//      cout << "step4"<< endl;
	p = belong[x];
	int _cnt = 0;
	for (int i = 0; i < top_k; i++)
	{
		/*
		while (anc.size() > 0){
			int p = anc[anc.size() - 1], t = OSS[p].a[0].next;
			if ((t == 0) || (tmp_dis[p] + OSS[p].a[t].dist > MAX_DIS))
			anc.pop_back();
			else break;
		}
		 */
		//--------------find minimum-------------
		int k = -1, dist_k = -1;
		for (int j = 0; j < anc.size(); j++)
		{
			_cnt++;
			//      cnt1++;
			int q = anc[j];
			//        printf("q: %d\n", q);
			if (OSS[q].current == 0)
				continue;
			int _dist = distanceQuery(uniqueVertex[q], x) + OSS[q].a[OSS[q].current].dist;
			if (k < 0 || (dist_k > _dist))
			{
				k = q;
				dist_k = _dist;
			}
		}
		//    printf("k dist_k: %d %d\n", k, dist_k);
		if (k < 0)
			break;
		int y = OSS[k].a[OSS[k].current].key;
		//   cout << y << " " << dist_k << " " << k << endl;
		result.push_back(make_pair(y, dist_k));

		//--------------delete and update ------------------

		//-----------------------------------

		for (int j = 0; j < anc.size(); j++)
		{
			_cnt++;
			int t = anc[j];

			if (OSS[t].a[OSS[t].current].key == y)
				OSS[t].current = OSS[t].a[OSS[t].current].next;
			int pos = hash.get_value(t, y);
			//         if (t == 1476 && y == 30137)
			//            cout << t << " " << pos << endl;

			if (pos >= 0)
			{
				OSS[t].a[OSS[t].a[pos].previous].next = OSS[t].a[pos].next;
				OSS[t].a[OSS[t].a[pos].next].previous = OSS[t].a[pos].previous;
				OSS[t].trash_can.push_back(pos);
			}
		}

		//    while (1);
		double current_time = GetTime();
		time_save.push_back(current_time - pre_time);
		if ((i + 1) % 5 == 0)
		{
			times_period[((i + 1) / 5) - 1].push_back(current_time - start_time);
			//      printf("%.6lf\n", current_time - start_time);
		}
		//    printf("%.6lf\n", current_time - pre_time);
		pre_time = current_time;
	}
	//    cout << "_cnt: " << _cnt << endl;
	//   get_min_double(time_save);

	for (int i = 0; i < anc.size(); i++)
	{
		int q = anc[i];
		for (int j = OSS[q].trash_can.size() - 1; j >= 0; j--)
		{
			int pos = OSS[q].trash_can[j];
			OSS[q].a[OSS[q].a[pos].next].previous = pos;
			OSS[q].a[OSS[q].a[pos].previous].next = pos;
		}
		OSS[q].trash_can.clear();
	}

	return result;
}
void TkNN::object_setting(int n)
{
	is_current_object = (int *)malloc(sizeof(int) * (n + 1));
	is_current_user = (int *)malloc(sizeof(int) * (n + 1));
	current_distance = (int *)malloc(sizeof(int) * (n + 1));
	group_height = (int *)malloc(sizeof(int) * (n + 1));
	current_state = (int *)malloc(sizeof(int) * (n + 1));
	current_stamp = 0;
	for (int i = 0; i <= n; i++)
		current_state[i] = 0;
}
void TkNN::print(int root)
{
	int cnt = 0;
	for (int i = OSS[root].a[0].next; i != 0; i = OSS[root].a[i].next)
	{
		printf("(%d, %d)", OSS[root].a[i].key, OSS[root].a[i].dist);
		cnt++;
	}
	printf("\n");

	if (cnt != OSS[root].size_num)
	{
		cout << "cnt: " << cnt << "     "
			 << "OSS[root].size_num: " << OSS[root].size_num << endl;
		while (1)
			;
	}
}
bool TkNN::double_objects(int p)
{
	for (int i = OSS[p].a[0].next; i != 0; i = OSS[p].a[i].next)
	{
		if (OSS[p].a[i].key == OSS[p].a[OSS[p].a[i].previous].key)
			return false;
		if ((OSS[p].a[i].next != 0) && (OSS[p].a[OSS[p].a[i].next].dist < OSS[p].a[i].dist))
			return false;
	}
	return true;
}
bool TkNN::check_everyone()
{
	for (int i = 0; i < n; i++)
		if (double_objects(i) == false)
			return false;
	return true;
}

// double get_mean_double(vector<double> &times)
// {
// 	double mean = 0.0;
// 	for (double val : times)
// 	{
// 		mean += val;
// 	}
// 	return mean / times.size();
// }
// double get_var_double(vector<double> &times)
// {
// 	double mean = get_mean_double(times);
// 	double var = 0.0;
// 	for (double val : times)
// 	{
// 		var += (val - mean) * (val - mean);
// 	}
// 	return var / times.size();
// }
// double get_max_double(vector<double> &times)
// {
// 	double max = 0.0;
// 	for (double val : times)
// 	{
// 		if (max < val)
// 			max = val;
// 	}
// 	return max;
// }

/// vector<int> list;
/// int vector<vector<int> > neighbor;
// void recover_number(int p){
//	list.push_back(p);

//}
//-- query NY-t.index NY-t.obj NY-t.query sort
/*
int main(int argc, char *argv[])
{
	srand((int)(time(0)));
	cout << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << endl;
	readIndex(argv[1]);

	//	neighbor.clear();
	//	vector<int> _v;
	//	_v.clear();
	// for (int i = 0; i <= n; i++)
	//		neighbor.push_back(_v);
	//	list.clear();
	//	recover_number(root);

	higher = (int *)malloc(sizeof(int) * (n + 1));
	cloest_higher = (int *)malloc(sizeof(int) * (n + 1));
	for (int i = 0; i < n; i++)
	{
		higher[i] = n;
		cloest_higher[i] = INT_MAX;
		for (int j = 0; j < posSize[i]; j++)
		{
			if (height[belong[pos[i][j]]] < higher[i])
				higher[i] = height[belong[pos[i][j]]];
			int _d = distanceQuery(uniqueVertex[i], pos[i][j]);
			if (cloest_higher[i] > _d)
				cloest_higher[i] = _d;
		}
	}
	//	readGraph(argv[2]);
	//	cout << "Load Graph Finished!" << endl;
	//	printf("height[root]: %d\n", height[root]);
	//	getDFSList();

	//----------------------prepare-------------

	LOG2 = (int *)malloc(sizeof(int) * (n * 2 + 10));
	LOGD = (int *)malloc(sizeof(int) * (n * 2 + 10));
	int k = 0, j = 1;
	for (int i = 0; i < n * 2 + 10; i++)
	{
		if (i > j * 2)
		{
			j *= 2;
			k++;
		}
		LOG2[i] = k;
		LOGD[i] = j;
	}
	insert_type = argv[5];
	query_type = argv[6];

	//---------------------------------
	//-------------test------------------
	kNN knn;
	knn.create_kNN_index();
	//	printf("finish create knn index\n");
	FILE *fobj = fopen(argv[2], "r");
	int number_object;
	fscanf(fobj, "%d", &number_object);
	double start_time = GetTime();
	knn.object_setting(n);

	for (int i = 0; i <= n; i++)
		is_current_object[i] = 0;
	for (int i = 0; i < number_object; i++)
	{
		int x;
		fscanf(fobj, "%d", &x);
		is_current_object[x] = 1;
	}
	int cnt_object = 0;
	for (int i = 1; i <= n; i++)
		if (is_current_object[i] == 1)
			cnt_object++;
	printf("cnt_object: %d\n", cnt_object);
	double bef_init = GetTime();
	knn.initialize_object();
	//  printf("initialize time: %.3lf second\n", GetTime() - bef_init);
	//  return 0;
	printf("knn.object_number[root]: %d\n", knn.object_number[root]);
	double end_time = GetTime();
	printf("object initialization time: %.6lf ms\n", (end_time - start_time) * 1e3);
	//	printf("finish insert object\n");

	//----------------------query------------------------------
	//	printf("hehehe\n");
	if (knn.check_everyone())
	{
		printf("right\n");
	}
	else
		printf("wrong\n");
	vector<double> times;
	times.clear();
	for (int i = 0; i < knn.period; i++)
		knn.times_period[i].clear();
	if (argv[3][1] == 'q')
	{

		//    ofstream out_result;
		//    out_result.open("result.txt");
		FILE *fquery = fopen(argv[4], "r");
		int q_n;
		fscanf(fquery, "%d", &q_n);
		//	printf("n: %d\n", n);
		knn.query_mark = (int *)malloc(sizeof(int) * (n + 1));
		knn.query_mark_stamp = 0;
		for (int i = 0; i <= n; i++)
			knn.query_mark[i] = 0;
		//	printf("start query processing\n");
		start_time = GetTime();
		vector<double> time_array;
		time_array.clear();
		double cnt_delay_time = 0;
		tmp_dis = (int *)malloc(sizeof(int) * (n + 1));
		vector<pair<int, int>> res;
		//   q_n = 1;
		for (int i = 0; i < q_n; i++)
		{

			int x, k;
			fscanf(fquery, "%d %d", &x, &k);
			//	printf("x, k: %d %d\n", x, k);
			double _start_time = GetTime();
			start_time = GetTime();
			if (strcmp(query_type, "delay") == 0)
			{
				res = knn.query_delay(x, k);
			}
			else if (strcmp(query_type, "naive") == 0)
			{
				//    cout << "naive" << endl;
				res = knn.query_naive(x, k);
			}
			else
			{
				res = knn.query(x, k);
			}
			//  res = knn.query_delay(x, k);
			times.push_back(GetTime() - start_time);
			for (int i = 0; i < knn.period; i++)
				if (knn.times_period[i].size() < times.size())
					knn.times_period[i].push_back(GetTime() - start_time);
			//      printf("%.6lf\n", GetTime() - start_time);
			//    printf("------------\n");
			//     if (!knn.check_knn(res, x, k)){
			//         cout << x " " << k << endl;
			//        stop();
			//     }
			//     if (false){
			//         for (int j = 0; j < k; j++)
			//             out_result << res[j].first << " " << res[j].second << " ";
			//        out_result << endl;
			//     }

			double _end_time = GetTime();
			time_array.push_back(_end_time - _start_time);

			bool to_print = false;
			if (to_print)
			{
				for (int j = 0; j < res.size(); j++)
					printf("(%d, %d) ", res[j].first, res[j].second);
				printf("\n");
			}

			double max_delay_time = 0.0;
			for (int j = 0; j < knn.time_save.size(); j++)
			{
				if (knn.time_save[j] > max_delay_time)
					max_delay_time = knn.time_save[j];
				//    printf("%.6lf ", knn.time_save[j]);
			}
			//   printf("\n");
			//   printf("result.size(): %d\n", res.size());
			cnt_delay_time += max_delay_time;
		}

		end_time = GetTime();
		double ave_time = (end_time - start_time) / q_n;
		printf("Average query time: %.6lf ms\n", get_mean_double(times) * 1e3);
		printf("Var query time: %.6lf ms\n", get_var_double(times) * 1e3);
		//    for (int i = 0; i < knn.period; i++)
		//         printf("Average query top %d time: %.6lf ms\n", i, get_mean_double(knn.times_period[i]) * 1e3);
		//	printf("Average query time: %.6lf ms\n", ave_time * 1e3);
		//    printf("Average pre time: %.6lf ms\n", cnt_pre_query_time / double(q_n) * 1e3);
		//   printf("Average delay query time: %.6lf ms\n", (cnt_delay_time / double(q_n)) * 1e3);
		double var_time;
		var_time = 0.0;

		//   out_result.close();
		//  for (int i = 0; i < time_array.size(); i++)
		//       var_time += (time_array[i] - ave_time) * (time_array[i] - ave_time);
		//  var_time /= q_n;
		//   printf("Variance query time: %.6lf ms \n", var_time * 1e3);
		//   knn.print(root);
	}
	else if (argv[3][1] == 'u')
	{

		FILE *fupdate = fopen(argv[4], "r");
		int u_n;
		fscanf(fupdate, "%d", &u_n);
		u_n = 100;
		char st[20];

		for (int i = 0; i < u_n; i++)
		{
			fscanf(fupdate, "%s", st);
			//	cout << st << endl;
			//	fscanf(fupdate, "%s", st);
			//	cout << st << endl;
			//	while (1);
			int x;
			fscanf(fupdate, "%d", &x);
			//	cout << st << " " << x << endl;
			//	while (1);

			start_time = GetTime();
			if (st[0] == 'i')
			{
				//		cout << "insert: " << x << endl;
				knn.insert(x);
			}
			else
			{
				//		cout << "delete: " << x << endl;
				knn.delete_object(x);
			}
			times.push_back(GetTime() - start_time);

			//	if (knn.check_everyone()){
			//		printf("right\n");
			//	}
			//	else {
			//		printf("wrong\n");
			//		while (1);
			//	}
			//	fscanf(fupdate, "%s", st);
		}

		end_time = GetTime();
		double ave_time = (end_time - start_time) / u_n;
		printf("Average update time: %.6lf ms\n", get_mean_double(times) * 1e3);
		printf("Var update time: %.6lf ms\n", get_var_double(times) * 1e3);
		cout << reset_times << endl;

		bool check = false;
		if (check)
		{
			kNN knn2;
			knn2.create_kNN_index();
			//  knn2.object_setting(n);
			knn2.initialize_object();
			for (int i = 0; i < n; i++)
			{
				int v = knn2.OSS[i].a[0].next;
				int cnt = 0;
				for (int j = knn.OSS[i].a[0].next; j != 0; j = knn.OSS[i].a[j].next)
				{
					cnt++;
					printf("i, j, v, cnt: %d %d %d %d\n", i, j, v, cnt);
					printf("(%d, %d) :(%d, %d)\n", knn.OSS[i].a[j].key,
						   knn.OSS[i].a[j].dist, knn2.OSS[i].a[v].key, knn2.OSS[i].a[v].dist);
					if (cnt >= MAX_K)
						break;
					if (knn.OSS[i].a[j].dist != knn2.OSS[i].a[v].dist)
					{
						printf("%d %d\n", knn.OSS[0].size_num, knn2.OSS[0].size_num);

						char c;
						cin >> c;
					}

					v = knn2.OSS[i].a[v].next;
				}
			}
		}
		//	knn.print(root);
		//	if (knn.check_everyone()){
		//	printf("right\n");
		//}
		// else printf("wrong\n");
	}
	else
	{
		FILE *fupdate = fopen(argv[4], "r");
		int u_n;
		fscanf(fupdate, "%d", &u_n);
		char st[20];

		start_time = GetTime();
		for (int i = 0; i < u_n; i++)
		{
			fscanf(fupdate, "%s", st);
			int x;
			fscanf(fupdate, "%d", &x);
			if (st[0] == 'i')
			{
				is_current_object[x] = 1;
			}
			else
			{
				is_current_object[x] = 0;
			}
		}
		for (int i = 1; i <= n; i++)
			if (is_current_object[i] == 1)
				cnt_object++;
		cout << "cnt_object: " << cnt_object << endl;
		//	cout << "is_current_object[190524]: " << is_current_object[190524] << endl;
		knn.initialize_object();
		//	knn.print(root);
		if (knn.check_everyone())
		{
			printf("right\n");
		}
		else
			printf("wrong\n");
	}
}
*/
