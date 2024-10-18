#ifndef GRAPH_H
#define GRAPH_H

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <algorithm>
#include <ctime>
#include <sys/time.h>
#define INFINITY 999999999
using namespace std;

// clock_t ct;

double GetTime(void);
// {
// 	struct timeval tv;
// 	gettimeofday(&tv, NULL);
// 	return tv.tv_sec + tv.tv_usec * 1e-6;
// }

struct Graph
{
	int n, m;
	vector<int> V;
	vector<map<int, int>> E;
	vector<vector<pair<int, int>>> Edge;
	vector<int> D;
	int *X, *Y;
	Graph();
	Graph(int tmp, char *file);
	Graph(char *file);
	~Graph();

	void EdgeInitialize();

	bool isEdgeExist(int u, int v);

	void insertEdge(int u, int v, int k);

	void deleteEdge(int u, int v);

	void read_coordinate(char *filename);
};

struct PT
{
	int dis;
	int x;
	PT();
	PT(int _dis, int _x);
	bool operator<(const PT _pt) const;

};

// struct TN
// {
// 	int height;
// 	int nodeid;
// 	TN(int _height, int _id){
// 		height = _height;
// 		nodeid = _id;
// 	}
// 	bool operator<(const TN _tn) const{
// 		return height < _tn.height;
// 	}
// };

// void readGraph(char *filename);
// void readIndex(char *file);

#endif // GRAPH_H