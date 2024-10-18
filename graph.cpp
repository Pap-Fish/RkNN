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
#include "graph.h"
using namespace std;

double GetTime(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

Graph::Graph()
{
    n = m = 0;
    V.clear();
    E.clear();
}

Graph::Graph(int tmp, char *file)
{
    Graph();
    FILE *fin = fopen(file, "r");
    fscanf(fin, "%d", &n);
    fscanf(fin, "%d", &m);
    for (int i = 0; i <= n; i++)
    {
        map<int, int> v;
        v.clear();
        E.push_back(v);
    }
    for (int i = 0; i < m; i++)
    {
        int x, y, z;
        fscanf(fin, "%d %d %d", &x, &y, &z);
        E[x].insert(make_pair(y, z));
    }
    D.clear();
    D.push_back(0);
    for (int i = 1; i <= n; i++)
        D.push_back(E[i].size());
}
Graph::Graph(char *file)
{
    Graph();
    FILE *fin = fopen(file, "r");
    fscanf(fin, "%d", &n);
    fscanf(fin, "%d", &m);
    for (int i = 0; i <= n; i++)
    {
        map<int, int> v;
        v.clear();
        E.push_back(v);
    }
    for (int i = 0; i < m; i++)
    {
        int x, y, z = 0;
        fscanf(fin, "%d%d%d", &x, &y, &z);
        if (E[x].find(y) != E[x].end())
        {
            if (E[x][y] > z)
            {
                E[x][y] = z;
                E[y][x] = z;
            }
        }
        else
        {
            E[x].insert(make_pair(y, z));
            E[y].insert(make_pair(x, z));
        }
    }
    D.clear();
    D.push_back(0);
    for (int i = 1; i <= n; i++)
        D.push_back(E[i].size());
}

Graph::~Graph(){
}

void Graph::EdgeInitialize()
{
    Edge.clear();
    for (int i = 0; i <= n; i++)
    {
        vector<pair<int, int>> Ed;
        Ed.clear();
        for (map<int, int>::iterator it = E[i].begin(); it != E[i].end(); it++)
        {
            Ed.push_back(*it);
        }
        Edge.push_back(Ed);
    }
}
bool Graph::isEdgeExist(int u, int v)
{
    if (E[u].find(v) == E[u].end())
        return false;
    else
        return true;
}
void Graph::insertEdge(int u, int v, int k)
{
    {
        if (E[u].find(v) != E[u].end())
            return;
        E[u].insert(make_pair(v, k));
        E[v].insert(make_pair(u, k));
        D[u]++;
        D[v]++;
    }
}
void Graph::deleteEdge(int u, int v)
{
    if (E[u].find(v) == E[u].end())
        return;
    E[u].erase(E[u].find(v));
    E[v].erase(E[v].find(u));
    D[u]--;
    D[v]--;
}
void Graph::read_coordinate(char *filename)
{
    printf("read coordinate %s\n", filename);
    X = (int *)malloc(sizeof(int) * (n + 1));
    Y = (int *)malloc(sizeof(int) * (n + 1));
    //	printf("kkk\n");
    FILE *fco = fopen(filename, "r");
    //		printf("kkk\n");
    int tmp;
    fscanf(fco, "%d", &tmp);
    //	printf("tmp n: %d %d\n", tmp, n);
    for (int i = 1; i <= n; i++)
    {
        int v, x, y;
        fscanf(fco, "%d %d %d", &v, &x, &y);
        //		printf("v x y: %d %d %d\n", v, x, y);
        X[v] = x;
        Y[v] = y;
    }
}

PT::PT()
{
}
PT::PT(int _dis, int _x)
{
    dis = _dis;
    x = _x;
}
bool PT::operator<(const PT _pt) const
{
    if (dis == _pt.dis)
        return x > _pt.x;
    return dis > _pt.dis;
}
