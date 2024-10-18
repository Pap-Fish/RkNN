#ifndef SMPV_H
#define SMPV_H

#include <cstdio>
#include <vector>
#include <map>
#include <set>
#include "graph.h"
using namespace std;


// struct Graph;
/*
SGid: subgraph id
vertex id: vertex id in original graph
borderid: border vertex id in current subgraph
innerid: inner vertex id in current subgraph
*/
struct SG
{
    int SGid;
    int num;
    vector<int> Vertices;
    vector<int> BorderList;
    vector<int> InnerList;
    // vector<int> users;
    map<int, int> belongBorder;
    map<int, int> belongInner;
    int **BBDT;
    int **IBDT;
    int *Border_mark;
    SG(int sgid);
    void initializeSG();
};


class SMPV
{
public:
    Graph G;
    vector<SG> SMPVG;
    int *belongSG;
    int num_sg;
    SMPV();
    bool initializeSMPV(char *partitionFile);
    void readIndex(char *indexFile);
    void readPartitionFile(char *filename);
    void computeSGTable(int sgid);
    void batchComputeSGTable(int sid, int eid);
    void printSMPV(char *filename);
    void computeSMPVSize();

    Graph inducedGr(int sgid);
    bool isAdjacent(int sg1, int sg2);
    bool shouldPatition(string grfile);
};

#endif // SMPV_H