#include <cstdio>
#include <cstdlib>
#include <vector>
#include <queue>
#include <algorithm>
#include <map>
#include <cmath>
#include <sstream>
#include <functional>
#include <thread>
#include "SMPV.h"
#include "graph.h"
using namespace std;
static const int SIZEOFINT = 4;
const int numThreads = 60;

void printIntVector(vector<int> &a, FILE *fout)
{
    int len = a.size();
    fwrite(&len, SIZEOFINT, 1, fout);
    if (len == 0)
    {
        return;
    }
    for (size_t i = 0; i < len; i++)
    {
        fwrite(&a[i], SIZEOFINT, 1, fout);
    }
}

void printIntMap(map<int, int> &a, FILE *fout)
{
    int len = a.size();
    fwrite(&len, SIZEOFINT, 1, fout);
    if (len == 0)
    {
        return;
    }
    for (size_t i = 0; i < len; i++)
    {
        auto it = next(a.begin(), i);
        fwrite(&it->first, SIZEOFINT, 1, fout);
        fwrite(&it->second, SIZEOFINT, 1, fout);
    }
}

void printIntArray(int **a, int line, int col, FILE *fout)
{
    for(size_t i=0; i<line; i++){
        int *al = a[i];
        fwrite(al, SIZEOFINT, col, fout);
    }
}

void printGr(Graph G, const char *grfile, const char *idMapfile){
    FILE *fout = fopen(grfile, "w");
    fprintf(fout, "%d %d\n", G.n, G.m);
    for(int i=1; i<=G.n; i++){
        for (auto it = G.E[i].begin(); it != G.E[i].end(); it++){
            fprintf(fout, "%d %d %d\n", i, it->first, it->second);
        }
    }
    printf("--------------output graph finish----------------\n");
    FILE *fidmap = fopen(idMapfile, "w");
    fprintf(fidmap, "%d\n", G.n);
    for(int i=1; i<=G.n; i++){
        fprintf(fidmap, "%d %d\n", G.V[i], i);
    }
    printf("--------------output idMap finish----------------\n");
    fclose(fidmap);
    fclose(fout);
}

SG::SG(int sgid)
{
    SGid = sgid;
    num = 0;
    Vertices.clear();
    BorderList.clear();
    InnerList.clear();
    belongBorder.clear();
    belongInner.clear();
    // users.clear();
}

SMPV::SMPV(){}

bool SMPV::initializeSMPV(char *partitionFile)
{
    SMPVG.clear();
    belongSG = (int *)malloc(sizeof(int) * (G.n + 1));
    readPartitionFile(partitionFile);
    if(shouldPatition(partitionFile)){
        return false;
    }
    // determine border vertex
    bool is_border;
    for (int i = 1; i <= G.n; i++)
    {
        SG *sgx = &SMPVG[belongSG[i]];
        is_border = false;
        for (map<int, int>::iterator it = G.E[i].begin(); it != G.E[i].end(); it++)
        {
            int nb = it->first;
            // border vertex connect vertex from different SG
            if (belongSG[i] != belongSG[nb])
            {
                sgx->belongBorder[i] = sgx->BorderList.size();
                sgx->BorderList.push_back(i);
                is_border = true;
                break;
            }
        }
        if (!is_border)
        {
            sgx->belongInner[i] = sgx->InnerList.size();
            sgx->InnerList.push_back(i);
        }
    }
    // build thread pool
    std::vector<std::thread> threads;
    int batch_size = static_cast<int>(ceil(num_sg*1.0/numThreads));
    for (int i = 0; i < numThreads; i++)
    {
        int start = i*batch_size;
        if(start >= num_sg)
            break;
        int end = min((i+1)*batch_size, num_sg)-1;
        threads.push_back(thread([this,start,end](){
            this->batchComputeSGTable(start, end);
        }));
    }
    for (auto& thread : threads) {
        thread.join();
    }

    return true;
}

void SMPV::readIndex(char *indexFile)
{
    FILE *fin = fopen(indexFile, "rb");
    fread(&num_sg, SIZEOFINT, 1, fin);
    SMPVG.clear();
    belongSG = (int *)malloc(sizeof(int) * (G.n + 1));
    for(size_t i=0; i< num_sg; i++){
        SMPVG.push_back(SG(i));
    }
    for(size_t i=0; i< num_sg; i++){
        int sg_id;
        fread(&sg_id, SIZEOFINT, 1, fin);
        SG *sg = &SMPVG[sg_id];
        fread(&sg->num, SIZEOFINT, 1, fin);
        for(size_t j=0; j<sg->num; j++){
            int v;
            fread(&v, SIZEOFINT, 1, fin);
            sg->Vertices.push_back(v);
            belongSG[v] = sg_id;
        }
        // read Border
        int _t;
        fread(&_t, SIZEOFINT, 1, fin);
        for(size_t j=0; j<_t; j++){
            int v;
            fread(&v, SIZEOFINT, 1, fin);
            sg->BorderList.push_back(v);
        }
        fread(&_t, SIZEOFINT, 1, fin);
        for(size_t j=0; j<_t; j++){
            int f,s;
            fread(&f, SIZEOFINT, 1, fin);
            fread(&s, SIZEOFINT, 1, fin);
            sg->belongBorder[f] = s;
        }
        // read Inner
        fread(&_t, SIZEOFINT, 1, fin);
        for(size_t j=0; j<_t; j++){
            int v;
            fread(&v, SIZEOFINT, 1, fin);
            sg->InnerList.push_back(v);
        }
        fread(&_t, SIZEOFINT, 1, fin);
        for(size_t j=0; j<_t; j++){
            int f,s;
            fread(&f, SIZEOFINT, 1, fin);
            fread(&s, SIZEOFINT, 1, fin);
            sg->belongInner[f] = s;
        }
        // read BBDT
        int line, col;
        line = col = sg->belongBorder.size();
        sg->BBDT = (int **)malloc(sizeof(int *) * line);
        for(size_t j=0; j<line; j++){
            sg->BBDT[j] = (int *)malloc(sizeof(int) * col);
            fread(sg->BBDT[j], SIZEOFINT, col, fin);
        }
        line = sg->belongInner.size();
        sg->IBDT = (int **)malloc(sizeof(int *) * line);
        for(size_t j=0; j<line; j++){
            sg->IBDT[j] = (int *)malloc(sizeof(int) * col);
            fread(sg->IBDT[j], SIZEOFINT, col, fin);
        }   
    }
}

void SMPV::readPartitionFile(char *filename)
{
    FILE *fpar = fopen(filename, "r");
    fscanf(fpar, "%d\n", &num_sg);
    for (size_t i = 0; i < num_sg; i++)
    {
        SMPVG.push_back(SG(i));
    }
    for (size_t i = 0; i < G.n; i++)
    {
        int x, sgid;
        fscanf(fpar, "%d %d\n", &x, &sgid);
        SMPVG[sgid].Vertices.push_back(x);
        SMPVG[sgid].num++;
        belongSG[x] = sgid;
    }
}

bool SMPV::shouldPatition(string parfile){
    int cnt_g = 0;
    size_t lastSlashPos = parfile.find_last_of('/');
    if (lastSlashPos != std::string::npos)
    {
        parfile = parfile.substr(lastSlashPos + 1);
    }
    size_t underline = parfile.find('_');
    parfile = parfile.substr(0, underline);
    size_t hyphenPos = parfile.find('-');
    size_t secondhyphenPos = parfile.find('-', hyphenPos + 1);
    if (secondhyphenPos != string::npos){
        hyphenPos = secondhyphenPos;
    }
    string grapgh_name = parfile.substr(0, hyphenPos);
    stringstream ss;
    // check if there has sgs such that num > 2*240
    for(size_t i=0; i< num_sg; i++){
        SG sg = SMPVG[i];
        if(sg.num < 2*240)
            continue;
        cnt_g++;
        ss.str("");
        ss << "./dataSet/" << grapgh_name << "/" << grapgh_name << "-"<< sg.SGid << "-d.gr";
        string outgr = ss.str();
        ss.str("");
        ss << "./dataSet/" << grapgh_name << "/" << grapgh_name << "-"<< sg.SGid << "-idMap.txt";
        string outidMap = ss.str();
        Graph id_G = inducedGr(sg.SGid);
        printGr(id_G, outgr.c_str(), outidMap.c_str());
    }

    if(cnt_g > 0)
        return true;
    else
        return false;
}

void SG::initializeSG()
{
    BBDT = (int **)malloc(sizeof(int *) * (BorderList.size()));
    IBDT = (int **)malloc(sizeof(int *) * (InnerList.size()));
    for (size_t i = 0; i < BorderList.size(); i++)
    {
        BBDT[i] = (int *)malloc(sizeof(int) * (BorderList.size()));
        for (size_t j = 0; j < BorderList.size(); j++)
            BBDT[i][j] = INFINITY;
    }
    for (size_t i = 0; i < InnerList.size(); i++)
    {
        IBDT[i] = (int *)malloc(sizeof(int) * (BorderList.size()));
        for (size_t j = 0; j < BorderList.size(); j++)
            IBDT[i][j] = INFINITY;
    }
}

void SMPV::batchComputeSGTable(int sid, int eid){
    for(int i=sid; i<= eid; i++)
        computeSGTable(i);
}

void SMPV::computeSGTable(int sgid)
{
    SG *sg = &SMPVG[sgid];
    sg->initializeSG();
    int remain;
    map<int,int> belongVr;
    for(size_t i=0; i<sg->num;i++){
        belongVr[sg->Vertices[i]] = i;
    }

    // int *dist = (int *)malloc(sizeof(int) * sg->num);

    // compute BBDT
    for (size_t i = 0; i < sg->BorderList.size(); i++)
    {
        priority_queue<PT> q;
        vector<bool> visited(sg->num, false);
        // for(size_t j=0; j<sg->num;j++)
        //     dist[j] = INFINITY;
        q.push(PT(0, sg->BorderList[i]));
        sg->BBDT[i][i] = 0;
        remain = sg->BorderList.size();
        while (!q.empty())
        {
            PT pt = q.top();
            q.pop();
            int current = pt.x;
            int dist = pt.dis;
            if(visited[belongVr[current]])
                continue;
            visited[belongVr[current]] = true;
            if (sg->belongBorder.find(current) != sg->belongBorder.end())
            {
                sg->BBDT[i][sg->belongBorder[current]] = dist;
                remain--;
                if (remain <= 0)
                    break;
            }
            for (map<int, int>::iterator it = G.E[current].begin(); it != G.E[current].end(); it++)
            {
                int nb = it->first;
                if (belongSG[nb] != sg->SGid)
                    continue;
                if (!visited[belongVr[nb]])
                    q.push(PT(dist + it->second, nb));
            }
        }
    }
    // compute IBDT
    for (size_t i = 0; i < sg->InnerList.size(); i++)
    {
        priority_queue<PT> q;
        vector<bool> visited(sg->num, false);
        q.push(PT(0, sg->InnerList[i]));
        remain = sg->BorderList.size();
        while (!q.empty())
        {
            PT pt = q.top();
            q.pop();
            int current = pt.x;
            int dist = pt.dis;
            if(visited[belongVr[current]])
                continue;
            visited[belongVr[current]] = true;
            if (sg->belongBorder.find(current) != sg->belongBorder.end())
            {
                sg->IBDT[i][sg->belongBorder[current]] = dist;
                remain--;
                if (remain <= 0)
                    break;
            }
            for (map<int, int>::iterator it = G.E[current].begin(); it != G.E[current].end(); it++)
            {
                int nb = it->first;
                if (belongSG[nb] != sg->SGid)
                    continue;
                if (!visited[belongVr[nb]])
                    q.push(PT(dist + it->second, nb));
            }
        }
    }
}

bool SMPV::isAdjacent(int sg1, int sg2){
    bool flag = false;
    SG sgx = SMPVG[sg1];
    for(int i=0; i<sgx.BorderList.size(); i++){
        int bd = sgx.BorderList[i];
        for (map<int, int>::iterator it = G.E[bd].begin(); it != G.E[bd].end(); it++){
            int nb = it->first;
            if(SMPVG[belongSG[nb]].belongBorder.find(nb) == SMPVG[belongSG[nb]].belongBorder.end())
                continue;
            if(belongSG[nb] == sg2){
                flag = true;
                break;
            }
        }
        if(flag)
            break;
    }
    return flag;
}

int findVertexId(Graph G, int v){
    for(int i=0; i<G.V.size(); i++){
        if(G.V[i] == v){
            return i;
        }
    }
    return -1;
}

Graph SMPV::inducedGr(int sgid){
    Graph id_G = Graph();
    SG sg = SMPVG[sgid];
    vector<bool> visited(G.n+1, false);
    int v, m=0;
    int start_v = sg.Vertices[0];

    id_G.V = sg.Vertices;
    id_G.V.insert(id_G.V.begin(), 0); 
    for(size_t i=0; i<id_G.V.size(); i++){
        map<int,int> tm;
        id_G.E.push_back(tm);
    }

    queue<int> Q;
    Q.push(start_v);
    while(!Q.empty()){
        v = Q.front();
        Q.pop();
        if(visited[v])
            continue;
        visited[v] = true;
        int vid = findVertexId(id_G, v);
        for (auto it = G.E[v].begin(); it != G.E[v].end(); it++){
            int nb = it->first;
            if(visited[nb] || belongSG[nb]!=sgid)
                continue;
            int w = it->second;
            int nid = findVertexId(id_G, nb);
            id_G.E[vid].insert(make_pair(nid,w));
            id_G.E[nid].insert(make_pair(vid,w));
            m += 2;
            Q.push(nb);
        }
    }
    id_G.n = id_G.V.size() - 1;
    id_G.m = m;
    return id_G;
}

void SMPV::printSMPV(char *filename)
{
    FILE *fout = fopen(filename, "wb");
    fwrite(&num_sg, SIZEOFINT, 1, fout);
    for (size_t i = 0; i < num_sg; i++)
    {
        SG sg = SMPVG[i];
        fwrite(&sg.SGid, SIZEOFINT, 1, fout);
        printIntVector(sg.Vertices, fout);
        printIntVector(sg.BorderList, fout);
        printIntMap(sg.belongBorder, fout);
        printIntVector(sg.InnerList, fout);
        printIntMap(sg.belongInner, fout);
        int line, col;
        line = col = sg.belongBorder.size();
        printIntArray(sg.BBDT, line, col, fout);
        line = sg.belongInner.size();
        printIntArray(sg.IBDT, line,col,fout);
    }
}

void SMPV::computeSMPVSize(){
    printf("subgraph num: %d\n", num_sg);
    long long cnt_data=0;
    long long cnt_adflist = 0;
    long long total_size = 0;
    int per_num = 0;

    vector<int> sgN;
    int max_sggn = 0;
    int min_sggn = 99999999;
    int null_sg = 0;

    for(size_t i=1; i<=G.n; i++){
        cnt_adflist += G.E[i].size() * 2;
    }
    for(size_t i=0; i< num_sg; i++){
        SG sg = SMPVG[i];
        if(sg.num>(2*240))
            sgN.push_back(sg.SGid);
        if(sg.num > max_sggn)
            max_sggn = sg.num;
        if(sg.num < min_sggn)
            min_sggn = sg.num;
        if(sg.num == 0){
            null_sg++;
            continue;
        }
        per_num += sg.num;
        // cnt_data += sg.num;
        cnt_data += sg.BorderList.size();
        cnt_data += (sg.belongBorder.size()*sg.belongBorder.size()); //BBDT size
        cnt_data += (sg.belongInner.size()*sg.belongBorder.size()); //IBDT size
        total_size += sizeof(int **) * 2;
        total_size += sizeof(int *) * (sg.belongBorder.size() + sg.belongInner.size());
    }

    total_size += (cnt_data + cnt_adflist) * 4.0;
    float denominator = 1024 * 1024 * 1024;
    float denominator2 = 1024 * 1024;
    printf("Valid sg num: %d, average vertices in sg: %d\n", num_sg-null_sg, per_num/(num_sg-null_sg));
    printf("SMPV Index size: %0.3lfGB = %0.3lfMB\n", total_size/denominator, total_size/denominator2);
    printf("Adjacency list size: %0.3lfGB = %0.3lfMB\n", cnt_adflist*4.0/denominator, cnt_adflist*4.0/denominator2);


    printf("max sg num: %d, min sg num: %d\n", max_sggn, min_sggn);
    printf("the num of sg that sgn >= 2*240: %d\n", sgN.size());
    if(sgN.size()>0){
        for(int sid : sgN){
            printf("SG id: %d\n", sid);
            printf("SG num: %d\n", SMPVG[sid].num);
        }
    }
}