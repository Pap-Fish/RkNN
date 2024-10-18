#include <cstdio>
#include <cstring>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <future>
#include <chrono>
#include <functional>
#include "graph.h"
#include "TkNN.h"
using namespace std;
#define MAX_K 100
#define INT_MAX 999999999
#define RESERVE_TIME 1

const int infinity = 999999999;
const int SIZEOFINT = 4;
char *insert_type, *query_type;
const char *BUTimeFile = "./BUIndex_time.txt";

// Specifies the upper limit of function execution time
const auto timeout = chrono::seconds(35);
bool stopflag;

extern int *belong, *height, *pa, *uniqueVertex, *chSize, *posSize;
extern int **ch, **pos;
extern int root, TreeSize;
extern int *is_current_object, *is_current_user;
extern int n;

int cnt_knns, cnt_visit;
int farest;
int far_v;
int *visit_mark, *is_rknn;
int query_mark_stamp;
TkNN tknn;
vector<vector<int>> BUTree;
int deepest = 0;

double get_mean_double(vector<double> &times)
{
    if (times.size() == 0)
    {
        return 0.0;
    }
    double mean = 0.0;
    for (double val : times)
    {
        mean += val;
    }
    return mean / times.size();
}
int get_mean_int(vector<int> &times)
{
    if (times.size() == 0)
    {
        return 0.0;
    }
    double mean = 0.0;
    for (double val : times)
    {
        mean += val;
    }
    return int(round(mean / times.size()));
}
double get_var_double(vector<double> &times)
{
    double mean = get_mean_double(times);
    double var = 0.0;
    for (double val : times)
    {
        var += (val - mean) * (val - mean);
    }
    return var / times.size();
}
double get_max_double(vector<double> &times)
{
    double max = 0.0;
    for (double val : times)
    {
        if (max < val)
            max = val;
    }
    return max;
}

FILE *frknn;
void printVector(int x, int k, vector<int> &a)
{
    fwrite(&x, SIZEOFINT, 1, frknn);
    fwrite(&k, SIZEOFINT, 1, frknn);
    int len = a.size();
    fwrite(&len, SIZEOFINT, 1, frknn);
    if (len == 0)
    {
        return;
    }
    for (int i = 0; i < len; i++)
    {
        fwrite(&a[i], SIZEOFINT, 1, frknn);
    }
}

int *Degree;
int **Neighbor, **Weight;
void readGraph(char *filename)
{
    FILE *file = fopen(filename, "r");
    int n, m;
    fscanf(file, "%d %d", &n, &m);
    Degree = (int *)malloc(sizeof(int) * (n + 1));
    vector<vector<pair<int, int>>> nb;
    vector<pair<int, int>> v;
    v.clear();
    for (int i = 0; i <= n; i++)
    {
        //	Degree[i] = 0;
        nb.push_back(v);
    }
    //	cout << n << " " << m << endl;
    for (int i = 0; i < m; i++)
    {
        int x, y, z;
        fscanf(file, "%d %d %d", &x, &y, &z);
        //		Degree[x]++;
        //		cout << x << " " << y << " " << z << endl;
        nb[x].push_back(make_pair(y, z));
    }
    Neighbor = (int **)malloc(sizeof(int *) * (n + 1));
    Weight = (int **)malloc(sizeof(int *) * (n + 1));
    for (int i = 1; i <= n; i++)
    {
        Degree[i] = nb[i].size();
        Neighbor[i] = (int *)malloc(sizeof(int) * nb[i].size());
        Weight[i] = (int *)malloc(sizeof(int) * nb[i].size());
        for (int j = 0; j < nb[i].size(); j++)
        {
            Neighbor[i][j] = nb[i][j].first;
            Weight[i][j] = nb[i][j].second;
        }
    }
}
void traversal_BUTree(int p)
{
    if (chSize[p] > 0)
        BUTree[p].resize(chSize[p], -1);
    for (int i = 0; i < chSize[p]; i++)
    {
        int child = ch[p][i];
        int chx = uniqueVertex[child];
        traversal_BUTree(child);
        if (chSize[child] <= 1)
        {
            if (is_current_user[chx] == 1)
                BUTree[p][i] = child;
            else if (chSize[ch[p][i]] != 0)
                BUTree[p][i] = BUTree[child][0];
        }
        else if (chSize[child] > 1)
            BUTree[p][i] = child;
    }
}
void construct_BUTree()
{
    for (int i = 0; i < TreeSize; i++)
    {
        vector<int> nextBU;
        nextBU.clear();
        BUTree.push_back(nextBU);
    }
    traversal_BUTree(root);
}

void get_user_subtree(int p, vector<int> &a)
{
    int x = uniqueVertex[p];
    visit_mark[p] = query_mark_stamp;
    if (is_current_user[x] == 1)
        a.push_back(x);
    for (int i = 0; i < chSize[p]; i++)
        if (tknn.user_number[ch[p][i]] > 0)
            get_user_subtree(ch[p][i], a);
}

bool knn_verify(int x, int target, int top_k)
{
    cnt_knns++;
    return tknn.Verify(x, target, top_k);
}

// bool verify_by_block2(int x, int target, int top_k)
// {
//     int spn = SPNQuery(x, target);
//     int x_node = belong[x];
//     int pre_spn = x;

//     while (is_rknn[spn] == -1)
//     {
//         if (spn == pre_spn)
//             break;
//         pre_spn = spn;
//         spn = SPNQuery(spn, target);
//     }
//     if (is_rknn[spn] == 0)
//         return 0;
//     else
//         return knn_verify(x, target, top_k);
// }

// p: tree node, x: query vertex, top_k: k
bool verify_VertexCut(int p, int x, int top_k)
{
    for (int i = 0; i < posSize[p]; i++)
    {
        int v = pos[p][i];
        if (is_rknn[v] == -1)
            is_rknn[v] = knn_verify(v, x, top_k) ? 1 : 0;
        if (is_rknn[v] == 1)
            return true;
    }
    return false;
}

vector<int> query_base(int x, int top_k)
{
    cnt_knns = 0;
    farest = 0;
    cnt_visit = 0;

    priority_queue<PT> H;
    query_mark_stamp++;
    vector<int> result;
    result.clear();
    if (is_current_object[x] != 1)
    {
        return result;
    }

    int *query_mark = (int *)malloc(sizeof(int) * (n + 1));
    for (int i = 0; i <= n; i++)
    {
        query_mark[i] = 0;
    }

    // result.push_back(x);
    visit_mark[x] = query_mark_stamp;
    // query_mark[x] = query_mark_stamp;

    for (int i = 0; i < Degree[x]; i++)
    {
        H.push(PT(Weight[x][i], Neighbor[x][i]));
    }
    bool is_potential;
    while (!H.empty() && !stopflag)
    {
        PT pa = H.top();
        int v = pa.x;
        int dis = pa.dis;
        is_potential = false;
        H.pop();
        if (visit_mark[v] == query_mark_stamp)
        {
            continue;
        }
        if (dis > farest)
            farest = dis;
        vector<pair<int, int>> range = tknn.query(v, top_k, dis);
        cnt_knns++;
        // if(is_current_object[v] == 0)
        // 	range.insert(range.begin(), make_pair(v, 0));
        if (dis <= range[range.size() - 1].second)
        {
            if (is_current_user[v] == 1)
                result.push_back(v);
            is_potential = true;
        }
        query_mark[v] = 1;
        visit_mark[v] = query_mark_stamp;
        for (const auto &element : range)
        {
            int v2 = element.first;
            if (query_mark[v2] != 1)
            {
                if (knn_verify(v2, x, top_k) && is_current_user[v2] == 1)
                {
                    result.push_back(v2);
                }
                query_mark[v2] = 1;
            }
        }
        if (is_potential)
        {
            for (int i = 0; i < Degree[v]; i++)
            {
                if (visit_mark[Neighbor[v][i]] != query_mark_stamp)
                    H.push(PT(dis + Weight[v][i], Neighbor[v][i]));
            }
        }
    }
    // cnt visited vertex
    for (int i = 1; i <= n; i++)
    {
        if (visit_mark[i] == query_mark_stamp)
            cnt_visit++;
    }
    return result;
}

vector<int> query_COREX(int x, int top_k)
{
    cnt_knns = 0;
    farest = 0;
    cnt_visit = 0;

    priority_queue<PT> H;
    query_mark_stamp++;
    vector<int> result;
    result.clear();
    if (is_current_object[x] != 1)
    {
        return result;
    }
    visit_mark[x] = query_mark_stamp;
    for (int i = 0; i < Degree[x]; i++)
    {
        H.push(PT(Weight[x][i], Neighbor[x][i]));
    }
    bool is_open;
    while (!H.empty() && !stopflag)
    {
        PT pa = H.top();
        int v = pa.x;
        int dis = pa.dis;
        H.pop();
        if (visit_mark[v] == query_mark_stamp)
        {
            continue;
        }
        is_open = knn_verify(v, x, top_k);
        visit_mark[v] = query_mark_stamp;
        if (is_open)
        {
            if (is_current_user[v] == 1)
                result.push_back(v);
            for (int i = 0; i < Degree[v]; i++)
            {
                if (visit_mark[Neighbor[v][i]] != query_mark_stamp)
                    H.push(PT(dis + Weight[v][i], Neighbor[v][i]));
            }
        }
    }

    for (int i = 1; i <= n; i++)
    {
        if (visit_mark[i] == query_mark_stamp)
            cnt_visit++;
    }
    return result;
}

void explore_progenyBU(int p, int x, int top_k, vector<int> &res)
{
    if (stopflag || p == -1)
        return;

    int v = uniqueVertex[p];
    visit_mark[p] = query_mark_stamp;
    bool passport = false;
    cnt_visit++;

    for (int i = 0; i < posSize[p]; i++)
    {
        if (is_rknn[pos[p][i]] == 1)
            passport = true;
    }
    if (is_current_user[v] == 1)
    {
        if (is_rknn[v] == -1)
            is_rknn[v] = knn_verify(v, x, top_k) ? 1 : 0;
        if (is_rknn[v] == 1)
        {
            passport = true;
            res.push_back(v);
        }
    }
    if (!passport)
    {
        if (is_rknn[v] == 1)
        {
            passport = true;
        }
        else
        {
            passport = verify_VertexCut(p, x, top_k);
        }
        if (!passport)
        {
            is_rknn[v] = knn_verify(v, x, top_k) ? 1 : 0;
            if (is_rknn[v] == 1)
                passport = true;
        }
    }
    if (passport)
    {
        for (int i = 0; i < chSize[p]; i++)
        {
            int next_branch = BUTree[p][i];
            // if (next_branch != -1 && visit_mark[next_branch] != query_mark_stamp)
            explore_progenyBU(next_branch, x, top_k, res);
        }
    }
}

vector<int> query_onTreeBranch(int x, int top_k)
{
    cnt_knns = 0;
    cnt_visit = 0;
    farest = 0;
    deepest = 0;
    vector<int> res;
    res.clear();
    query_mark_stamp++;

    for (int i = 0; i <= n; i++)
    {
        is_rknn[i] = -1;
    }
    int query_p = belong[x];
    visit_mark[query_p] = query_mark_stamp;
    for (int i = 0; i < chSize[query_p]; i++)
        explore_progenyBU(BUTree[query_p][i], x, top_k, res);

    int q = pa[query_p];
    bool passport;

    while (q != -1)
    {
        if (stopflag)
            break;
        if (visit_mark[q] == query_mark_stamp)
        {
            q = pa[q];
            continue;
        }
        cnt_visit++;
        int v = uniqueVertex[q];
        visit_mark[q] = query_mark_stamp;
        passport = false;
        if (is_current_user[v] == 1)
        {
            if (is_rknn[v] == -1)
                is_rknn[v] = knn_verify(v, x, top_k) ? 1 : 0;
            if (is_rknn[v] == 1)
            {
                res.push_back(v);
            }
        }
        if (chSize[q] > 1)
        {
            for (int i = 0; i < posSize[q]; i++)
            {
                if (is_rknn[pos[q][i]] == 1)
                {
                    passport = true;
                    break;
                }
            }
            if (!passport)
            {
                if (is_rknn[v] == -1)
                    is_rknn[v] = knn_verify(v, x, top_k) ? 1 : 0;
                if (is_rknn[v] == 1)
                {
                    passport = true;
                }
                else
                    passport = verify_VertexCut(q, x, top_k);
            }
            if (passport)
            {
                for (int j = 0; j < chSize[q]; j++)
                {
                    if (visit_mark[ch[q][j]] == query_mark_stamp)
                        continue;
                    explore_progenyBU(BUTree[q][j], x, top_k, res);
                }
            }
            else
                break;
        }
        q = pa[q];
    }

    return res;
}

void explore_progenyBUOpt(int p, int x, int top_k, vector<int> &res)
{
    if (stopflag || p == -1 || tknn.user_number[p] == 0)
        return;
    if (height[p] > deepest)
        deepest = height[p];

    visit_mark[p] = query_mark_stamp;
    int cmp = 0;
    int v = uniqueVertex[p];
    bool passport = false;
    cnt_visit++;
    for (int i = 0; i < posSize[p]; i++)
    {
        if (is_rknn[pos[p][i]] == -1)
            cmp++;
        if (is_rknn[pos[p][i]] == 1)
            passport = true;
    }
    if (tknn.user_number[p] <= cmp)
    {
        vector<int> sub_users;
        sub_users.clear();
        get_user_subtree(p, sub_users);

        for (int i = 0; i < sub_users.size(); i++)
        {
            if (knn_verify(sub_users[i], x, top_k))
                res.push_back(sub_users[i]);
        }
        cnt_visit += sub_users.size();
        return;
    }
    if (is_current_user[v] == 1)
    {
        if (is_rknn[v] == -1)
            is_rknn[v] = knn_verify(v, x, top_k) ? 1 : 0;
        if (is_rknn[v] == 1)
        {
            passport = true;
            res.push_back(v);
        }
    }
    // int cnt_haveUser = 0;
    // for (int i = 0; i < chSize[p]; i++)
    // {
    //     if (tknn.user_number[ch[p][i]] != 0)
    //         cnt_haveUser++;
    // }
    // if (cnt_haveUser <= 1)
    // {
    //     passport = true;
    // }
    if (!passport)
    {
        if (is_rknn[v] == 1)
        {
            passport = true;
        }
        else
        {
            passport = verify_VertexCut(p, x, top_k);
        }
        if (!passport)
        {
            is_rknn[v] = knn_verify(v, x, top_k) ? 1 : 0;
            if (is_rknn[v] == 1)
                passport = true;
        }
    }
    if (passport)
    {
        for (int i = 0; i < chSize[p]; i++)
        {
            int next_branch = BUTree[p][i];
            // if (next_branch != -1 && visit_mark[next_branch] != query_mark_stamp)
            explore_progenyBUOpt(next_branch, x, top_k, res);
        }
    }
}

vector<int> query_onTreeBranchOpt(int x, int top_k)
{
    cnt_knns = 0;
    cnt_visit = 0;
    farest = 0;
    deepest = 0;
    vector<int> res;
    res.clear();
    query_mark_stamp++;

    for (int i = 0; i <= n; i++)
    {
        is_rknn[i] = -1;
    }

    int query_p = belong[x];
    visit_mark[query_p] = query_mark_stamp;
    is_rknn[x] = 1;
    for (int i = 0; i < chSize[query_p]; i++)
        explore_progenyBUOpt(BUTree[query_p][i], x, top_k, res);

    int q = pa[query_p];
    bool passport;
    while (q != -1)
    {
        if (stopflag)
            break;
        if (visit_mark[q] == query_mark_stamp)
        {
            q = pa[q];
            continue;
        }
        cnt_visit++;
        int v = uniqueVertex[q];
        visit_mark[q] = query_mark_stamp;
        passport = false;
        if (is_current_user[v] == 1)
        {
            if (is_rknn[v] == -1)
                is_rknn[v] = knn_verify(v, x, top_k) ? 1 : 0;
            if (is_rknn[v] == 1)
            {
                res.push_back(v);
            }
        }
        if (chSize[q] > 1)
        {
            // int cnt_haveUser=0;
            // for (int j = 0; j < chSize[q]; j++)
            // {
            //     if (visit_mark[ch[q][j]] == query_mark_stamp)
            //         continue;
            //     if (tknn.user_number[ch[q][j]] != 0)
            //         cnt_haveUser++;
            // }
            // if(cnt_haveUser==0)
            //     passport = true;
            for (int i = 0; i < posSize[q]; i++)
            {
                if (is_rknn[pos[q][i]] == 1)
                {
                    passport = true;
                    break;
                }
            }
            if (!passport)
            {
                if (is_rknn[v] == -1)
                    is_rknn[v] = knn_verify(v, x, top_k) ? 1 : 0;
                if (is_rknn[v] == 1)
                {
                    passport = true;
                }
                else
                    passport = verify_VertexCut(q, x, top_k);
            }
            if (passport)
            {
                for (int j = 0; j < chSize[q]; j++)
                {
                    if (visit_mark[ch[q][j]] == query_mark_stamp)
                        continue;
                    explore_progenyBUOpt(BUTree[q][j], x, top_k, res);
                }
            }
            else
                break;
        }
        q = pa[q];
    }

    return res;
}

void explore_TD(int p, int x, int top_k, vector<int> &res)
{
    if (stopflag || p == -1)
        return;
    int v = uniqueVertex[p];
    visit_mark[p] = query_mark_stamp;
    bool passport = false;
    cnt_visit++;

    if (is_rknn[v] == -1)
        is_rknn[v] = knn_verify(v, x, top_k) ? 1 : 0;
    if (is_rknn[v] == 1)
    {
        passport = true;
        if (is_current_user[v] == 1)
            res.push_back(v);
    }
    else
    {
        passport = verify_VertexCut(p, x, top_k);
    }
    if (passport)
    {
        for (int i = 0; i < chSize[p]; i++)
        {
            int next_child = ch[p][i];
            // if (next_branch != -1 && visit_mark[next_branch] != query_mark_stamp)
            explore_TD(next_child, x, top_k, res);
        }
    }
}

vector<int> query_onTreeTD(int x, int top_k)
{
    cnt_knns = 0;
    cnt_visit = 0;
    farest = 0;
    deepest = 0;
    vector<int> res;
    res.clear();
    query_mark_stamp++;

    for (int i = 0; i <= n; i++)
    {
        is_rknn[i] = -1;
    }
    int query_p = belong[x];
    is_rknn[x] = 1;
    visit_mark[query_p] = query_mark_stamp;
    for (int i = 0; i < chSize[query_p]; i++)
        explore_TD(BUTree[query_p][i], x, top_k, res);

    int q = pa[query_p];
    bool passport;

    while (q != -1)
    {
        if (stopflag)
            break;
        if (visit_mark[q] == query_mark_stamp)
        {
            q = pa[q];
            continue;
        }
        cnt_visit++;
        int v = uniqueVertex[q];
        visit_mark[q] = query_mark_stamp;
        passport = false;
        if (is_current_user[v] == 1)
        {
            if (is_rknn[v] == -1)
                is_rknn[v] = knn_verify(v, x, top_k) ? 1 : 0;
            if (is_rknn[v] == 1)
            {
                res.push_back(v);
            }
        }
        if (chSize[q] > 1)
        {
            if (is_rknn[v] == -1)
            {
                is_rknn[v] = knn_verify(v, x, top_k) ? 1 : 0;
            }
            if (is_rknn[v] == 1)
            {
                passport = true;
            }
            else
            {
                passport = verify_VertexCut(q, x, top_k);
            }
            if (passport)
            {
                for (int j = 0; j < chSize[q]; j++)
                {
                    if (visit_mark[ch[q][j]] == query_mark_stamp)
                        continue;
                    explore_TD(ch[q][j], x, top_k, res);
                }
            }
            else
                break;
        }
        q = pa[q];
    }

    return res;
}

void computeIndexSize()
{
    long long BU_data = 0;
    long long td_data = 0;
    long long H2H_data = 0;
    int max_w = 0;

    for (int i = 0; i < TreeSize; i++)
    {
        td_data += posSize[i];
        td_data += chSize[i];
        BU_data += BUTree[i].size();
        H2H_data += (height[i] - 1);
        if (posSize[i] - 1 > max_w)
            max_w = posSize[i] - 1;
    }
    td_data += TreeSize; // the size of parent array
    long long total_data = td_data + BU_data;
    long long total_H2H = total_data + H2H_data;
    float denominator = 1024 * 1024 * 1024;
    float denominator2 = 1024 * 1024;
    printf("BU-Index size: %0.3lfGB = %0.3lfMB\n", total_data * 4.0 / denominator, total_data * 4.0 / denominator2);
    printf("BU-Index with H2H size: %0.3lfGB = %0.3lfMB\n", total_H2H * 4.0 / denominator, total_H2H * 4.0 / denominator2);
    printf("max tree width: %d\n", max_w);
}

void printBUTree()
{
    char *out_put = "./BUTree.index";
    FILE *fTree = fopen(out_put, "w");
    for (int i = 0; i < TreeSize; i++)
    {
        int BUsize = BUTree[i].size();
        fprintf(fTree, "%d %d\n", i, BUsize);
        for (int j = 0; j < BUsize; j++)
            fprintf(fTree, "%d ", BUTree[i][j]);
        fprintf(fTree, "\n");
    }
    fclose(fTree);
}

void printTree()
{
    for (int i = 0; i < TreeSize; i++)
    {
        printf("----------------------------------------\n");
        printf("Tree node: %d\n", i);
        printf("uniqueVertex: %d, ", uniqueVertex[i]);
        printf("vert: ");
        for (int j = 0; j < posSize[i]; j++)
        {
            printf("%d, ", pos[i][j]);
        }
        printf("\n");
        printf("child: ");
        for (int j = 0; j < chSize[i]; j++)
        {
            printf("%d, ", uniqueVertex[ch[i][j]]);
        }
        printf("\n");
    }
}

void test()
{
    int x = 80125;
    int b_dist = 389584, e_dist = 272823;
    int cnt_bs = 0, cnt_es = 0;
    for (int i = 1; i <= n; i++)
    {
        int p = belong[i];
        int _dist = distanceQueryFull(x, i);
        if (_dist < b_dist)
        {
            if (is_current_user[i])
                cnt_bs++;
        }
        if (_dist < e_dist)
        {
            if (is_current_user[i])
                cnt_es++;
        }
    }
    printf("branch scope: %d\n", cnt_bs);
    printf("eager scope: %d\n", cnt_es);
}

int main(int argc, char *argv[])
{
    srand((int)(time(0)));
    cout << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " ";
    cout << " " << argv[5] << " " << argv[6] << " " << argv[7] << " normal" << endl;
    // cout << argv[8] << endl;
    readGraph(argv[1]);
    tknn.initialize_knn(argv[2], argv[3]);

    // int spn = SPNQuery(7, 12);

    // test construct time
    // double start_time = GetTime();
    construct_BUTree();
    // double end_time = GetTime();
    // printf("BU index initialization time: %.6lf ms\n", (end_time - start_time) * 1e3);
    // FILE *fId = fopen(BUTimeFile, "a");
    // fprintf(fId, "DataSet: %s, BU-Index Build time: %.6lf ms\n", argv[1], (end_time - start_time) * 1e3);
    // fclose(fId);
    // return 0;

    // printTree();
    // printBUTree();
    // test();
    // computeIndexSize();

    is_rknn = (int *)malloc(sizeof(int) * (n + 1));

    query_type = argv[6];

    if (argv[4][1] == 'q')
    {

        //    ofstream out_result;
        //    out_result.open("result.txt");
        FILE *fquery = fopen(argv[5], "r");
        int q_n;
        fscanf(fquery, "%d", &q_n);
        //	printf("n: %d\n", n);
        // cknn.query_mark = (int *)malloc(sizeof(int) * (n + 1));
        // for (int i = 0; i <= n; i++)
        //     cknn.query_mark[i] = 0;
        query_mark_stamp = 0;
        visit_mark = (int *)malloc(sizeof(int) * (n + 1));
        for (int i = 0; i <= n; i++)
            visit_mark[i] = 0;
        printf("--------------start query processing----------------\n");
        // frknn = fopen(argv[8], "wb");
        // fwrite(&q_n, SIZEOFINT, 1, frknn);
        printf("qn: %d\n", q_n);

        vector<double> time_array;
        time_array.clear();
        vector<int> knncalls_cnt;
        knncalls_cnt.clear();
        vector<int> ave_size;
        ave_size.clear();
        vector<int> pfarest;
        pfarest.clear();
        vector<int> visited_node;
        visited_node.clear();
        vector<int> res;
        res.clear();

        double start_time = GetTime();
        int cnt_oot = 0;
        int barWidth = 100;
        int cnt_ttt = 0;
        auto boundQuery = [&](int x, int k)
        {
            if (strcmp(query_type, "branch") == 0)
                return query_onTreeBranch(x, k);
            else if (strcmp(query_type, "branchopt") == 0)
                return query_onTreeBranchOpt(x, k);
            else if (strcmp(query_type, "td") == 0)
                return query_onTreeTD(x, k);
            else if (strcmp(query_type, "corex") == 0)
                return query_COREX(x, k);
            else
                return query_base(x, k);
        };

        for (int i = 0; i < q_n; i++)
        {
            // Progress bar
            float progress = (float)(i + 1) / q_n;
            int pos = barWidth * progress;
            cout << " [";
            for (int b = 0; b < barWidth; b++)
            {
                if (b < pos)
                    cout << "=";
                else if (b == pos)
                    cout << ">";
                else
                    cout << " ";
            }
            cout << "] " << int(progress * 100.0) << "% / " << i + 1 << "\r";
            cout.flush();
            stopflag = false;

            int x, k;
            fscanf(fquery, "%d %d", &x, &k);

            // vector<pair<int, int>> rr = tknn.query(80709, 10);

            future<vector<int>> resultFuture = async(launch::async, boundQuery, x, k);
            // Gets the current time
            auto cstart_time = chrono::steady_clock::now();
            double _start_time = GetTime();
            while (true)
            {
                if (resultFuture.wait_for(chrono::seconds(0)) == future_status::ready)
                {
                    res = resultFuture.get();
                    break;
                }
                // compute function execution time
                auto cend_time = chrono::steady_clock::now();
                auto elapsed_time = chrono::duration_cast<std::chrono::seconds>(cend_time - cstart_time);
                // check if exceed the specified time
                if (elapsed_time >= timeout)
                {
                    stopflag = true;
                    // Wait for the asynchronous task to complete
                    resultFuture.wait();
                    break;
                }
            }
            if (stopflag)
            {
                cnt_oot++;
                continue;
            }

            double _end_time = GetTime();

            // int far_n = belong[far_v];
            // int lcca = LCAQuery(far_n, belong[x]);
            // if(lcca == far_n)
            //     cnt_ttt++;

            pfarest.push_back(farest);

            time_array.push_back(_end_time - _start_time);
            knncalls_cnt.push_back(cnt_knns);
            ave_size.push_back(res.size());
            visited_node.push_back(cnt_visit);

            bool to_print = false;
            // if(i >= q_n - 5)
            // 	to_print = true;
            if (to_print)
            {
                printf("Vertex %d's r%dnn: ", x, k);
                for (int j = 0; j < res.size(); j++)
                    printf("(%d) ", res[j]);
                printf("\n");
            }
             // printVector(x, k, res);
        }
        printf("\n");

        printf("Average query time: %.6lf ms\n", get_mean_double(time_array) * 1e3);
        printf("Var query time: %.6lf ms\n", get_var_double(time_array) * 1e3);
        printf("Average knn calls: %d\n", get_mean_int(knncalls_cnt));
        printf("Average visited vertex or TreeNode: %d\n", get_mean_int(visited_node));
        printf("Average farest visited vertex distance: %d\n", get_mean_int(pfarest));
        printf("out of time: %d\n", cnt_oot);
        printf("Average result size: %d\n", get_mean_int(ave_size));

        // output part
        FILE *fre;
        fre = fopen(argv[7], "a");
        fprintf(fre, "%s %s %s %s %s\n", argv[1], argv[3], argv[5], argv[6], "normal");

        fprintf(fre, "Average query time: %.6lf ms\n", get_mean_double(time_array) * 1e3);
        fprintf(fre, "Var query time: %.6lf ms\n", get_var_double(time_array) * 1e3);
        fprintf(fre, "Average knn calls: %d\n", get_mean_int(knncalls_cnt));
        fprintf(fre, "Average visited vertex or TreeNode: %d\n", get_mean_int(visited_node));
        fprintf(fre, "Average farest visited vertex distance: %d\n", get_mean_int(pfarest));
        fprintf(fre, "out of time: %d\n", cnt_oot);
        fprintf(fre, "Average result size: %d\n", get_mean_int(ave_size));
        /*
                fprintf(fre, "cct: %d\n", cct);
                fprintf(fre, "base visited: %.2lf, tree visited: %.2lf, tree2 visited: %.2lf\n", vbase*1.0/q_n, vtree*1.0/q_n, vtree2*1.0/q_n);
                fprintf(fre, "base time: %.6lf, tree time: %.6lf, tree2 time: %.6lf\n", get_mean_double(t1)* 1e3, get_mean_double(t2)* 1e3, get_mean_double(t3)* 1e3);
        */
        fclose(fre);

        double var_time;
        var_time = 0.0;

        //   out_result.close();
        //  for (int i = 0; i < time_array.size(); i++)
        //       var_time += (time_array[i] - ave_time) * (time_array[i] - ave_time);
        //  var_time /= q_n;
        //   printf("Variance query time: %.6lf ms \n", var_time * 1e3);
        //   knn.print(root);
    }
}