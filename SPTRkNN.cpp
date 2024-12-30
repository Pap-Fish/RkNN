#include <cstdio>
#include <cstdlib>
#include <vector>
#include <stack>
#include <algorithm>
// #include <map>
#include <unordered_map>
#include <list>
#include <cmath>
#include <sstream>
#include <future>
#include <chrono>
#include "graph.h"
#include "TkNN.h"
using namespace std;

int *Degree;
int **Neighbor, **Weight;
int N, M;
// Specifies the upper limit of function execution time
const auto timeout = chrono::seconds(35);
const int SIZEOFINT = 4;

int mark_stamp;
int *visit_mark;
TkNN tknn;
extern int *is_current_object, *is_current_user;
int farest;
int far_v;
int cnt_knns, cnt_visit;
int Cnt_User;
bool stopflag;

void readGraph(char *filename)
{
    FILE *file = fopen(filename, "r");
    fscanf(file, "%d %d", &N, &M);
    Degree = (int *)malloc(sizeof(int) * (N + 1));
    vector<vector<pair<int, int>>> nb;
    vector<pair<int, int>> v;
    v.clear();
    for (int i = 0; i <= N; i++)
    {
        //	Degree[i] = 0;
        nb.push_back(v);
    }
    //	cout << n << " " << m << endl;
    for (int i = 0; i < M; i++)
    {
        int x, y, z;
        fscanf(file, "%d %d %d", &x, &y, &z);
        //		Degree[x]++;
        //		cout << x << " " << y << " " << z << endl;
        nb[x].push_back(make_pair(y, z));
    }
    Neighbor = (int **)malloc(sizeof(int *) * (N + 1));
    Weight = (int **)malloc(sizeof(int *) * (N + 1));
    for (int i = 1; i <= N; i++)
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

vector<double> knn_query_times;
bool knn_verify(int x, int target, int top_k)
{
    cnt_knns++;
    double start_time = GetTime();
    bool flag = tknn.Verify(x, target, top_k);
    knn_query_times.push_back(GetTime() - start_time);
    return flag;
}

class NSPT
{ // necessary shortest path tree
public:
    int vnum;
    int root;
    double BuildTime; // record the timecost of BuildTree
    // vector<int> Vertex2Id;
    vector<int> DIST;
    unordered_map<int, int> Vertex2Id;
    vector<vector<int>> CH; // children
    vector<int> Prev;

    NSPT(int _root) : vnum(0), root(_root)
    {
        DIST.clear();
        CH.clear();
        Prev.clear();
    }

    bool ExistNode(int v)
    {
        if (Vertex2Id.find(v) != Vertex2Id.end())
            return true;
        else
            return false;
    }

    void BuildNode(int v, int dist)
    {
        Vertex2Id[v] = vnum++;
        DIST.push_back(dist);
        vector<int> ch;
        ch.clear();
        CH.push_back(ch);
    }

    void BUildTree()
    {
        double start_time = GetTime();
        vector<int> Dist(N + 1, INFINITY);
        Prev.resize(N + 1, -1);
        priority_queue<PT> Q;
        int cnt_user = Cnt_User;

        Dist[root] = 0;
        Q.push(PT(0, root));
        while (!Q.empty())
        {
            PT pt = Q.top();
            Q.pop();
            int v = pt.x;
            int dis = pt.dis;
            if (dis > Dist[v])
            {
                continue;
            }
            for (int i = 0; i < Degree[v]; i++)
            {
                int nb = Neighbor[v][i];
                int ndist = dis + Weight[v][i];
                if (ndist < Dist[nb])
                {
                    Dist[nb] = ndist;
                    Prev[nb] = v;
                    Q.push(PT(ndist, nb));
                }
            }
            if (is_current_user[v] == 1)
            {
                cnt_user--;
                if (cnt_user == 0)
                    break;
            }
        }

        // prune the SPT
        mark_stamp++;
        visit_mark[root] = mark_stamp;
        BuildNode(root, 0);
        // traverse the users vertices
        for (int v = 1; v <= N; v++)
        {
            if (is_current_user[v] != 1 || visit_mark[v] == mark_stamp)
                continue;
            stack<int> path;
            for (int i = v; i != -1; i = Prev[i])
            {
                path.push(i);
                if (visit_mark[i] == mark_stamp)
                    break;
            }

            int pre = root;
            while (!path.empty())
            {
                int cv = path.top();
                path.pop();
                if (visit_mark[cv] == mark_stamp)
                {
                    pre = cv;
                    // if(!ExistNode(cv))
                    //     BuildNode(cv, Dist[cv]);
                    continue;
                }
                BuildNode(cv, Dist[cv]);
                CH[Vertex2Id[pre]].push_back(cv);
                visit_mark[cv] = mark_stamp;
                pre = cv;
            }
        }
        BuildTime = GetTime() - start_time;
    }

    void PrintTree()
    {
        int pathI = 0;
        for (int v = 1; v <= N; v++)
        {
            if (Vertex2Id.find(v) == Vertex2Id.end())
                continue;
            if (CH[Vertex2Id[v]].size() != 0)
                continue;
            vector<int> path;
            for (int i = v; i != -1; i = Prev[i])
            {
                path.push_back(i);
            }
            reverse(path.begin(), path.end());
            pathI++;
            cout << "The " << pathI << " path:" << endl;
            for (int i = 0; i < path.size(); i++)
            {
                int node = path[i];
                cout << node;
                if (i == path.size() - 1)
                    continue;
                int nnode = path[i + 1];
                int dis = DIST[Vertex2Id[nnode]] - DIST[Vertex2Id[node]];
                cout << " --(" << dis << ")--> ";
            }
            cout << endl;
        }
    }

    void PrintSize()
    {
        int vector_data = 0;
        int map_data = 0;

        vector_data += DIST.size();
        vector_data += Prev.size();
        for (int i = 0; i < vnum; i++)
        {
            vector_data += CH[i].size();
        }
        map_data += Vertex2Id.size();
        int total_data = vector_data + map_data;
        float denominator = 1024 * 1024 * 1024;
        float denominator2 = 1024 * 1024;
        printf("NSPT v%d size: %0.3lfGB = %0.3lfMB\n", root, total_data * 4.0 / denominator, total_data * 4.0 / denominator2);
    }

    // query vertex is root
    vector<int> queryRknn(int k)
    {
        cnt_knns = 0;
        farest = 0;
        cnt_visit = 0;
        knn_query_times.clear();

        vector<int> res;
        res.clear();
        stack<int> S;
        S.push(root);

        while (!S.empty() && !stopflag)
        {
            int v = S.top();
            int vid = Vertex2Id[v];
            S.pop();
            cnt_visit++;
            if (DIST[vid] > farest)
                farest = DIST[vid];
            if (knn_verify(v, root, k))
            {
                if (is_current_user[v] == 1)
                    res.push_back(v);
                for (int ch : CH[vid])
                {
                    S.push(ch);
                }
            }
        }
        return res;
    }
};

class NSPT_Pool
{
public:
    int capacity;
    list<shared_ptr<NSPT>> pool;
    unordered_map<int, list<shared_ptr<NSPT>>::iterator> map;

    NSPT_Pool(int _capacity) : capacity(_capacity) {}

    bool isExist(int v)
    {
        if (map.find(v) != map.end())
            return true;
        return false;
    }

    shared_ptr<NSPT> create(int v)
    {
        ;
        if (isExist(v))
            return *map[v];
        if (pool.size() >= capacity)
        {
            auto oldest = pool.back();
            map.erase(oldest->root);
            pool.pop_back();
        }
        auto nspt = make_shared<NSPT>(v);
        pool.push_front(nspt);
        map[v] = pool.begin();
        return nspt;
    }

    shared_ptr<NSPT> get(int v)
    {
        if (!isExist(v))
            return nullptr;
        return *map[v];
    }

    void remove(int v)
    {
        if (isExist(v))
        {
            auto it = map[v];
            pool.erase(it);
            map.erase(v);
        }
    }
};

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

/*
argv[1]: graph file
argv[2]: tree decomposition file
argv[3]: obj file
argv[4]: -q
argv[5]: query file
argv[6]: log file
argv[7]: result file
*/
int main(int argc, char *argv[])
{
    cout << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " ";
    cout << argv[5] << " " << argv[6] << endl;
    readGraph(argv[1]);
    tknn.initialize_knn(argv[2], argv[3]);
    visit_mark = (int *)malloc(sizeof(int) * (N + 1));
    for (int i = 0; i <= N; i++)
        visit_mark[i] = 0;
    mark_stamp = 0;
    Cnt_User = 0;
    for (int i = 1; i <= N; i++)
    {
        if (is_current_user[i] == 1)
            Cnt_User++;
    }

    // NSPT nspt = NSPT(11);
    // nspt.BUildTree();
    // nspt.PrintTree();
    // nspt.PrintSize();
    // NSPT_Pool NSPTs = NSPT_Pool(50);
    // for(int i=1; i<=100; i++){
    //     auto nspt = NSPTs.create(i);
    //     nspt->BUildTree();
    // }
    // auto nspt2 = NSPTs.get(50);
    // if(nspt2 == nullptr)
    //     printf("is null!\n");
    // nspt2->PrintSize();

    if (argv[4][1] == 'q')
    {
        FILE *fquery = fopen(argv[5], "r");
        int q_n;
        fscanf(fquery, "%d", &q_n);
        // NSPT_Pool NSPTs = NSPT_Pool(5000);
        printf("--------------start query processing----------------\n");
        printf("qn: %d\n", q_n);
        // frknn = fopen(argv[7], "wb");
        // fwrite(&q_n, SIZEOFINT, 1, frknn);

        vector<double> qtime;
        qtime.clear();
        vector<double> build_spt_time;
        build_spt_time.clear();
        vector<double> total_time;
        total_time.clear();
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

            NSPT xspt = NSPT(x);
            xspt.BUildTree();
            auto boundQuery = [&xspt](int k)
            { return xspt.queryRknn(k); };
            future<vector<int>> resultFuture = async(launch::async, boundQuery, k);
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

            build_spt_time.push_back(xspt.BuildTime);
            qtime.push_back(_end_time - _start_time);
            total_time.push_back(qtime.back()+build_spt_time.back());
            pfarest.push_back(farest);
            knncalls_cnt.push_back(cnt_knns);
            ave_size.push_back(res.size());
            visited_node.push_back(cnt_visit);

            // printVector(x, k, res);
        }
        printf("\n");
        printf("Average query time: %.6lf ms\n", get_mean_double(qtime) * 1e3);
        printf("Var query time: %.6lf ms\n", get_var_double(qtime) * 1e3);
        printf("Average Necessary Tree build time: %.6lf ms\n", get_mean_double(build_spt_time) * 1e3);
        printf("Average total query time: %.6lf ms\n", get_mean_double(total_time) * 1e3);
        printf("Var total query time: %.6lf ms\n", get_var_double(total_time) * 1e3);
        printf("Average knn calls: %d\n", get_mean_int(knncalls_cnt));
        printf("Average visited vertex or TreeNode: %d\n", get_mean_int(visited_node));
        printf("Average farest visited vertex distance: %d\n", get_mean_int(pfarest));
        printf("out of time: %d\n", cnt_oot);
        printf("Average result size: %d\n", get_mean_int(ave_size));

        // output part
        FILE *fre;
        fre = fopen(argv[6], "a");
        fprintf(fre, "%s %s %s %s\n", argv[1], argv[3], argv[5], "nspt");

        fprintf(fre, "Average query time: %.6lf ms\n", get_mean_double(qtime) * 1e3);
        fprintf(fre, "Var query time: %.6lf ms\n", get_var_double(qtime) * 1e3);
        fprintf(fre, "Average Necessary Tree build time: %.6lf ms\n", get_mean_double(build_spt_time) * 1e3);
        fprintf(fre, "Average total query time: %.6lf ms\n", get_mean_double(total_time) * 1e3);
        fprintf(fre, "Var total query time: %.6lf ms\n", get_var_double(total_time) * 1e3);
        fprintf(fre, "Average knn calls: %d\n", get_mean_int(knncalls_cnt));
        fprintf(fre, "Average visited vertex or TreeNode: %d\n", get_mean_int(visited_node));
        fprintf(fre, "Average farest visited vertex distance: %d\n", get_mean_int(pfarest));
        fprintf(fre, "out of time: %d\n", cnt_oot);
        fprintf(fre, "Average result size: %d\n", get_mean_int(ave_size));

        fclose(fre);
    }
}
