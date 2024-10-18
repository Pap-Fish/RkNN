#include <cstdio>
#include <cstdlib>
#include <vector>
#include <queue>
#include <algorithm>
#include <map>
#include <cmath>
#include <sstream>
#include <future>
#include <chrono>
#include "SMPV.h"
#include "graph.h"
#include "TkNN.h"
using namespace std;
static const int SIZEOFINT = 4;
const char* BUTimeFile = "./SMPVIndex_time.txt";

SMPV smpv;
TkNN tknn;
extern int *is_current_object, *is_current_user;
int *sg_visit_mark, *query_mark;
int visit_mark_stamp;
int farest, far_v;
int cnt_knns, cnt_visit;
int cnt_vsg, cnt_range;

// Specifies the upper limit of function execution time
const auto timeout = chrono::seconds(35);
bool stopflag;

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
int get_max_int(vector<int> &times)
{
    int max = 0.0;
    for (int val : times)
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

void showSMPV()
{
    for (SG sg : smpv.SMPVG)
    {
        printf("----------------------------------------------------------------\n");
        printf("sg id:%d\n", sg.SGid);
        printf("vertices: ");
        for (auto v : sg.Vertices)
            printf("%d ", v);
        printf("\n");
        printf("Border Vertex: ");
        for (auto v : sg.BorderList)
            printf("%d ", v);
        printf("\n");
        printf("Inner Vertex: ");
        for (auto v : sg.InnerList)
            printf("%d ", v);
        printf("\n");
        printf("BBDT: \n");
        for (size_t i = 0; i < sg.belongBorder.size(); i++)
        {
            for (size_t j = 0; j < sg.belongBorder.size(); j++)
                printf("%d ", sg.BBDT[i][j]);
            printf("\n");
        }
        printf("IBDT: \n");
        for (size_t i = 0; i < sg.belongInner.size(); i++)
        {
            for (size_t j = 0; j < sg.belongBorder.size(); j++)
                printf("%d ", sg.IBDT[i][j]);
            printf("\n");
        }
    }
    smpv.computeSMPVSize();
}

vector<double> knn_query_times;

bool knn_verify(int x, int target, int top_k)
{
    cnt_knns++;
    double start_time = GetTime();
    bool flag = tknn.Verify(x, target, top_k);
    knn_query_times.push_back(GetTime()-start_time);
    return flag;
}

vector<int> findUsersInSG(int sg_id)
{
    vector<int> users;
    users.clear();
    SG sg = smpv.SMPVG[sg_id];
    for (size_t i = 0; i < sg.num; i++)
    {
        int u = sg.Vertices[i];
        if (is_current_user[u] == 1)
            users.push_back(u);
    }
    return users;
}

vector<int> query_SMPV(int x, int top_k)
{
    cnt_knns = 0;
    farest = 0;
    cnt_visit = 0;
    cnt_vsg = 0;
    cnt_range=0;
    visit_mark_stamp++;
    knn_query_times.clear();

    priority_queue<PT> PQ;
    set<int> resset;
    resset.clear();
    if (is_current_object[x] != 1)
    {
        vector<int> res;
        return res;
    }

    // StartSG
    int sg_id = smpv.belongSG[x];
    sg_visit_mark[sg_id] = visit_mark_stamp;
    vector<int> users = findUsersInSG(sg_id);
    for (int user : users)
    {
        if (knn_verify(user, x, top_k))
            resset.insert(user);
        query_mark[user] = visit_mark_stamp;
        cnt_visit++;
    }
    for (size_t i = 0; i < smpv.SMPVG[sg_id].BorderList.size(); i++)
    {
        int border = smpv.SMPVG[sg_id].BorderList[i];
        int dist = distanceQueryFull(x, border);
        PQ.push(PT(dist, border));

        smpv.SMPVG[sg_id].Border_mark[i] = visit_mark_stamp;
    }

    bool passport;
    while (!PQ.empty() && !stopflag)
    {
        PT pt = PQ.top();
        PQ.pop();
        int border = pt.x;
        int dis = pt.dis;
        int sg_id = smpv.belongSG[border];
        cnt_visit++;

        vector<pair<int, int>> range = tknn.query(border, top_k, dis);
        cnt_range++;
        cnt_knns++;
        if (dis <= range[range.size() - 1].second)
            passport = true;
        else
            passport = false;
        // for (const auto &element : range)
        // {
        //     int v2 = element.first;
        //     // if(query_mark[v2] == visit_mark_stamp)
        //     //     continue;
        //     if (knn_verify(v2, x, top_k) && is_current_user[v2] == 1)
        //         resset.insert(v2);
        //     query_mark[v2] = visit_mark_stamp;
        // }
        if (passport)
        {
            if (sg_visit_mark[sg_id] != visit_mark_stamp)
            {
                users = findUsersInSG(sg_id);
                for (int user : users)
                {
                    // if(query_mark[user] == visit_mark_stamp)
                    //     continue;
                    if (knn_verify(user, x, top_k))
                        resset.insert(user);
                    query_mark[user] = visit_mark_stamp;
                    cnt_visit++;
                }
                for (size_t i = 0; i < smpv.SMPVG[sg_id].BorderList.size(); i++)
                {
                    int _b = smpv.SMPVG[sg_id].BorderList[i];
                    if (smpv.SMPVG[sg_id].Border_mark[i] == visit_mark_stamp)
                        continue;
                    int dist = distanceQueryFull(x, _b);
                    PQ.push(PT(dist, _b));
                    smpv.SMPVG[sg_id].Border_mark[i] = visit_mark_stamp;
                }
                sg_visit_mark[sg_id] = visit_mark_stamp;
            }
            for (auto it = smpv.G.E[border].begin(); it != smpv.G.E[border].end(); it++)
            {
                int nb = it->first;
                if(smpv.SMPVG[smpv.belongSG[nb]].belongBorder.find(nb) == smpv.SMPVG[smpv.belongSG[nb]].belongBorder.end())
                    continue;
                int bid = smpv.SMPVG[smpv.belongSG[nb]].belongBorder[nb];
                if (sg_visit_mark[smpv.belongSG[nb]] == visit_mark_stamp || smpv.SMPVG[smpv.belongSG[nb]].Border_mark[bid] == visit_mark_stamp)
                    continue;
                int dist = distanceQueryFull(x, nb);
                PQ.push(PT(dist, nb));
                smpv.SMPVG[smpv.belongSG[nb]].Border_mark[bid] = visit_mark_stamp;
            }
        }
    }

    vector<int> res(resset.begin(), resset.end());
    return res;
}

/*
argv[1]: graph file
argv[2]: -b or -r   // build or read
if "-b":
    argv[3]: METIS partition file
else if "-r"
    argv[3]: SMPV Index file

argv[4]: -q
argv[5]: obj file
argv[6]: tree decomposition file
argv[7]: query file
argv[8]: log file
argv[9]: result file
*/
int main(int argc, char *argv[])
{
    smpv.G = Graph(argv[1]);
    printf("--------------read graph finish----------------\n");
    if (argv[2][1] == 'b')
    {
        double st = GetTime();
        bool flag = smpv.initializeSMPV(argv[3]);
        double dt = GetTime();
        if(!flag){
            printf("need further partition!\n");
            return 0;
        }
        printf("build smpv index time: %.6f ms\n", (dt-st)*1e3);
        // FILE *fId = fopen(BUTimeFile, "a");
        // fprintf(fId, "DataSet: %s, SMPV-Index Build time: %.6lf ms\n", argv[1], (dt-st) * 1e3);
        // fclose(fId);
        // return 0;
        
        // save smpv
        string str = argv[1];
        size_t lastSlashPos = str.find_last_of('/');
        if (lastSlashPos != std::string::npos)
        {
            str = str.substr(lastSlashPos + 1);
        }
        size_t hyphenPos = str.find('-');
        size_t secondhyphenPos = str.find('-', hyphenPos + 1);
        if (secondhyphenPos != string::npos)
        {
            hyphenPos = secondhyphenPos;
        }
        string grapgh_name = str.substr(0, hyphenPos);
        stringstream ss;
        ss << "./dataSet/" << grapgh_name << "/" << grapgh_name << "-d.smpv";
        string out_file = ss.str();
        smpv.printSMPV(const_cast<char *>(out_file.c_str()));
    }
    else if (argv[2][1] == 'r')
    {
        double st = GetTime();
        smpv.readIndex(argv[3]);
        double dt = GetTime();
        printf("load smpv index time: %.6f ms\n", (dt-st)*1e3);
    }
    smpv.computeSMPVSize();
    // showSMPV();
    if (argc > 4)
    {
        if (argv[4][1] == 'q')
        {
            // loading tree decomposition
            tknn.initialize_knn(argv[6], argv[5]);
            visit_mark_stamp = 0;
            sg_visit_mark = (int *)malloc(sizeof(int) * (smpv.num_sg));
            query_mark = (int *)malloc(sizeof(int) * (smpv.G.n+1));
            for (size_t i = 0; i < smpv.num_sg; i++)
                sg_visit_mark[i] = 0;
            for(size_t i = 0; i <= smpv.G.n; i++)
                query_mark[i] = 0;
            
            for (SG &sg : smpv.SMPVG)
            {
                sg.Border_mark = (int *)malloc(sizeof(int) * (sg.BorderList.size()));
                for (size_t i = 0; i < sg.BorderList.size(); i++)
                    sg.Border_mark[i] = 0;
            }
            FILE *fquery = fopen(argv[7], "r");
            int q_n;
            fscanf(fquery, "%d", &q_n);
            printf("--------------start query processing----------------\n");
            // frknn = fopen(argv[9], "wb");
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
            vector<int> visited_sg;
            visited_sg.clear();
            vector<double> visited_user;
            visited_user.clear();
            vector<double> border_rate;
            border_rate.clear();
            vector<double> mean_knnTimes;
            mean_knnTimes.clear();
            vector<int> res;
            res.clear();
            int max_vsg=0;
            double max_users=0;
            int barWidth = 100;
            int cnt_oot = 0;
            for (size_t i = 0; i < q_n; i++)
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

                future<vector<int>> resultFuture = async(launch::async, query_SMPV, x, k);
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

                pfarest.push_back(farest);
                time_array.push_back(_end_time - _start_time);
                knncalls_cnt.push_back(cnt_knns);
                ave_size.push_back(res.size());
                visited_node.push_back(cnt_visit);

                int num_user = 0;
                double bors = 0;
                for(int i=0; i<smpv.num_sg; i++){
                    if(sg_visit_mark[i] == visit_mark_stamp){
                        cnt_vsg++;
                        vector<int> u = findUsersInSG(i);
                        num_user += u.size();
                        bors += (smpv.SMPVG[i].BorderList.size()) * 1.0 /smpv.SMPVG[i].num;
                    }
                }
                double per = num_user*1.0/cnt_vsg;
                bors = bors/cnt_vsg;
                if(cnt_vsg>max_vsg)
                    max_vsg = cnt_vsg;
                if(per > max_users)
                    max_users = per;
                visited_sg.push_back(cnt_vsg);
                visited_user.push_back(per);
                border_rate.push_back(bors);
                mean_knnTimes.push_back(get_mean_double(knn_query_times));

                bool to_print = false;
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
            printf("Average visited sg: %d\n", get_mean_int(visited_sg));
            printf("Average users per sg: %.3lf\n", get_mean_double(visited_user));
            printf("Max visited sg: %d\n", max_vsg);
            printf("Max users per sg: %.3lf\n", max_users);
            printf("Ave border rate: %.3lf\n", get_mean_double(border_rate));
            printf("Ave knn time: %.6lf\n", get_mean_double(mean_knnTimes)* 1e3);

            FILE *fre;
            fre = fopen(argv[8], "a");
            fprintf(fre, "%s %s %s %s\n", argv[1], argv[5], argv[7], "smpv");

            fprintf(fre, "Average query time: %.6lf ms\n", get_mean_double(time_array) * 1e3);
            fprintf(fre, "Var query time: %.6lf ms\n", get_var_double(time_array) * 1e3);
            fprintf(fre, "Average knn calls: %d\n", get_mean_int(knncalls_cnt));
            fprintf(fre, "Average visited vertex or TreeNode: %d\n", get_mean_int(visited_node));
            fprintf(fre, "Average farest visited vertex distance: %d\n", get_mean_int(pfarest));
            fprintf(fre, "out of time: %d\n", cnt_oot);
            fprintf(fre, "Average result size: %d\n", get_mean_int(ave_size));
            fprintf(fre, "Average visited sg: %d\n", get_mean_int(visited_sg));
            fprintf(fre, "Average users per sg: %.3lf\n", get_mean_double(visited_user));
            fprintf(fre,"Max visited sg: %d\n", max_vsg);
            fprintf(fre,"Max users per sg: %.3lf\n", max_users);
            fprintf(fre,"Ave border rate: %.3lf\n", get_mean_double(border_rate));
            fprintf(fre,"Ave knn time: %6lf\n", get_mean_double(mean_knnTimes)* 1e3);

            fclose(fre);
        }
    }
}