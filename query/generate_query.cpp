#include <cstdio>
#include <cstring>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <queue>
#include <sys/time.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <set>
#include <fstream>
using namespace std;

const int SIZEOFINT = 4;
int *generated_object;
// char *obj_file = "../NY-d.gr.0.1.obj";
// char *index_file = "../NY-d.index";
// char *output_dir = "NY-d.gr.rknn-";
int query_times = 10000;
int reverse_k = 10;
int threshold = 1;

int *Degree;
int **Neighbor, **Weight;
int n, m;
void readGraph(char *filename)
{
    FILE *file = fopen(filename, "r");
    fscanf(file, "%d %d", &n, &m);
    printf("n: %d\n", n);
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

/*
argv[1]: dataSet name
argv[2]: obejects percentage
argv[3]: users percentage
argv[4]: reverse_k
argv[5]: query times
argv[6]: random or special
argv[7]: important degree
argv[8]: query percentage of important degree
argv[9]/[7]: output dir
*/
int main(int argc, char *argv[])
{
    srand((int)(time(0)));
    // FILE *fin = fopen(index_file, "rb");
    // fread(&n, SIZEOFINT, 1, fin);
    // printf("n: %d\n", n);

    // reverse_k = atoi(argv[1]);
    char *endPtr; // 用于存储转换结束位置的指针
    reverse_k = strtol(argv[4], &endPtr, 10);
    if (endPtr == argv[4] || *endPtr != '\0')
    {
        printf("-----------parameter error-------------\n");
        exit(1);
    }

    query_times = strtol(argv[5], &endPtr, 10);
    if (endPtr == argv[5] || *endPtr != '\0')
    {
        printf("-----------parameter error-------------\n");
        exit(1);
    }


    stringstream ss;
    ss << "../dataSet/" << argv[1] << "/" << argv[1] << "-d." << argv[2] << "-" << argv[3] << ".obj";
    string obj_file = ss.str();
    ss.str("");
    FILE *fobj = fopen(obj_file.c_str(), "r");
    int object_number, user_number;
    fscanf(fobj, "%d %d\n", &object_number, &user_number);
    printf("object_number: %d\n", object_number);

    // double start_time = GetTime();

    vector<int> objects;

    generated_object = (int *)malloc(sizeof(int) * object_number);
    for (int i = 0; i < object_number; i++)
        generated_object[i] = 0;
    for (int i = 0; i < object_number; i++)
    {
        int x;
        char c;
        fscanf(fobj, "%c %d\n", &c, &x);
        if(c == 'u')
            continue;
        objects.push_back(x);
    }
    cout << "------------------object loading finish----------------\n" << endl;

    if (strcmp(argv[6], "special") == 0)
    {
        ss << "../dataSet/" << argv[1] <<  "-d.gr";
        string g_file = ss.str();
        ss.str("");
        readGraph((char *)g_file.c_str());
        int im_degree = strtol(argv[7], &endPtr, 10);
        if (endPtr == argv[7] || *endPtr != '\0')
        {
            printf("-----------parameter error-------------\n");
            exit(1);
        }
        printf("im_degree: %d\n", im_degree);
        double im_per = strtod(argv[8], &endPtr);
	    if (endPtr == argv[8] || *endPtr != '\0')
	    {
	    	printf("-----------parameter error-------------");
	    	exit(1);
        }
        printf("im_per: %.2f\n", im_per);

        int cnt_od = 0;
        for(int i = 1; i <= n; i++){
            if(Degree[i] >= im_degree && find(objects.begin(), objects.end(), i) != objects.end())
                cnt_od++;
        }
        int im_query_times = (int)(query_times * im_per);
        if(im_query_times > cnt_od){
            threshold += im_query_times / cnt_od;
        }

        ss << argv[9] << "/"  << argv[1] << "-" << argv[2] << "-" << argv[3] << ".rknn-" << im_degree << "-" << im_per << "-" << reverse_k << ".query";
        string output_file = ss.str();
        FILE *fout = fopen(output_file.c_str(), "w");
        fprintf(fout, "%d\n", query_times);

        for(int i = 0; i < query_times-im_query_times; i++){
            int randomNum = rand() % object_number;
            while(Degree[objects[randomNum]] >= im_degree || generated_object[randomNum] >= threshold)
                randomNum = rand() % object_number;
            fprintf(fout, "%d %d\n", objects[randomNum], reverse_k);
            generated_object[randomNum] += 1;
        }
        // generate the querys of special degree node
        for(int i = 0; i < im_query_times; i++){
            int randomNum = rand() % object_number;
            while (Degree[objects[randomNum]] < im_degree || generated_object[randomNum] >= threshold)
                randomNum = rand() % object_number;
            fprintf(fout, "%d %d\n", objects[randomNum], reverse_k);
            generated_object[randomNum] += 1;
        }
        fclose(fout);
    }
    else
    {
        
        ss << argv[7] << "/" << argv[1] << "-" << argv[2] << "-" << argv[3] << ".rknn-" << reverse_k << ".query";
        string output_file = ss.str();
        FILE *fout = fopen(output_file.c_str(), "w");
        fprintf(fout, "%d\n", query_times);

        if (query_times > object_number) 
        {
            threshold += query_times / object_number;
        }

        // generate query

        for (int i = 0; i < query_times; i++)
        {
            int randomNum = rand() % object_number;

            while (generated_object[randomNum] >= threshold)
                randomNum = rand() % object_number;
            int cc = objects[randomNum];
            fprintf(fout, "%d %d\n", objects[randomNum], reverse_k);
            generated_object[randomNum] += 1;
        }
        fclose(fout);
    }

    cout << "------------------query generating finish----------------\n" << endl;

    // FILE *fout = fopen(output_file, "r");

    // int n;
    // fscanf(fout, "%d", &n);
    // cout << "n: " << n << endl;
    // for(int i=0; i<n; i++){
    //     int x, k;
    //     fscanf(fout, "%d %d", &x, &k);
    //     cout << ": " << x << " " << k << endl;
    // }

    fclose(fobj);
}
