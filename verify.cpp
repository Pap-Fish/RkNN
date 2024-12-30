#include <stdio.h>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <map>
#include <algorithm>
#include <ctime>
#include "TkNN.h"
#define MAX_K 100
using namespace std;

static const int SIZEOFINT = 4;
int farest = 0;
int far_v = 0;
extern int *toRMQ, *height, *pa, *uniqueVertex, **RMQIndex;
extern int n;
extern int *belong, *is_current_object, *is_current_user;
extern int root, TreeSize;
extern int **dis, **pos, **pos2;
extern int *posSize, *pos2Size;
extern int *chSize;
extern int **ch;
TkNN knn;
int q_n;
vector<vector<pair<int, int>>> knn_res;
// store the mapping of vertex x to res's position
int *PosTOVertex;
// store the reverse k of each vertex x's rknn
int *res_k;
vector<vector<int>> rknn_res;

char *index_file = "./dataSet/NY-d.index";
char *obj_file = "./dataSet/NY/NY-d.0.1.obj";
char *knn_res_file = "./result/s-r.knn";
char *rknn_res_file = "./result/Tree2.rknn";

string verify_query_file = "";

static bool pair_object_compare(pair<int, int> a, pair<int, int> b)
{
    return a.second < b.second;
}

void print_knn(){
    for(int i = 1; i <= n; i++){
        printf("Vertex %d's %dnn: \n", i, MAX_K);
        for(const auto &pa : knn_res[i]){
            printf("(%d, %d)", pa.first, pa.second);
        }
        printf("\n");
    }
}

void print_rknn(int x){
    int pos;
    for(int i = 0; i < q_n; i++){
        if(PosTOVertex[i] == x){
            pos = i;
            break;
        }
    }
    for(const auto &value : rknn_res[pos]){
        printf("(%d) ", value);
    }
    printf("\n");
}

void graph_knn()
{
    knn.initialize_knn(index_file, obj_file);

    // knn_res.resize(n + 1);


    // for (int i = 1; i <= n; i++)
    // {
        
    //     vector<pair<int, int>> res = knn.query(i, MAX_K);
    //     sort(res.begin(), res.end(), pair_object_compare);
    //     knn_res[i] = res;
    // }

}

// void saveRes(){
//     FILE *fknn = fopen(knn_res_file, "wb");
//     int k = MAX_K;
//     fwrite(&n, SIZEOFINT, 1, fknn);
//     fwrite(&k, SIZEOFINT, 1, fknn);
//     for(int i = 1; i <= n; i++){
//         fwrite(&i, SIZEOFINT, 1, fknn);
//         int len = knn_res[i].size();
//         fwrite(&len, SIZEOFINT, 1, fknn);
//         for(int j = 0; j < len; j++){
//             fwrite(&knn_res[i][j].first, SIZEOFINT, 1, fknn);
//             fwrite(&knn_res[i][j].second, SIZEOFINT, 1, fknn);
//         }

//     }
//     fclose(fknn);
// }

void readKNNRes(){
    FILE *fknn = fopen(knn_res_file, "rb");
    int k;
    fread(&n, SIZEOFINT, 1, fknn);
    printf("n: %d\n", n);
    fread(&k, SIZEOFINT, 1, fknn);
    printf("knn_max_k: %d\n", k);

    // initialize knn_res vector
    knn_res.resize(n + 1);
    for(int i = 0; i <= n; i++){
        knn_res[i].clear();
    }

    for(int i = 1; i <= n; i++){
        int x, len;
        fread(&x, SIZEOFINT, 1, fknn);
        fread(&len, SIZEOFINT, 1, fknn);
        for(int j = 0; j < len; j++){
            int q, dis;
            fread(&q, SIZEOFINT, 1, fknn);
            fread(&dis, SIZEOFINT, 1, fknn);
            pair<int, int> tmp_pair(q, dis);
            knn_res[x].push_back(tmp_pair);
        }
    }
    fclose(fknn);
}

void readReadRKNNRes(){
    FILE *frknn = fopen(rknn_res_file, "rb");
    int r_n;
    fread(&r_n, SIZEOFINT, 1, frknn);
    printf("result num: %d, verify num: %d\n", r_n, q_n);

    rknn_res.clear();
    PosTOVertex = (int *)malloc(sizeof(int) * q_n);
    res_k = (int *)malloc(sizeof(int) * q_n);
    for(int i = 0; i < q_n; i++){
        int x, reverse_k, len;
        fread(&x, SIZEOFINT, 1, frknn);
        fread(&reverse_k, SIZEOFINT, 1, frknn);
        fread(&len, SIZEOFINT, 1, frknn);
        PosTOVertex[i] = x;
        res_k[i] = reverse_k;
        vector<int> ans;
        ans.clear();
        for(int j = 0; j < len; j++){
            int q;
            fread(&q, SIZEOFINT, 1, frknn);
            ans.push_back(q);
        }
        rknn_res.push_back(ans);
    }
    fclose(frknn);
}

/*
argv[1]: index file
argv[2]: obj file
argv[3]: rknn result file
argv[4]: q_n. the number of result want to verify
*/
int main(int argc, char *argv[])
{
    // test if knn_res_file exist
    // ifstream fileStream(knn_res_file);
    // if(!fileStream){
    //     graph_knn();
    //     saveRes();
    // }
    // else{
    //     readKNNRes();
    // } 
    cout << argv[1] << " " << argv[2] << " " << argv[3] << endl;
    index_file = argv[1];
    obj_file = argv[2];
    rknn_res_file = argv[3];

    char *endPtr; // 用于存储转换结束位置的指针
    q_n = strtol(argv[4], &endPtr, 10);
    if (endPtr == argv[4] || *endPtr != '\0')
    {
        printf("-----------parameter error-------------\n");
        exit(1);
    }

    graph_knn();
    // print_knn();
    readReadRKNNRes();
    int cnt_lack = 0;
    vector<pair<int,vector<int>>> ErrorV;
    ErrorV.clear();
    map<int, vector<int>> lack_map;
    lack_map.clear();
    vector<int> x_error;
    // q_n = 100;
    for(int i = 0; i < q_n; i++){
        // printf("i: %d", i);
        int x = PosTOVertex[i];
        int reverse_k = res_k[i];
        int cnt_error = 0;
        x_error.clear();
        if(is_current_object[x] != 1){
            printf("query vertex %d is not object!\n", &x);
            continue;
        }
        if(reverse_k > MAX_K){
            printf("vertex %d's k is too large than %d, can't verify\n", x, reverse_k);
            continue;
        }
        // verify if all the result is correct rknn
        int len = rknn_res[i].size();
        for(int j = 0; j < len; j++){
            int q = rknn_res[i][j];
            if(q == x){
                continue;
            }
            vector<pair<int, int>> q_knn = knn.query(q, MAX_K);
            int pos = reverse_k - 1;
            while(pos < q_knn.size() - 1){
                if(q_knn[pos+1].second!=q_knn[pos].second){
                    break;
                }
                pos++;
            }
            auto end = q_knn.end();
            if(pos != q_knn.size() - 1){
                end = q_knn.begin() + pos + 1;
            }
  
            auto found = find_if(q_knn.begin(), end, [x](const pair<int, int>& e) {
                                    return e.first == x;
                                });
            if(found == end){
                cnt_error++;
                x_error.push_back(q);
                continue;
            }
            // if(i % 100 == 0){
            //     printf("x: %d, q: %d's knn: \n", x, q);
            //     for(int t = 0; t < pos; t++){
            //         printf("vertex: %d, dis: %d\n", knn_res[q][t].first, knn_res[q][t].second);
            //     }
            // }
        }
        if(cnt_error == 0){
            printf("vertex %d's every rknn result is true!\n", x);
        }
        else{
            printf("vertex %d's rknn result has wrong!\n", x);
            ErrorV.push_back(make_pair(x, x_error));
        }

        // verify if all the correct answer is contained
        vector<int> lacks;
        lacks.clear();
        int cnt_size = 0;
        for(int j = 1; j <= n; j++){
            if(is_current_user[j] != 1)
                continue;
            vector<pair<int, int>> j_knn = knn.query(j, MAX_K);
            int pos = reverse_k - 1;
            while(pos < j_knn.size() - 1){
                if(j_knn[pos+1].second!=j_knn[pos].second){
                    break;
                }
                pos++;
            }
            auto end = j_knn.end();
            if(pos != j_knn.size() - 1){
                end = j_knn.begin() + pos + 1;
            }
            auto found = find_if(j_knn.begin(), end, [x](const pair<int, int>& e) {
                                    return e.first == x;
                                });
            if(found == end){
                continue;
            }
            cnt_size++;
            auto it = find(rknn_res[i].begin(), rknn_res[i].end(), j);
            if(it == rknn_res[i].end()){
                lacks.push_back(j);
            }

        }
        if(lacks.size() > 0){
            cnt_lack++;
            lack_map[x] = lacks;
            printf("vertex %d lack %d rknn !\n", x, lacks.size());
        }
        else{
            printf("vertex %d contain all rknn !\n", x);
        }
        printf("vertex %d's rknn result size: %d, correct size: %d\n", x, len, cnt_size);
    }
    

    cout << "---------------------------------------------------------------------" << endl;

    if(ErrorV.size() == 0){
        printf("Correct rate: 100% !\n");
    }
    else{
        double error_rate = ErrorV.size() / static_cast<double>(q_n);
        printf("Error rate: %.4f \n", error_rate);
        printf("The Error query vertex: \n");
        for(int i = 0; i < ErrorV.size(); i++){
            printf("%d: ", ErrorV[i].first);
            // for(const auto &e : ErrorV[i].second){
            //     printf("(%d) ", e);
            // }
            // printf("\n");
        }
        printf("\n");
    }
    cout << "---------------------------------------------------------------------" << endl;
    if(lack_map.size() == 0){
        printf("Lack rate: 0% !\n");
    }
    else{
        double lack_rate = cnt_lack / static_cast<double>(q_n);
        printf("Lack rate: %.4f \n", lack_rate);
        printf("The lack query vertex: \n"); 
        for(const auto &pair : lack_map){
            int q = pair.first;
            printf("vertex %d's rknn: ", q);
            print_rknn(q);
            printf("vertex %d miss: ", q);
            for(int value : pair.second){
                printf("(%d) ", value);
            }
            printf("\n");
        }
    }
}