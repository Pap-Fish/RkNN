#include <cstdio>
#include <cstring>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <queue>
#include <sys/time.h>
#include <vector>
#include <cmath>
#include <set>
#include <fstream>
#include "../SMPV.h"
#include "../graph.h"
using namespace std;

/*
argv[1]: graph file
argv[2]: queries set or objects set percentage
argv[3]: users set or POIs set percentage
argv[4]: output dir
// argv[5]: smpv file. It means the users set percentage is for subgraph.
*/
int main(int argc, char *argv[])
{
    srand((int)(time(0)));
    FILE *file = fopen(argv[1], "r");
    int n, m;
    fscanf(file, "%d %d", &n, &m);
    fclose(file);

    char *endPtr;
    double querySet_percentage = strtod(argv[2], &endPtr);
    if (endPtr == argv[2] || *endPtr != '\0')
    {
        printf("-----------parameter error-------------");
        exit(1);
        ;
    }
    int querySet_number = (int)(round(n * querySet_percentage));
    printf("querySet number: %d\n", querySet_number);

    int *is_object = (int *)malloc(sizeof(int) * (n + 1));
    for (int i = 0; i <= n; i++)
    {
        is_object[i] = 0;
    }

    double userSet_percentage = strtod(argv[3], &endPtr);
    if (endPtr == argv[3] || *endPtr != '\0')
    {
        printf("-----------parameter error-------------");
        exit(1);
        ;
    }
    int userSet_number = 0;
    // if (argc <= 4)
    // {
        userSet_number = (int)(round(n * userSet_percentage));
        printf("userSet number: %d\n", userSet_number);

        cout << "-------------------------start generate--------------------" << endl;

        string str = argv[1];
        size_t hyphenPos = str.find('-');
        size_t secondhyphenPos = str.find('-', hyphenPos + 1);
        if (secondhyphenPos != string::npos)
        {
            hyphenPos = secondhyphenPos;
        }
        string grapgh_name = str.substr(0, hyphenPos);
        stringstream ss;
        ss << argv[4] << "/" << grapgh_name << "-d." << querySet_percentage << "-" << userSet_percentage << ".obj";
        string out_file = ss.str();
        FILE *fout = fopen(out_file.c_str(), "w");
        fprintf(fout, "%d %d\n", querySet_number, userSet_number);
        for (int i = 0; i < querySet_number; i++)
        {
            int randomOb = (rand() % n) + 1;
            while (is_object[randomOb] != 0)
                randomOb = (rand() % n) + 1;
            is_object[randomOb] = 1;
            fprintf(fout, "%c %d\n", 'q', randomOb);
        }
        for (int i = 0; i < userSet_number; i++)
        {
            int randomUser = (rand() % n) + 1;
            while (is_object[randomUser] != 0)
                randomUser = (rand() % n) + 1;
            is_object[randomUser] = 2;
            fprintf(fout, "%c %d\n", 'u', randomUser);
        }
        fclose(fout);
    // }
    // else
    // {
    //     SMPV smpv;
    //     smpv.G = Graph();
    //     smpv.G.n = n;
    //     smpv.G.m = m;
    //     smpv.readIndex(argv[4]);

    //     vector<int> perUsers(smpv.num_sg, 0);

    //     for(int i=0; i<smpv.num_sg; i++){
    //         perUsers[i] = (int)(round(smpv.SMPVG[i].num * userSet_percentage));
    //         userSet_number += perUsers[i];
    //     }
    //     printf("userSet number: %d, global user rate: %.2lf\n", userSet_number, userSet_number*1.0/n);

    //     cout << "-------------------------start generate--------------------" << endl;
    //     string str = argv[1];
    //     size_t hyphenPos = str.find('-');
    //     size_t secondhyphenPos = str.find('-', hyphenPos + 1);
    //     if (secondhyphenPos != string::npos)
    //     {
    //         hyphenPos = secondhyphenPos;
    //     }
    //     string grapgh_name = str.substr(0, hyphenPos);
    //     stringstream ss;
    //     ss << "./" << grapgh_name << "/" << grapgh_name << "-d." << querySet_percentage << "-smpv" << userSet_percentage << ".obj";
    //     string out_file = ss.str();
    //     FILE *fout = fopen(out_file.c_str(), "w");
    //     fprintf(fout, "%d %d\n", querySet_number, userSet_number);
    //     for (int i = 0; i < querySet_number; i++)
    //     {
    //         int randomOb = (rand() % n) + 1;
    //         while (is_object[randomOb] != 0)
    //             randomOb = (rand() % n) + 1;
    //         is_object[randomOb] = 1;
    //         fprintf(fout, "%c %d\n", 'q', randomOb);
    //     }

    //     for(int i=0; i<smpv.num_sg; i++){
    //         SG sg = smpv.SMPVG[i];
    //         for(int j=0; j<perUsers[i]; j++){
    //             int upos = (rand() % sg.num);
    //             int randomUser = sg.Vertices[upos];
    //             while (is_object[randomUser] != 0){
    //                 upos = (rand() % sg.num);
    //                 randomUser = sg.Vertices[upos];
    //             }
    //             is_object[randomUser] = 2;
    //             fprintf(fout, "%c %d\n", 'u', randomUser);
    //         } 
    //     }
    //     fclose(fout);
    // }

    cout << "-------------------------generate finish--------------------" << endl;
}