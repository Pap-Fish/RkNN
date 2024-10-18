#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <queue>
#include <map>
#include "Grid.h"
#include "graph.h"

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

using namespace std;
Graph G;
int *visit_mark;
int visit_mark_stamp;
map<int,int> idMap;

// int findVertexId(Graph G, int v){
//     for(int i=0; i<G.V.size(); i++){
//         if(G.V[i] == v){
//             return i;
//         }
//     }
//     return -1;
// }

Grid::Grid()
{
    minx = 999999999;
    maxx = -999999999;
    miny = 999999999;
    maxy = -999999999;
}

Grid::Grid(Graph G, int _row, int _col) : Grid()
{
    colSize = _col;
    rowSize = _row;
    buildGrid(G);
}

Grid::Grid(int x1, int x2, int y1, int y2) : Grid()
{
    minx = x1;
    maxx = x2;
    miny = y1;
    maxy = y2;
}

int Grid::width()
{
    return maxx - minx;
}

int Grid::height()
{
    return maxy - miny;
}

void Grid::calBorder(Graph G)
{
    for (int i = 1; i <= G.n; i++)
    {
        int lon = G.X[i];
        int lat = G.Y[i];
        minx = MIN(minx, lon);
        maxx = MAX(maxx, lon);
        miny = MIN(miny, lat);
        maxy = MAX(maxy, lat);
    }
}

int Grid::getCellId(int lon, int lat)
{
    int col = int((lon - minx) / colLength);
    int row = int((lat - miny) / rowLength);
    return row * colSize + col;
}

void Grid::buildGrid(Graph G)
{
    calBorder(G);
    colLength = width() * 1.0 / colSize;
    rowLength = height() * 1.0 / rowSize;
}

bool Grid::isContained(int lon, int lat)
{
    if (lon < minx || lon > maxx || lat < miny || lat > maxy)
        return false;
    else
        return true;
}

Graph Grid::inducedGr(Graph og, int colnum, int rownum)
{
    if (colnum == colSize && rownum == rowSize)
        return og;

    visit_mark_stamp++;
    // compute the colnum*rownum middle grid
    int cur_minx, cur_maxx, cur_miny, cur_maxy;
    int leftUp_col = (colSize - colnum) / 2;
    int leftUp_row = (rowSize - rownum) / 2;
    int rightDown_col = leftUp_col + colnum;
    int rightDown_row = leftUp_row + rownum;

    cur_minx = minx + leftUp_col * colLength;
    cur_maxy = maxy - leftUp_row * rowLength;
    cur_maxx = minx + rightDown_col * colLength;
    cur_miny = maxy - rightDown_row * rowLength;

    Grid middle_grid = Grid(cur_minx, cur_maxx, cur_miny, cur_maxy);

    // find first vertex of induced graph
    int start_v = 0;
    for (int i = 1; i <= og.n; i++)
    {
        if (middle_grid.isContained(og.X[i], og.Y[i]))
        {
            start_v = i;
            break;
        }
    }

    for(int i=0; i<=og.n; i++){
        idMap[i] = -1;
    }
    Graph id_G = Graph();
    // id start from 1, pad the index 0
    id_G.V.push_back(0);
    map<int,int> tm0;
    id_G.E.push_back(tm0);
    // the index of v in V is its new id in induced graph
    idMap[start_v] = id_G.V.size();
    id_G.V.push_back(start_v);
    map<int,int> tm;
    tm.clear();
    id_G.E.push_back(tm);
    int v, vid;
    int m=0;
    queue<int> Q;
    Q.push(start_v);

    while (!Q.empty())
    {
        v = Q.front();
        Q.pop();
        if(visit_mark[v] == visit_mark_stamp)
            continue;

        visit_mark[v] = visit_mark_stamp;
        vid = idMap[v];
        for (auto it = og.E[v].begin(); it != og.E[v].end(); it++)
        {
            int nb = it->first;
            if(visit_mark[nb] == visit_mark_stamp)
                continue;
            if (!middle_grid.isContained(og.X[nb], og.Y[nb]))
                continue;
            int w = it->second;
            int nid = idMap[nb];
            if(nid == -1){
                nid = id_G.V.size();
                idMap[nb] = nid;
                id_G.V.push_back(nb);
                map<int,int> tm;
                id_G.E.push_back(tm);
            }
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

void printGr(const char *file, Graph G){
    FILE *fout = fopen(file, "w");
    fprintf(fout, "%d %d\n", G.n, G.m);
    for(int i=1; i<=G.n; i++){
        for (auto it = G.E[i].begin(); it != G.E[i].end(); it++){
            fprintf(fout, "%d %d %d\n", i, it->first, it->second);
        }
    }
    fclose(fout);
}

/**
 * argv[1]: graph file
 * argv[2]: co file
 * argv[3]: col size
 * argv[4]: row size
 *  */

int main(int argc, char *argv[])
{
    cout << "--------------- read graph ---------------" << endl;
    G = Graph(argv[1]);
    cout << "--------------- read co file ---------------" << endl;
    G.read_coordinate(argv[2]);
    cout << "--------------- finish ---------------" << endl;

    char *endPtr;
    int colSize = strtol(argv[3], &endPtr, 10);
    if (endPtr == argv[3] || *endPtr != '\0')
    {
        printf("-----------parameter error-------------\n");
        exit(1);
    }

    int rowSize = strtol(argv[4], &endPtr, 10);
    if (endPtr == argv[4] || *endPtr != '\0')
    {
        printf("-----------parameter error-------------\n");
        exit(1);
    }

    visit_mark = (int *)malloc(sizeof(int)*(G.n+1));
    for(int i=0; i<=G.n; i++)
        visit_mark[i] = 0;
    visit_mark_stamp = 0;

    cout << "--------------- build grid ---------------" << endl;
    Grid grid = Grid(G, rowSize, colSize);
    cout << "--------------- finish ---------------" << endl;

    stringstream ss;
    string dir = "./dataSet/USA";

    for(int i=1; i<10; i++){
        cout << "--------------- start build " << i << " x " << i << " graph ---------------" << endl;
        Graph id_G = grid.inducedGr(G, i, i);
        ss.str("");
        ss << dir << "G" << i << "-d.gr";
        string output = ss.str();
        printGr(output.c_str(), id_G);
    }
}