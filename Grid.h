#include <vector>
#include "graph.h"
using namespace std;

class Grid
{
public:
    int minx, maxx, miny, maxy;
    int colSize;	
	int rowSize;	
	double colLength;
	double rowLength;

    Grid();
    Grid(Graph G, int _row, int _col);
    Grid(int x1, int x2, int y1, int y2);
    int width();
    int height();
    void calBorder(Graph G);
    int getCellId(int lon, int lat);
    void buildGrid(Graph G);
    bool isContained(int lon, int lat);

    Graph inducedGr(Graph og, int colnum, int rownum); // build the induced graph from colnum*rownum grid of middle
};