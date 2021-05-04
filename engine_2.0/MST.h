#ifndef MST_H
#define MST_H

#include "PrimParams.h"

using namespace std; 

// set to 0 to disable debugging prints
#define DEBUG 0
#define VERBOSE_DEBUG 0


void createWeightGraphs(int *graph_D, int *graph_y, int *graph_x, int *x_pts, int *y_pts, int N);
void writePrimResults(int *parent_nodes, int *child_nodes, int N);
void runPrim(int *graph_D, int *graph_y, int *graph_x, int parent_nodes[], int child_nodes[], int N);

#endif
