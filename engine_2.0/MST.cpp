#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <cstdio>
#include <vector>
#include <stack>
#include <math.h>
#include <algorithm> //for min/max

#include "MST.h"
#include "PrimParams.h"

using namespace std; 

// set to 0 to disable debugging prints
#define DEBUG 0
#define VERBOSE_DEBUG 0


// update graph 2D arrays with weights of Manhattan Distance, y difference, x maximum
// row and column dimensions should be the same for all arrays
void createWeightGraphs(int *graph_D, int *graph_y, int *graph_x, int *x_pts, int *y_pts, int N)
{
    int x1, x2, y1, y2; // x and y coordinates of points to compare
    int ptr_val = 0;
    for (int i=0; i < N; ++i) {
        x1 = *(x_pts + i);
        y1 = *(y_pts + i);
        for (int j=0; j < N; ++j) {
            ptr_val = i*N + j;  // used to access next array element
            x2 = *(x_pts + j);
            y2 = *(y_pts + j);
            *(graph_D + ptr_val) = abs(x1 - x2) + abs(y1 - y2); // Manhattan distance
            *(graph_y + ptr_val) = -abs(y1 - y2); // difference in y points
            // x max
            if (x1 > x2) {
                *(graph_x + ptr_val) = -1 * x1;
            } else {
                *(graph_x + ptr_val) = -1 * x2;
            }

            // Print for debug 
            if (VERBOSE_DEBUG) {
                cout << "i: " << i << ", j: " << j << endl;
                cout << "<" << *(graph_D + ptr_val) << ", " << *(graph_y + ptr_val) << ", " << *(graph_x + ptr_val) << ">" << endl;
            }
        }
    }
    return;
}


// write results of Prim's Algorithm to file
// each line of format: parent child 
// -1 is parent of root
void writePrimResults(int *parent_nodes, int *child_nodes, int N)
{
    ofstream results;
    results.open("primResults.txt");
    for (int i = 0; i < N; ++i) { 
        results << *(parent_nodes + i) << " " << *(child_nodes + i) << "\n";
    }
    results.close();
    return;
}


// generate Minimum Spanning Tree
void runPrim(int *graph_D, int *graph_y, int *graph_x, int parent_nodes[], int child_nodes[], int N)
{
    // initialize MST parameters
    PrimParams params(N);
    params.refresh_all();

    int numberConnected = 1;
    parent_nodes[0] = -1;
    child_nodes[0] = 0;
    while (numberConnected != N) {
        // add node to tree
        params.add_node(graph_D, graph_y, graph_x, N);
        params.set_connected(params.get_child());

        parent_nodes[numberConnected] = params.get_parent();
        child_nodes[numberConnected] = params.get_child();
        numberConnected++;

        params.refresh_all(); // reset params for next node
    }

    writePrimResults(parent_nodes, child_nodes, N);

    // Print for debug
    if (DEBUG) {
        for (int i = 0; i < N; ++i) {
            cout << "Parent: " << parent_nodes[i] << ", Child: " << child_nodes[i] << endl;
        }
    }

    return;
}








