#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstdio>
#include <vector>
#include <stack>
#include <math.h>
#include <algorithm> //for min/max

#include "PrimParams.h"
using namespace std;

// constructor for initialization
PrimParams::PrimParams(int N)
{
    for(int i = 0; i<N; ++i) {
        connectedPoints.push_back(false);
    }

    connectedPoints[0] = true;
}

// reset parameters for new run
void PrimParams::refresh_all()
{
    dist = 999999;
    y_diff = 999999;
    x_max = 999999;
    parent = 0;
    child = 0;
}


// set all major parameters at once
void PrimParams::set_all(int d, int y, int x, int p, int c)
{
    dist = d;
    y_diff = y;
    x_max = x;
    parent = p;
    child = c;
}


// check weights and update values if necessary
// used to determine the proper parent/child relationships
void PrimParams::check_params(int d, int y, int x, int p, int c)
{
    if (d < dist) {
        set_all(d, y, x, p, c);
    } else if (d == dist) {
        // 1st check tie. check 2nd (y diff)
        if (y < y_diff) {
            set_all(d, y, x, p, c);
        } else if (y == y_diff) {
            // 1st and 2nd tie. check 3rd (x max)
            if (x < x_max) {
                set_all(d, y, x, p, c);
            }
        }
    }
}

// add node to MST with correct parent/child relationship
void PrimParams::add_node(int *graph_D, int *graph_y, int *graph_x, int N)
{
    int d, y, x, reference, endPoint; // temp parameters

    for (int j = 0; j < N; j++) {
        if (connectedPoints[j]) { //if already connected
            reference = j;
        } else { //if not connected, skip
            continue;
        }

        // we now have the reference. now we need an unconnected node
        for (int k = 0; k < N; k++) {
            if (!connectedPoints[k]) {
                endPoint = k; // we now have an unconnected node
            } else {
                continue; // keep searching
            }

            // nodes found, check parameters 
            d = *(graph_D + N*reference + endPoint); // distance
            y = *(graph_y + N*reference + endPoint); // y difference
            x = *(graph_x + N*reference + endPoint); // x max

            check_params(d, y, x, reference, endPoint); // check against current best weights and update if necessary
        }
    }
}

