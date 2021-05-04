#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <cstdio>
#include <vector>
#include <stack>
#include <math.h>
#include <algorithm> //for min/max
#include <ctime>
#include <chrono>

// set to 0 to disable debugging prints
#define DEBUG 0
#define VERBOSE_DEBUG 0

using namespace std;

#include "MST.h"
#include "LRST.h"
#include "KR.h"
#include "results.h"




int main(int argc, char** argv)
{
    const char *filename = argv[1];
    std::string temp_filename = argv[1];
   
    int N = getNumNodes(temp_filename);
    exportFilename(filename); // to be used by plotter

    cout << "NUM NODES: " << N << endl;

    // N is number of nodes in graph
    // x and y coordinates for nodes
    int x[N];
    int y[N];
    
    parseFile(x, y, N, filename);
    
    // init graphs for Manhattan Distance, y difference, and x maximum
    int graph_D[N][N];
    int graph_y[N][N];
    int graph_x[N][N];

    int parent_nodes[N];
    int child_nodes[N];

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    auto c_start = high_resolution_clock::now(); // Capture start time of the algos

    cout << endl << "============== Create Weight Graph =================" << endl << endl;
    // create 2D arrays with values for Prim
    createWeightGraphs(*graph_D, *graph_y, *graph_x, x, y, N); 

    // Run Prim's Algorithm
    cout << endl << "============== Generating MST =================" << endl << endl;
    runPrim(*graph_D, *graph_y, *graph_x, parent_nodes, child_nodes, N);
    int initialWL = calculateWL(parent_nodes, child_nodes, x, y, N);
    auto c_endPrim = high_resolution_clock::now(); // Capture end time of Prim

    // Run L-RST Algorithm
    cout << endl << "============== Running L-RST =================" << endl << endl;
    int overlap_LRST = runLRST(*graph_D, *graph_y, *graph_x, x, y, parent_nodes, child_nodes, N);
    auto c_endLRST = high_resolution_clock::now(); // Capture end time of L-RST

    // Run Kahng / Robins Algorithm 
    cout << endl << "============== Running Kahng / Robins =================" << endl << endl;
    int gridSize = getGridSize(temp_filename);
    init_KR();
    hananGrid(parent_nodes, child_nodes, *graph_D, *graph_y, *graph_x, x, y, N, gridSize, initialWL);
    int KR_WL = get_KR_WL();
    writeKRResults();

    auto c_endKR = high_resolution_clock::now(); // Capture end time of KR and all Algos

    // calculate Wirelength
    int LRST_WL = initialWL - overlap_LRST;

    writeWirelengths(initialWL, LRST_WL, KR_WL);

    auto c_end = high_resolution_clock::now();
    duration<double, std::milli> time_elapsed_ms = c_end - c_start;
    duration<double, std::milli> time_elapsed_MST_ms = c_endPrim - c_start;
    duration<double, std::milli> time_elapsed_LRST_ms = c_endLRST - c_endPrim;
    duration<double, std::milli> time_elapsed_KR_ms = c_endKR - c_endLRST;

    cout << "===================================" << endl;
    cout << "Original Wirelength: " << initialWL << endl;
    cout << "Grid Size: " << gridSize << "x" << gridSize << ", Number of Points: " << N << endl;
    cout << "Wirelength after L-RST: " << LRST_WL << endl;
    cout << "Wirelength after Kahng / Robins: " << KR_WL << endl;
    cout << "CPU time used for generating MST: " << time_elapsed_MST_ms.count() << " ms\n";
    cout << "CPU time used for running L-RST: " << time_elapsed_LRST_ms.count() << " ms\n";
    cout << "CPU time used for running Kahng / Robins: " << time_elapsed_KR_ms.count() << " ms\n";
    cout << "CPU time used for all 3 algorithms: " << time_elapsed_ms.count() << " ms\n";
    cout << "===================================" << endl;


    return 0;
}
