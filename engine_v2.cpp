#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cstdio>

using namespace std;
#define SIZE 5


// update graph 2D array with weights of Manhattan Distance
// row and column dimensions should be the same for all arrays
template <size_t rowsD, size_t colsD, size_t rowsy, size_t colsy, size_t rowsx, size_t colsx>
void createWeightGraphs(int (&graph_D)[rowsD][colsD], int (&graph_y)[rowsy][colsy], int (&graph_x)[rowsx][colsx], int *x_pts, int *y_pts)
{
    int x1, x2, y1, y2;
    for (int i=0; i < rowsD; ++i) {
        x1 = *(x_pts + i);
        y1 = *(y_pts + i);
        for (int j=0; j < colsD; ++j) {
            x2 = *(x_pts + j);
            y2 = *(y_pts + j);
            graph_D[i][j] = abs(x1 - x2) + abs(y1 - y2); 
            graph_y[i][j] = -abs(y1 - y2);
            if (x1 > x2) {
                graph_x[i][j] = -1 * x1;
            } else {
                graph_x[i][j] = -1 * x2;
            }
            
            /* Prints for Debugging
            cout << "i: " << i << ", j: " << j << endl;
            cout << graph_D[i][j] << endl;
            cout << graph_y[i][j] << endl;
            cout << graph_x[i][j] << endl;
            */
        } 
    }
    return;
}
            

template <size_t num1, size_t num2>
int findMinKeyIndex(int (&key_arr)[num1], bool (&connected_arr)[num2])
{
    int min_val = 99999999;
    int min_index;

    for (int i = 0; i < num1; ++i) {
        if (connected_arr[i] == false && key_arr[i] < min_val) {
            min_val = key_arr[i];
            min_index = i;
        }
    }
    /* Print for Debugging
    cout << min_index << endl;
    */
    return min_index;
}

template <size_t rowsD, size_t colsD, size_t rowsy, size_t colsy, size_t rowsx, size_t colsx>
void runPrim(int (&graph_D)[rowsD][colsD], int (&graph_y)[rowsy][colsy], int (&graph_x)[rowsx][colsx])
{
    // init MST arrays 
    bool connected[rowsD] = {false};
    int node_keys[rowsD];
    for (int x = 0; x < rowsD; ++x)
        node_keys[x] = 99999999;
    node_keys[0] = 0;
    int node_parents[rowsD];
    node_parents[0] = -1;
    
    int min_index, dist;
    for (int i = 0; i < rowsD; ++i) {
        min_index = findMinKeyIndex(node_keys, connected);
        
        // add node to tree
        connected[min_index] = true;
        for (int j = 0; j < colsD; ++j) {
            dist = graph_D[min_index][j]; 
            if (dist && !connected[j]) {
                if (dist < node_keys[j]) {
                    cout << "Min index: " << min_index << endl;
                    cout << "Dist: " << dist << endl;
                    node_parents[j] = min_index;
                    node_keys[j] = dist;
                } else if (dist == node_keys[j]) {
                    // add additional checks here
                }
            } 
        }
    }
    for (int i = 0; i < rowsD; ++i) {
        cout << "Parent: " << node_parents[i] << ", Child: " << i << endl;
    }
    return;
}


template <size_t N1, size_t N2>
void parseFile(int (&x)[N1], int (&y)[N2], const char filename[]) 
{
    ifstream benchmark(filename, ios::in);
    std::string word;
    for (int i = 0; i < N1; ++i) {
        benchmark >> word;
        x[i] = stoi(word);
        benchmark >> word;
        y[i] = stoi(word);
    }
    return;
}

int main()
{
    // TODO get filename from command line parsing
    const char filename[] = "Points/points_10_5.pts";

    // TODO get number of nodes using strtok from filename
    // N is number of nodes in graph
    const int N = SIZE;
    int x[N];
    int y[N];
    
    parseFile(x, y, filename);
    /* Print for Debugging 
    for (int i = 0; i < N; ++i) {
        cout << x[i] << ", " << y[i] << endl;
    }
    */

    // init graphs for Manhattan Distance, y difference, and x maximum
    int graph_D[N][N];
    int graph_y[N][N];
    int graph_x[N][N];

    // create 2D arrays with values for Prim
    createWeightGraphs(graph_D, graph_y, graph_x, &x[0], &y[0]); 
    
    // Run Prim's Algorithm
    runPrim(graph_D, graph_y, graph_x);

    return 0;
}
