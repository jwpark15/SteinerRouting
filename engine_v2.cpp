#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <cstdio>

using namespace std;
#define SIZE 5
#define DEBUG 1

// update graph 2D array with weights of Manhattan Distance
// row and column dimensions should be the same for all arrays
void createWeightGraphs(int *graph_D, int *graph_y, int *graph_x, int *x_pts, int *y_pts, int N)
{
    int x1, x2, y1, y2;
    int ptr_val = 0;
    for (int i=0; i < N; ++i) {
        x1 = *(x_pts + i);
        y1 = *(y_pts + i);
        for (int j=0; j < N; ++j) {
            ptr_val = i*N + j; 
            x2 = *(x_pts + j);
            y2 = *(y_pts + j);
            *(graph_D + ptr_val) = abs(x1 - x2) + abs(y1 - y2); 
            *(graph_y + ptr_val) = -abs(y1 - y2);
            if (x1 > x2) {
                *(graph_x + ptr_val) = -1 * x1;
            } else {
                *(graph_x + ptr_val) = -1 * x2;
            }
            
            // Print for debug 
            if (DEBUG) {
                cout << "i: " << i << ", j: " << j << endl;
                cout << *(graph_D + ptr_val) << endl;
                cout << *(graph_y + ptr_val) << endl;
                cout << *(graph_x + ptr_val) << endl;
            }
        } 
    }
    return;
}
            

int findMinKeyIndex(int *key_arr, bool *connected_arr, int N)
{
    int min_val = 99999999;
    int min_index;

    for (int i = 0; i < N; ++i) {
        if (*(connected_arr + i) == false && *(key_arr + i) < min_val) {
            min_val = *(key_arr + i);
            min_index = i;
        }
    }
    // Print for Debug
    if (DEBUG)
        cout << "MIN INDEX: " << min_index << endl;

    return min_index;
}

void writePrimResults(int *parent_nodes, int N)
{
    ofstream results;
    results.open("primResults.txt");
    for (int i = 0; i < N; ++i) {
        results << *(parent_nodes + i) << " " << i << "\n";
    }
    results.close();
    return;
}

void runPrim(int *graph_D, int *graph_y, int *graph_x, int N)
{
    // init MST arrays 
    bool connected[N];
    int node_keys[N];
    int node_keys_y[N];
    int node_keys_x[N];
    for (int x = 0; x < N; ++x) {
        connected[x] = false;
        node_keys[x] = 99999999;
        node_keys_y[x] = 99999999;
        node_keys_x[x] = 99999999;
    }
    node_keys[0] = 0;
    int node_parents[N];
    node_parents[0] = -1;
    
    int min_index, dist, y_diff, x_max;
    for (int i = 0; i < N; ++i) {
        min_index = findMinKeyIndex(node_keys, connected, N);
        
        // add node to tree
        connected[min_index] = true;
        for (int j = 0; j < N; ++j) {
            dist = *(graph_D + N*min_index + j); 
            if (dist && !connected[j]) {
                y_diff = *(graph_y + N*min_index + j); 
                x_max = *(graph_x + N*min_index + j); 
                if (dist < node_keys[j]) {
                    //cout << "Min index: " << min_index << endl;
                    //cout << "Dist: " << dist << endl;
                    node_parents[j] = min_index;
                    node_keys[j] = dist;
                    node_keys_y[j] = y_diff;
                    node_keys_x[j] = x_max;
                } else if (dist == node_keys[j]) {
                    if (y_diff < node_keys_y[j]) {
                        node_parents[j] = min_index;
                        node_keys[j] = dist;
                        node_keys_y[j] = y_diff;
                        node_keys_x[j] = x_max;
                    } else if (y_diff == node_keys_y[j]) {
                        if (x_max < node_keys_x[j]) {
                            node_parents[j] = min_index;
                            node_keys[j] = dist;
                            node_keys_y[j] = y_diff;
                            node_keys_x[j] = x_max;
                        } 
                    }
                }
            } 
        }
    }
    
    // Print for debug
    if (DEBUG) {
        for (int i = 0; i < N; ++i) {
            cout << "Parent: " << node_parents[i] << ", Child: " << i << endl;
        }
    }

    writePrimResults(node_parents, N);
    return;
}


void parseFile(int *x, int *y, int N, const char filename[]) 
{
    ifstream benchmark(filename, ios::in);
    std::string word;
    for (int i = 0; i < N; ++i) {
        benchmark >> word;
        *(x + i) = stoi(word);
        benchmark >> word;
        *(y + i) = stoi(word);
    }
    return;
}

int getNumNodes(std::string filepath)
{
    cout << "FILE: " << filepath << endl;
    int index1 = filepath.find_last_of('_');
    int index2 = filepath.find_last_of('.');
    std::string temp_str = filepath.substr(index1 + 1, index2 - index1 - 1);
    int num_nodes = stoi(temp_str);

    //cout << num_nodes << endl;
    
    return num_nodes;
}

int main(int argc, char** argv)
{
    cout << argv[1] << endl;
    const char *filename = argv[1];
    std::string temp_filename = argv[1];
   
    int N = getNumNodes(temp_filename);

    cout << "NUM NODES: " << N << endl;
    // TODO get number of nodes using strtok from filename
    // N is number of nodes in graph
    //const int N = SIZE;
    int x[N];
    int y[N];
    
    parseFile(x, y, N, filename);
    
    // Print for debug
    if (DEBUG) {
        for (int i = 0; i < N; ++i) {
            cout << x[i] << ", " << y[i] << endl;
        }
    } 

    // init graphs for Manhattan Distance, y difference, and x maximum
    int graph_D[N][N];
    int graph_y[N][N];
    int graph_x[N][N];

    // create 2D arrays with values for Prim
    createWeightGraphs(*graph_D, *graph_y, *graph_x, x, y, N); 
    // Run Prim's Algorithm
    runPrim(*graph_D, *graph_y, *graph_x, N);
    return 0;
}
