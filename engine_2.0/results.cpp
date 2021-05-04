#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <cstdio>
#include <math.h>
#include "results.h"

using namespace std;
// calculate the Manhattan distance wirelength for the MST
int calculateWL(int parent_nodes[], int child_nodes[], int *x_pts, int *y_pts, int N) {
    int wireLength = 0; 
    int pointA = 0; 
    int pointB = 0; 
    int distance = 0; 

    for (int i = 1; i < N; i++) {
        pointA = parent_nodes[i];
        pointB = child_nodes[i];
        distance = abs(*(x_pts + pointA) - *(x_pts + pointB)) + abs(*(y_pts + pointA) - *(y_pts + pointB)); 
        wireLength += distance;
    }    

    return wireLength;
}

// parse input file to get x and y coordinates for each node
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

// get number of nodes in benchmark by parsing filename
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


// write filename of benchmark to file to be used by python plotter
void exportFilename(const char *filename)
{
    ofstream results;
    results.open("benchmark.txt");
    results << filename;
    results.close();
    return;

}


// write wirelength results to file
void writeWirelengths(int mst, int lrst, int kr)
{
    ofstream results;
    results.open("wirelengthResults.txt");
    results << "MST: " << mst << endl;
    results << "LRST: " << lrst << endl;
    results << "KR: " << kr << endl;
    results.close();
    return;

}

int getGridSize(std::string filepath) {
    int index1 = filepath.find_last_of('_');
    int index2 = filepath.find_last_of('.');
    std::string temp_str = filepath.substr(index1 + 1, index2 - index1 - 1);

    int index3 = filepath.find_first_of('_');
    int index4 = filepath.find_last_of('_');
    temp_str = filepath.substr(index3 + 1, index4 - index3 - 1);

    int gridSize = stoi(temp_str);

    return gridSize;
}


