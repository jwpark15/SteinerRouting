#ifndef RESULTS_H
#define RESULTS_H

#include <iostream>
#include <fstream>
#include <string.h>
#include <string>
#include <cstdio>

using namespace std;

int calculateWL(int parent_nodes[], int child_nodes[], int *x_pts, int *y_pts, int N);
void parseFile(int *x, int *y, int N, const char filename[]);
int getNumNodes(std::string filepath);
void exportFilename(const char *filename);
void writeWirelengths(int mst, int lrst, int kr);
int getGridSize(std::string filepath);

#endif
