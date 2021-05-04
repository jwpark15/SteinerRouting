#ifndef KAHNGROBINS_H
#define KAHNGROBINS_H

#include <vector>

// Global variables for KR
/* ---------------------------------------------------------
std::vector<int> finalSteinerParents;
std::vector<int> finalSteinerChildren;

std::vector<int> wireLengths; // Keeps track of wirelengths per iteration
std::vector<int> listOfSteiners; // Keeps track of steiner points added per iteration
// these 2 arrays will have the parents and children created from a possible steiner point
int tempParents[9999999];
int tempChildren[9999999];

// These 2 arrays will have the current best grid
// This current best grid will be used as a reference to compare with potential graphs
int currentGridParents[9999999];
int currentGridChildren[9999999];

int finalIndex = 1;
int indexOfTemp = 0;
// ---------------------------------------------------------
*/


// KR Functions
void init_KR(void);
void writeKRResults();
int steinerWL(int parent_nodes[], int child_nodes[], int gridSize);
int pointsToGridNumbers(int parent_nodes[], int child_nodes[], int *x_pts, int *y_pts, int gridSize, int N);
int getGraph(int currentGridParents[], int currentGridChildren[], int potentialSteiner, int gridSize, int indexOfTemp, int bestWL);
void updateGrid();
void emptyTemp();
void updateFinalSteiners();
void updateArrays();
void runKR(bool hananGrid[], int gridSize, int N, int bestWL);
void hananGrid(int parent_nodes[], int child_nodes[], int *graph_D, int *graph_y, int *graph_x, int *x_pts, int *y_pts, int N, int gridSize, int bestWL);
int get_KR_WL(void);




#endif
