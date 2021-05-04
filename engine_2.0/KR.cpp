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

#include "KR.h"
using namespace std;

// Global variables for KR
// ---------------------------------------------------------
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

// Write the results of Kahng/Robins to file for use by plotter
void writeKRResults() {

    ofstream results;
    results.open("KRResults.txt");
    for (int i = 1; i < finalSteinerParents.size(); ++i) {
        results << finalSteinerParents[i] << " " << finalSteinerChildren[i] << "\n";
    }
    results.close();

    results.open("KRWirelengths.txt");

    for (int i = 0; i < wireLengths.size(); ++i) {
        results << wireLengths[i] << endl;
    }

    results.close();

    results.open("KRSteiners.txt");
    for (int i = 0; i < listOfSteiners.size(); ++i) {
        results << listOfSteiners[i] << "\n";
    }
    results.close();
}


// Finds the wirelength of a steiner graph
int steinerWL(int parent_nodes[], int child_nodes[], int gridSize) {
    int currentWL = 0;
    bool duplicate = false;
    for (int i = 0; i < 999999; i++) {
        int currentParent = parent_nodes[i];
        int currentChild = child_nodes[i];
        if (parent_nodes[i] == 0 && parent_nodes[i + 1] == 0) {
            break; // end of array
        }
        for (int j = 0; j < i; j++) { // need to check for duplicate connections
            if ((currentParent == parent_nodes[j] && currentChild == child_nodes[j]) || (currentParent == child_nodes[j] && currentChild == parent_nodes[j])) {
                duplicate = true; // duplicate connection has been found
                break;
            }
        }
        if (duplicate) {
            duplicate = false; // dont count the duplicate's wirelength
            continue;
        } else {
            int parentY = currentParent / (gridSize + 1);
            int parentX = currentParent - (parentY*gridSize + parentY);
            int childY = currentChild / (gridSize + 1);
            int childX = currentChild - (childY*gridSize + childY);
            int wlToAdd = (abs(parentY - childY) + abs(parentX - childX));
            currentWL += wlToAdd;
        }
    }
    return currentWL;
}


// KR Functions
// Converts the MST points to grid numbers
// REFERENCE
// grid number = (y*gridSize + y) + x
// y = grid number / (gridSize + 1)
// x = grid number - (y*gridSize+y)
// (x, y) is the coordinate
int pointsToGridNumbers(int parent_nodes[], int child_nodes[], int *x_pts, int *y_pts, int gridSize, int N) {
    for (int i = 1; i < N; i++) {
        currentGridParents[i] = (*(y_pts + parent_nodes[i]) * gridSize) + *(y_pts + parent_nodes[i]) + *(x_pts + parent_nodes[i]);
        currentGridChildren[i] = (*(y_pts + child_nodes[i]) * gridSize) + *(y_pts + child_nodes[i]) + *(x_pts + child_nodes[i]);
    }
    return 0;
}

// Creates the steiner graph with a given steiner point
// Gain of the new graph is returned
int getGraph(int currentGridParents[], int currentGridChildren[], int potentialSteiner, int gridSize, int indexOfTemp, int bestWL) {
    for (int i = 0; i < 999999; i++) { //STARTS at 1
        if ((currentGridParents[i] == 0) && (currentGridChildren[i] == 0) && currentGridParents[i+1] == 0 && currentGridChildren[i+1] == 0) {
            break;
        }
        int parentY = currentGridParents[i] / (gridSize + 1);
        int parentX = currentGridParents[i] - (parentY*gridSize + parentY);
        int childY = currentGridChildren[i] / (gridSize + 1);
        int childX = currentGridChildren[i] - (childY*gridSize + childY);
        int steinerY = potentialSteiner / (gridSize + 1);
        int steinerX = potentialSteiner - (steinerY*gridSize + steinerY);
        int currentEdge = abs(parentY - childY) + abs(parentX - childX); // calculates distance of current existing edge
        int newEdge1 = abs(parentY - steinerY) + abs(parentX - steinerX); // calculates distance of 1st new edge created by steiner point
        int newEdge2 = abs(childY - steinerY) + abs(childX - steinerX); // calculates distance of 2nd new edge created by steiner point
        if ((newEdge1 < currentEdge) && (newEdge2 < currentEdge)) { // wire length becomes smaller if both new edges are smaller
            tempParents[indexOfTemp] = currentGridParents[i];
            tempChildren[indexOfTemp] = potentialSteiner;
            indexOfTemp++;
            tempParents[indexOfTemp] = potentialSteiner;
            tempChildren[indexOfTemp] = currentGridChildren[i];
            indexOfTemp++;
        } else { // wire length will increase or not change
            tempParents[indexOfTemp] = currentGridParents[i];
            tempChildren[indexOfTemp] = currentGridChildren[i];
            indexOfTemp++;
        }
    }
    int wireLength = steinerWL(tempParents, tempChildren, gridSize); // use newly created edges to calculate the WL
    int gain = bestWL - wireLength; // get the gain
    return gain;
}


// The best steiner point to insert is found, so the current best graph is updated
void updateGrid() {
    bool duplicate = false;
    int index = 0;
    for (int i = 0; i < 999999; i++) {
        int parent = tempParents[i];
        int child = tempChildren[i];
        if ((parent == 0 && child == 0) && (tempParents[i+1] == 0 && tempChildren[i+1] == 0)) {
            break; // end of the array
        }
        for (int j = 0; j < i; j++) {
            if ((parent == tempParents[j] && child == tempChildren[j]) || (parent == tempChildren[j] && child == tempParents[j])) {
                duplicate = true;
                break;
            }
        }
        if (duplicate) {
            duplicate = false;
            continue;
        }
        currentGridParents[index] = parent;
        currentGridChildren[index] = child;
        // need to empty out currentGridParents
        index++;
        duplicate = false;
    }
}

// gets called if a new temporary graph needs to be built from a steiner point
void emptyTemp() {
    for (int i = 0; i < 999999; i++) {
        if ((tempParents[i] == 0 && tempChildren[i] == 0) && (tempParents[i+1] == 0 && tempChildren[i+1] == 0)) {
            break; // end of the array
        } else {
            tempParents[i] = 0;
            tempChildren[i] = 0;
        }
    }
}

// Copies the contents from the current best graph
void updateFinalSteiners() {
    int counter = 0;
    for (int i = finalIndex; i < 999999; i++) {
        if ((currentGridParents[counter] == 0 && currentGridChildren[counter] == 0) && (currentGridParents[counter+1] == 0 && currentGridChildren[counter+1] == 0)) {
            finalIndex = i;
            break; // end of the array
        } else {
            finalSteinerParents.push_back(currentGridParents[counter]);
            finalSteinerChildren.push_back(currentGridChildren[counter]);
            counter++; // counter keeps looping after every iteration
        }
    }
}

void updateArrays() {
    updateGrid();
    updateFinalSteiners();
    emptyTemp();
}

// This function takes in the hanan grid and then runs the rest of the KR algo
void runKR(bool hananGrid[], int gridSize, int N, int bestWL) {
    int bestSteiner = -1;
    int finalGain = 1;
    int bestCurrentGain = -1;
    while (finalGain > 0) { // keep looping until there is no more positive gain to be had
        for (int i = 0; i < gridSize*(gridSize+1)+gridSize + 1; i++) {
            if (!hananGrid[i]) { // this is not a steiner
                continue;
            }
            int potentialSteiner = i; // found a steiner
            int currentGain = getGraph(currentGridParents, currentGridChildren, potentialSteiner, gridSize, indexOfTemp, bestWL);
            if (currentGain > 0 && currentGain > bestCurrentGain) { // we find a new potential best steiner
                bestSteiner = potentialSteiner;
                bestCurrentGain = currentGain;
                emptyTemp();
            }
        } //end of for loop
        if (bestCurrentGain > 0) { // the best gain that is found is positive. So, we call update arrays func
            getGraph(currentGridParents, currentGridChildren, bestSteiner, gridSize, indexOfTemp, bestWL);        
            updateArrays();
            hananGrid[bestSteiner] = false; //no longer a steiner point, so remove it
            finalGain = bestCurrentGain; // this is the gain obtained from the new steiner point
            bestCurrentGain = -1; // reset for next iteration
            bestWL = bestWL - finalGain; // this is the new WL from the new steiner point
            listOfSteiners.push_back(bestSteiner); // add to list of steiners
            wireLengths.push_back(bestWL); // add to list of WLs
        } else {
            emptyTemp();
            break; // exit loop
        }
    }
}



// This function creates the Hanan grid
// An MST point is identified, and then the graph is read from left to right to find all hanan points

// REFERENCE
// grid number = (y*gridSize + y) + x
// y = grid number / (gridSize + 1)
// x = grid number - (y*gridSize+y)
// (x, y) is the coordinate
void hananGrid(int parent_nodes[], int child_nodes[], int *graph_D, int *graph_y, int *graph_x, int *x_pts, int *y_pts, int N, int gridSize, int bestWL) {
    int numberOfPoints = (gridSize*(gridSize+1)+gridSize) + 1;
    bool isSteiner[numberOfPoints]; // the indexes of the array correspond to the grid number. if value == true, then its a steiner point
    for (int i = 0; i < gridSize*(gridSize+1)+gridSize + 1; i++) {
        isSteiner[i] = false;
    }
    for (int i = 0; i < N; i++) {
        int currentPoint = child_nodes[i];
        int currentX = *(x_pts + currentPoint);
        int steinerY = *(y_pts + currentPoint);
        for (int j = 0; j < N; j++) {
            if (j == i) {
                continue;
            }
            int destination = child_nodes[j];
            if (*(x_pts + currentPoint) == *(x_pts + destination) || *(y_pts + currentPoint) == *(y_pts + destination)) {
                continue;
            }
            int steinerX = *(x_pts + destination);
            isSteiner[(steinerY*gridSize + steinerY) + steinerX] = true;
        }
    }
    pointsToGridNumbers(parent_nodes, child_nodes, x_pts, y_pts, gridSize, N);
    runKR(isSteiner, gridSize, N, bestWL);
}

void init_KR(void)
{
    finalSteinerParents.push_back(0);
    finalSteinerChildren.push_back(0);
    return;
}

int get_KR_WL(void)
{
    return wireLengths.back();
}














