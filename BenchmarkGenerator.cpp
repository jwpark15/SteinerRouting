#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <cstdio>
#include <ctime>
#include <time.h>
#include <vector>
#include <cstdlib>
using namespace std;


int main() {
    srand(time(NULL));
    int uniquePoints = 0;
    int goalNumberOfPoints  = (rand() % (200 + 1 - 100)) + 100;
    int gridSize = goalNumberOfPoints * 0.3;

    int test = 5;
    vector<int> initialX;
    vector<int> initialY;
    vector<int> finalX;
    vector<int> finalY;
    bool duplicate = false;

    while (uniquePoints < goalNumberOfPoints) {
        int randomX = rand() % gridSize + 1;
        int randomY = rand() % gridSize + 1;

        // need to see if these are within range and are not duplicate

        if (uniquePoints == 0 && randomX <= gridSize && randomY <= gridSize) {
            finalX.push_back(randomX);
            finalY.push_back(randomY);
            uniquePoints++;
            //uniquePoints = 999999;
            continue;
        }

        for (int j = 0; j < sizeof(finalX); j++) {
            if (finalX[j] == randomX && finalY[j] == randomY) {
                duplicate = true;
                break;
            }
        }

        if (duplicate || randomX > gridSize || randomY > gridSize) {
            // This is a bad point. Do not use it
            duplicate = false;
            continue;
        } else {
            // This is a non-duplicate that is within range of the grid
            // Add it to the final file and increment the number of points in the file
            finalX.push_back(randomX);
            finalY.push_back(randomY);
            uniquePoints++;
            duplicate = false;
        }

    }

    ofstream results;
    results.open("points_" + to_string(gridSize) + "_" + to_string(uniquePoints) + ".pts");
    for (int i = 0; i < uniquePoints; i++) {
        if (i == uniquePoints - 1) {
            results << to_string(finalX[i]) << " " << to_string(finalY[i]);
        } else {
            results << to_string(finalX[i]) << " " << to_string(finalY[i]) << endl;
        }
    }

    results.close();



    return 0;
}
