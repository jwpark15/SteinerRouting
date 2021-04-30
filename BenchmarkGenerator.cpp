#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <cstdio>
#include <ctime>
#include <time.h>
#include <vector>

using namespace std;


int main() {
    srand(time(NULL));
    int goalNumberOfPoints = 3000;
    int uniquePoints = 0;
    int gridSize = 100;

    /*
    ofstream results;
    results.open("points_" + to_string(gridSize) + "_" + to_string(numberOfPoints) + ".pts");

    for (int i = 0; i < numberOfPoints; i++) {
        int randomX = rand() % gridSize + 1;
        int randomY = rand() % gridSize + 1;
        if (i == numberOfPoints - 1) {
            results << to_string(randomX) << " " << to_string(randomY);
        } else {
            results << to_string(randomX) << " " << to_string(randomY) << endl;
        }
    }

    results.close();
    */
    int test = 5;
    vector<int> initialX;
    vector<int> initialY;
    vector<int> finalX;
    vector<int> finalY;
    bool duplicate = false;

    while (uniquePoints < goalNumberOfPoints) {
        cout << "New Point" << endl;
        int randomX = rand() % gridSize + 1;
        int randomY = rand() % gridSize + 1;

        cout << "X: " << randomX << ", Y: " << randomY << endl;
        cout << "Number of points: " << uniquePoints << endl;

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
        cout << "Here" << endl;

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

    cout << uniquePoints << endl;

    for (int i = 0; i < uniquePoints; i++) {
        cout << "I: " << i << ", X: " << finalX[i] << ", Y: " << finalY[i] << endl;
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