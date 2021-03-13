#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cstdio>

using namespace std;


int main() {
    
    ifstream benchmark("points_30_100.txt", ios::in); //open file
    string sourceX;
    string sourceY;
    int inputSize; // grid size
    int inputTotal; // number of points

    int test[100] = {0};

    test[1] = 0;
    fill(begin(test), end(test), 0);

    // not sure how to read grid size and nodes from the file size
    // manual input for now
    cout << "Enter gride size: ";
    cin >> inputSize;

    cout << "Enter number of nodes: ";
    cin >> inputTotal;

    int size = inputSize;
    int total = inputTotal;

    // first line is always the source point
    benchmark >> sourceX;
    benchmark >> sourceY;

    // making arrays to keep store coordinates with int data types
    int xCoordinates [total];
    int yCoordinates [total];

    // storing source point's data
    xCoordinates[0] = stoi(sourceX);
    yCoordinates[0] = stoi(sourceY);

    // stores current coordinate
    string currentX;
    string currentY;

    // this loop takes in the coordinates from txt file and stores into int array
    for (int k = 1; k < total; k++) {
        benchmark >> currentX;
        benchmark >> currentY;

        xCoordinates[k] = stoi(currentX);
        yCoordinates[k] = stoi(currentY);
    }

    // Manhattan Distance Loop

    // array that stores manhattan distance from source node to another node
    int distances [total];

    // 2nd value in prim formula
    int yDistances [total];

    // 3rd value in prim formula
    int xDistances [total];

    bool connected[100] = {false}; //stores connected REFERENCE points
    connected[0] = true; //this is the source
    int smallestManhattan = 2*size + 1; //largest theoretical manhattan

    // NEED TO DYNAMICALLY ALLOCATE THESE
    // Reference = point A
    // Potentual = point B
    int manhattanTiesReference[100] = {0};
    int manhattanTiesPotential[100] = {0};

    // referencesToConnect[0] and potentialsToConnect[0] are the starting and ending points of line 1
    // referencesToConnect[1] and potentialsToConnect[1] are the starting and ending points of line 2
    int referencesToConnect[100] = {0};
    int potentialsToConnect[100] = {0};

    int numberOfTies = 0; //keeps track of manhattan ties
    int numberOfConnectedPoints = 0;
    int closestReference;
    int closestPotential;

    // This outer for loop and the code within will find the smallest manhattan distance between 2 points that need to be connected
    // The outer for loop finds an already connected point and treats that as a reference.
    // The inner for loop then finds an unconnected point (potential) and calculates the manhattan distance to the current reference point
    for (int i = 0; i < total; i++) {
        if (!connected[i]) { // look for points that havent been connected
            continue;
        }

        int referencePointX = xCoordinates[i];
        int referencePointY = yCoordinates[i];


        for (int j = 0; j < total; j++) {

            if (connected[j]) { // look for points that need to be connected
                continue;
            }

            int potentialPointX = xCoordinates[j];
            int potentialPointY = yCoordinates[j];

            int manhattanDistance = abs(referencePointX - potentialPointX) + abs(referencePointY - potentialPointY); //manhattan formula
            //printf("%d \n", manhattanDistance); used foe debugging

            if (manhattanDistance < smallestManhattan) { // if there is a new small distance
                smallestManhattan = manhattanDistance; // update the new smallest

                // any previous manhattan ties are voided
                fill(begin(manhattanTiesReference), end(manhattanTiesReference), 0);
                fill(begin(manhattanTiesPotential), end(manhattanTiesPotential), 0);

                closestReference = i; //i is from outer loop
                closestPotential = j;
                numberOfTies = 0; // reset
            } else if (manhattanDistance == smallestManhattan) { //if there is a tie

                //record points a and b
                manhattanTiesReference[numberOfTies] = i;
                manhattanTiesPotential[numberOfTies] = j;
                numberOfTies++;
            }
        }
    }

    //if there are no ties
    if (numberOfTies == 0 && numberOfConnectedPoints < size) {
        referencesToConnect[numberOfConnectedPoints] = closestReference + 1;
        potentialsToConnect[numberOfConnectedPoints] = closestPotential + 1;
        connected[closestReference] = true;
        numberOfConnectedPoints++;
    }

    //used for debugging
    //printf("%d", referencesToConnect[0]);
    //printf("%d", potentialsToConnect[0]);

    //---------------------------------------------------------------------------------

    // TO-DO:
    // 1. Make lines 105 - 144 create starting and ending points for more than 1 line.
    // 2. Handle ties in manhattan distance
    // 3. Dynamically allocate the arrays from lines 77, 78, 82, and 83
    
    

    return 0;
    
}