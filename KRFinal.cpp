#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <cstdio>
#include <ctime>
#include <vector>

using namespace std;
#define SIZE 5
#define DEBUG 0

// GLOBAL VARIABLES TO ADD
//-------------------------------------------------------------------------
// Vectors that contain the final connections
vector<int> finalSteinerParents;
vector<int> finalSteinerChildren;

vector<int> wireLengths; // Keeps track of wirelengths per iteration
vector<int> listOfSteiners; // Keeps track of steiner points added per iteration

// these 2 arrays will have the parents and children created from a possible steiner point
int tempParents[9999999];
int tempChildren[9999999];

// These 2 arrays will have the current best grid
// This current best grid will be used as a reference to compare with potential graphs
int currentGridParents[9999999];
int currentGridChildren[9999999];

int finalIndex = 1;
int indexOfTemp = 0;
//-------------------------------------------------------------------------

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

int findMinKeyIndexY(int *key_arr_y, bool *connected_arr, int N)
{
    int min_val_y = 99999999;
    int min_index_y;

    for (int i = 0; i < N; ++i) {

        if (*(connected_arr + i) == false && *(key_arr_y + i) < min_val_y) {
            min_val_y = *(key_arr_y + i);
            min_index_y = i;
        }
    }
    // Print for Debug
    if (DEBUG)
        cout << "MIN INDEX Y: " << min_index_y << endl;

    return min_index_y;
}

int findMinKeyIndexX(int *key_arr_x, bool *connected_arr, int N)
{
    int min_val_x = 99999999;
    int min_index_x;

    for (int i = 0; i < N; ++i) {

        if (*(connected_arr + i) == false && *(key_arr_x + i) < min_val_x) {
            min_val_x = *(key_arr_x + i);
            min_index_x = i;
        }
    }
    // Print for Debug
    if (DEBUG)
        cout << "MIN INDEX X: " << min_index_x << endl;

    return min_index_x;
}

void writePrimResults(int *parent_nodes, int *child_nodes, int N)
{
    ofstream results;
    results.open("primResults.txt");
    for (int i = 0; i < N; ++i) {
        results << *(parent_nodes + i) << " " << *(child_nodes + i) << "\n";
    }
    results.close();
    return;
}

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

void runPrim(int *graph_D, int *graph_y, int *graph_x, int parent_nodes[], int child_nodes[], int N)
{
    // init MST arrays
    // make inits into 1 line
    bool connected[N];
    int node_keys[N];
    int node_keys_y[N];
    int node_keys_x[N];
    bool connectedPoints[N];
    connectedPoints[0] = true;
    
    for (int x = 0; x < N; ++x) {
        connected[x] = false;
        node_keys[x] = 99999999;
        node_keys_y[x] = 99999999;
        node_keys_x[x] = 99999999;
    }
    node_keys[0] = 0;
    int node_parents[N];
    node_parents[0] = -1;
    
    int min_index, min_index_y, min_index_x, dist, y_diff, x_max;
    
    for (int k = 0; k < N; k++) {
        connectedPoints[k] = false;
    }

    connectedPoints[0] = true;

    int numberConnected = 1;
    int reference = 12;
    int endPoint = 9;
    int smallestDistance = 99999;
    int smallestY = 99999;
    int smallestX = 99999;
    int currentDistance = 0;
    int totalEquivalentDistances = 0;
    int locationOfDistanceTies[N];
    int parentNode = 0;
    int childNode = 0;

    parent_nodes[0] = -1;
    child_nodes[0] = 0;
    while (numberConnected != N) {

        // find a point that is connected

        for (int j = 0; j < N; j++) {
            if (connectedPoints[j]) { //if already connected
                reference = j;
            } else { //if not connected, skip
                continue;
            }

            // we now have the reference. now we need an unconnected node
            for (int k = 0; k < N; k++) {
                //cout << "connectedPoints: " << connectedPoints[k] << endl;
                if (!connectedPoints[k]) {
                    endPoint = k; // we now have an unconnected node
                } else {
                    continue; // keep searching
                }

                // now we have a connected node and an unconnected node
                // we can now calculate distance, y diff, and x max

                dist = *(graph_D + N*reference + endPoint); // distance
                y_diff = *(graph_y + N*reference + endPoint); // y difference
                x_max = *(graph_x + N*reference + endPoint); // x max

                // we need to see if distance is tied

                if (dist < smallestDistance) { //all new record
                    // make helper func
                    smallestDistance = dist;
                    smallestY = y_diff;
                    smallestX = x_max;
                    parentNode = reference;
                    childNode = endPoint;
                } else if (dist == smallestDistance) {
                    // manhattan distances are tied
                    // now check the y differences

                    if (y_diff < smallestY) { //if there is a record in y differences
                        // new closest point
                        smallestDistance = dist;
                        smallestY = y_diff;
                        smallestX = x_max;
                        parentNode = reference;
                        childNode = endPoint;
                    } else if (y_diff == smallestY) {
                        // manhattan distances and y distances are tied
                        // now check the x

                        if (x_max < smallestX) { // if there is a record in x max
                            smallestDistance = dist;
                            smallestY = y_diff;
                            smallestX = x_max;
                            parentNode = reference;
                            childNode = endPoint;
                        } else if (x_max == smallestX) {
                            int random = rand() % 100;
                            if (random < 50) {
                                smallestDistance = dist;
                                smallestY = y_diff;
                                smallestX = x_max;
                                parentNode = reference;
                                childNode = endPoint;
                            } else {
                                smallestDistance = dist;
                                smallestY = y_diff;
                                parentNode = reference;
                                childNode = endPoint;
                            }
                        }

                    }
                }


            }
        }

        //cout << "reference: " << reference << endl;

        connectedPoints[childNode] = true;
        //cout << "this point is now connected: " << endPoint << endl;
        //cout << "Parent: " << parentNode << endl;
        //cout << "Child: " << childNode << endl;
        smallestDistance = 99999;
        smallestY = 99999;
        smallestX = 99999;
        //cout << "number connected: " << numberConnected << endl;

        parent_nodes[numberConnected] = parentNode;
        child_nodes[numberConnected] = childNode;
        numberConnected++;

    }

    writePrimResults(parent_nodes, child_nodes, N);


    for (int t = 0 ; t < N; t++) {
        //cout << connectedPoints[t] << endl;
    }
    
    // Print for debug
    if (DEBUG) {
        for (int i = 0; i < N; ++i) {
            cout << "Parent: " << parent_nodes[i] << ", Child: " << child_nodes[i] << endl;
        }
    }

    //writePrimResults(node_parents, N);
    return;
}

int calculateWL(int parent_nodes[], int child_nodes[], int *graph_D, int *graph_y, int *graph_x, int *x_pts, int *y_pts, int N) {
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

    cout << "Length: " << wireLength << endl;

    return wireLength;
}

// Finds the wirelength of a steiner graph
int steinerWL(int parent_nodes[], int child_nodes[], int gridSize) {
    int currentWL = 0;
    bool duplicate = false;
    for (int i = 0; i < 999999; i++) {
        int currentParent = parent_nodes[i];
        int currentChild = child_nodes[i];
        if (parent_nodes[i] == 0 && parent_nodes[i + 1] == 0) {
            break;
        }
        for (int j = 0; j < i; j++) {
            if ((currentParent == parent_nodes[j] && currentChild == child_nodes[j]) || (currentParent == child_nodes[j] && currentChild == parent_nodes[j])) {
                duplicate = true;
                break;
            }
        }
        if (duplicate) {
            duplicate = false;
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
        int currentEdge = abs(parentY - childY) + abs(parentX - childX);
        int newEdge1 = abs(parentY - steinerY) + abs(parentX - steinerX);
        int newEdge2 = abs(childY - steinerY) + abs(childX - steinerX);
        if ((newEdge1 < currentEdge) && (newEdge2 < currentEdge)) {
            tempParents[indexOfTemp] = currentGridParents[i];
            tempChildren[indexOfTemp] = potentialSteiner;
            indexOfTemp++;
            tempParents[indexOfTemp] = potentialSteiner;
            tempChildren[indexOfTemp] = currentGridChildren[i];
            indexOfTemp++;
        } else {
            tempParents[indexOfTemp] = currentGridParents[i];
            tempChildren[indexOfTemp] = currentGridChildren[i];
            indexOfTemp++;
        }
    }
    int wireLength = steinerWL(tempParents, tempChildren, gridSize);
    int gain = bestWL - wireLength;
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
    int new_parents[N]; // need new parents because steiner points will become parents
    for (int i = 0; i < N; i++) {
        new_parents[i] = parent_nodes[i];
    }
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

int main(int argc, char** argv)
{
    //algorithm
    cout << argv[1] << endl;
    const char *filename = argv[1];
    std::string temp_filename = argv[1];

    cout << "NAME: " << temp_filename << endl;
   
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

    int parent_nodes[N];
    int child_nodes[N];

    // create 2D arrays with values for Prim
    createWeightGraphs(*graph_D, *graph_y, *graph_x, x, y, N);

    // Run Prim's Algorithm
    // Get MST
    runPrim(*graph_D, *graph_y, *graph_x, parent_nodes, child_nodes, N);

    // CODE TO ADD
    //-----------------------------------------------------------------------
    // beginning of KR
    finalSteinerParents.push_back(0);
    finalSteinerChildren.push_back(0);
    int gridSize = getGridSize(temp_filename);
    int initialWL = calculateWL(parent_nodes, child_nodes, *graph_D, *graph_y, *graph_x, x, y, N);
    hananGrid(parent_nodes, child_nodes, *graph_D, *graph_y, *graph_x, x, y, N, gridSize, initialWL);
    writeKRResults();
    //------------------------------------------------------------------------

    // Add Global Variables as well
    
    return 0;
}
