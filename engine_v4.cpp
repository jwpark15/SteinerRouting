#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <cstdio>
#include <vector>
#include <stack>

using namespace std;
#define SIZE 5
#define DEBUG 1

int totalOverlap = 0;

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
                cout << "<" << *(graph_D + ptr_val) << ", " << *(graph_y + ptr_val) << ", " << *(graph_x + ptr_val) << ">" << endl;
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

class PrimParams {
    public:
        PrimParams(int N);
        int get_dist()   {return dist;}
        int get_y_diff() {return y_diff;}
        int get_x_max()  {return x_max;}
        int get_parent() {return parent;}
        int get_child()  {return child;} 

        void set_dist(int d)   {dist = d;}
        void set_y_diff(int y) {y_diff = y;}
        void set_x_max(int x)  {x_max = x;}
        void set_parent(int p) {parent = p;}
        void set_child(int c)  {child = c;}
        void set_connected(int index) {connectedPoints[index] = true;}
        
        void refresh_all();
        void set_all(int d, int y, int x, int p, int c);
        void check_params(int d, int y, int x, int p, int c);
        void add_node(int *graph_D, int *graph_y, int *graph_x, int N);
    private:
        int dist;
        int y_diff;
        int x_max;
        int parent;
        int child;
        vector<bool> connectedPoints;
};

PrimParams::PrimParams(int N)
{
    for(int i = 0; i<N; ++i) {
        connectedPoints.push_back(false);
    }

    connectedPoints[0] = true;
}

void PrimParams::refresh_all()
{
    dist = 999999;
    y_diff = 999999;
    x_max = 999999;
    parent = 0;
    child = 0;
}


void PrimParams::set_all(int d, int y, int x, int p, int c) 
{
    dist = d;
    y_diff = y;
    x_max = x;
    parent = p;
    child = c;
}

void PrimParams::check_params(int d, int y, int x, int p, int c)
{
    if (d < dist) {
        set_all(d, y, x, p, c);
    } else if (d == dist) {
        // 1st check tie. check 2nd (y diff)
        if (y < y_diff) {
            set_all(d, y, x, p, c);
        } else if (y == y_diff) {
            // 1st and 2nd tie. check 3rd (x max)
            if (x < x_max) {
                set_all(d, y, x, p, c);
            }
        }
    }
}

void PrimParams::add_node(int *graph_D, int *graph_y, int *graph_x, int N)
{
    int d, y, x;

    int reference = 12;
    int endPoint = 9;

    for (int j = 0; j < N; j++) {
        if (connectedPoints[j]) { //if already connected
            reference = j;
        } else { //if not connected, skip
            continue;
        }

        // we now have the reference. now we need an unconnected node
        for (int k = 0; k < N; k++) {
            if (!connectedPoints[k]) {
                endPoint = k; // we now have an unconnected node
            } else {
                continue; // keep searching
            }

            // nodes found, check parameters 
            d = *(graph_D + N*reference + endPoint); // distance
            y = *(graph_y + N*reference + endPoint); // y difference
            x = *(graph_x + N*reference + endPoint); // x max

            check_params(d, y, x, reference, endPoint);
        }
    }
}

void runPrim(int *graph_D, int *graph_y, int *graph_x, int parent_nodes[], int child_nodes[], int N)
{
    // initialize MST parameters
    PrimParams params(N);
    params.refresh_all();
     
    int numberConnected = 1;
    parent_nodes[0] = -1;
    child_nodes[0] = 0;
    while (numberConnected != N) {
        // add node to tree
        params.add_node(graph_D, graph_y, graph_x, N);

        params.set_connected(params.get_child());

        cout << "Parent: " << params.get_parent() << endl;
        cout << "Child: " << params.get_child() << endl;
        cout << "number connected: " << numberConnected << endl;

        parent_nodes[numberConnected] = params.get_parent();
        child_nodes[numberConnected] = params.get_child();
        numberConnected++;

        params.refresh_all();
    }

    writePrimResults(parent_nodes, child_nodes, N);
    
    // Print for debug
    if (DEBUG) {
        for (int i = 0; i < N; ++i) {
            cout << "Parent: " << parent_nodes[i] << ", Child: " << child_nodes[i] << endl;
        }
    }

    return;
}

class Node {
    public:
        Node(int val, int x_coord, int y_coord);
        // basic setters/getters
        void set_upper_Z(int z) {upper_Z = z;}
        void set_lower_Z(int z) {lower_Z = z;}
        void set_parent(int p) {parent = p;}
        void set_parent_x(int px) {parent_x = px;}
        void set_parent_y(int py) {parent_y = py;}
        void add_child(int c) {children.push_back(c);}
        void add_child_x(int cx) {children_x.push_back(cx);}
        void add_child_y(int cy) {children_y.push_back(cy);}
        void set_x(int coord) {x = coord;}
        void set_y(int coord) {y = coord;}
        int get_x(void) {return x;}
        int get_y(void) {return y;}
        
        // funcs for tree traversal 
        int pick_next_child(void);
        void partialLRST(void);
        void calculate_upperL_overlap(void);
        void calculate_lowerL_overlap(void);

        //debugging
        void print_children();
    private:
        int tag;
        int x;
        int y;
        int lower_Z;
        int upper_Z;
        std::vector<int> upper; //vector of orientations for children if node is upper
        std::vector<int> lower; //vector of orientations for children if node is lower
        std::vector<int> children;
        std::vector<int> children_x;
        std::vector<int> children_y;
        int parent;
        int parent_x;
        int parent_y;
        int child_pointer; // keep track of which branches have been traversed
};

void Node::print_children()
{
    cout << "printing children: " << endl;
    for(auto it = begin(children); it != end(children); ++it)
        cout << *it << endl;
}
  
Node::Node(int val, int x_coord, int y_coord)
{
    tag = val;
    lower_Z = 0;
    upper_Z = 0;
    x = x_coord;
    y = y_coord;
    child_pointer = 0;
}

// return -1 if bottom node or all children have been traversed
// else return tag of next child
int Node::pick_next_child(void)
{
    if (children.size() == 0 || child_pointer >= children.size()) {
        return -1;
    } else {
        int temp_pointer = child_pointer++;
        return children[temp_pointer]; 
    }
}


void Node::calculate_upperL_overlap(void)
{
    int corner_x = (parent_y > y) ? x : parent_x;
    int corner_y = (parent_y > y) ? parent_y : y;
}

void Node::calculate_lowerL_overlap(void)
{
    int corner_x = (parent_y < y) ? x : parent_x;
    int corner_y = (parent_y < y) ? parent_y : y;
}


void Node::partialLRST(void)
{
    // skip bottom nodes
    if (children.size() != 0) {
        // upper L
        calculate_upperL_overlap();

        // lower L
        calculate_lowerL_overlap();
    }

}


void createTree(vector<Node> &tree, int parent_nodes[], int child_nodes[], int x[], int y[], int N)
{
    // each vector index is a Node object
    // index corresponds to node ID
    for (int i = 0; i<N; ++i) {
        Node node(i, x[i], y[i]);
        tree.push_back(node);
    }

    // add parents and children to Node objects
    int p, c;
    for (int j = 0; j<N; ++j) {
        p = parent_nodes[j];
        c = child_nodes[j];
        tree[c].set_parent(p);
        tree[c].set_parent_x(x[p]);
        tree[c].set_parent_y(y[p]);
        if (p != -1) {
            tree[p].add_child(c);
            tree[p].add_child_x(x[c]);
            tree[p].add_child_y(y[c]);
        }
    }

    return;
}


void calculateOverlapBottomUp(vector<Node> &tree, int N)
{
    stack<int> node_stack;
    int nodesProcessed = 0;
    int current_node = 0;
    int next_node = 0;
    node_stack.push(current_node);
    tree[current_node].print_children();
      
    while(nodesProcessed < N) {
        current_node = node_stack.top();
        next_node = tree[current_node].pick_next_child();
        if(next_node == -1) {
            node_stack.pop(); 
            nodesProcessed++;
            tree[current_node].partialLRST();
            cout << "popping: " << current_node << "...Nodes processed: " << nodesProcessed << endl;
        } else {
            cout << "next node is: " << next_node << endl;
            node_stack.push(next_node);
        }
    }
}

void runLRST(int *graph_D, int *graph_y, int *graph_x, int x[], int y[], int parent_nodes[], int child_nodes[], int N)
{
    // initialize tree
    vector<Node> tree;
    createTree(tree, parent_nodes, child_nodes, x, y, N); 

    calculateOverlapBottomUp(tree, N);
    return;
}



int steinerCalculation(int *graph_D, int *graph_y, int *graph_x, int *x_pts, int *y_pts, int bottom, int middle, int top)  {

    // function needs to return the proper orientation of connections between nodes
    // 1 = lower for bottom connection, lower for upper connection
    // 2 = lower for bottom connection, upper for upper connection
    // 3 = upper for bottom connection, lower for upper connection
    // 4 = upper for bottom connection, upper for upper connection

    int V1X, V1Y, V2X, V2Y;
    bool llOverlap = false, luOverlap = false, ulOverlap = false, uuOverlap = false;
    int value;

    //------------------------------------------------------------------------------------------
    // LOWER, LOWER

    if (*(x_pts + bottom) < *(x_pts + middle) && *(y_pts + bottom) > *(y_pts + middle)) {
        V1X = *(x_pts + bottom);
        V1Y = *(y_pts + middle);
    } else if (*(x_pts + bottom) < *(x_pts + middle) && *(y_pts + bottom) < *(y_pts + middle)) {
        V1X = *(x_pts + middle);
        V1Y = *(y_pts + bottom);
    } else if (*(x_pts + bottom) > *(x_pts + middle) && *(y_pts + bottom) > *(y_pts + middle)) {
        V1X = *(x_pts + bottom);
        V1Y = *(y_pts + middle);
    } else if (*(x_pts + bottom) > *(x_pts + middle) && *(y_pts + bottom) < *(y_pts + middle)) {
        V1X = *(x_pts + middle);
        V1Y = *(y_pts + bottom);
    }

    if (*(x_pts + middle) < *(x_pts + top) && *(y_pts + middle) > *(y_pts + top)) {
        V2X = *(x_pts + middle);
        V2Y = *(y_pts + top);
    } else if (*(x_pts + middle) < *(x_pts + top) && *(y_pts + middle) < *(y_pts + top)) {
        V2X = *(x_pts + top);
        V2Y = *(y_pts + middle);
    } else if (*(x_pts + middle) > *(x_pts + top) && *(y_pts + middle) > *(y_pts + top)) {
        V2X = *(x_pts + middle);
        V2Y = *(y_pts + top);
    } else if (*(x_pts + middle) > *(x_pts + top) && *(y_pts + middle) < *(y_pts + top)) {
        V2X = *(x_pts + top);
        V2Y = *(y_pts + middle);
    }

    cout << V1X << V1Y << V2X << V2Y << *(x_pts + middle) << *(y_pts + middle) << endl;

    if (V1X == V2X) {
        llOverlap = true;
    }

    //------------------------------------------------------------------------------------------

    //------------------------------------------------------------------------------------------
    // LOWER, UPPER

    if (*(x_pts + bottom) < *(x_pts + middle) && *(y_pts + bottom) > *(y_pts + middle)) {
        V1X = *(x_pts + bottom);
        V1Y = *(y_pts + middle);
    } else if (*(x_pts + bottom) < *(x_pts + middle) && *(y_pts + bottom) < *(y_pts + middle)) {
        V1X = *(x_pts + middle);
        V1Y = *(y_pts + bottom);
    } else if (*(x_pts + bottom) > *(x_pts + middle) && *(y_pts + bottom) > *(y_pts + middle)) {
        V1X = *(x_pts + bottom);
        V1Y = *(y_pts + middle);
    } else if (*(x_pts + bottom) > *(x_pts + middle) && *(y_pts + bottom) < *(y_pts + middle)) {
        V1X = *(x_pts + middle);
        V1Y = *(y_pts + bottom);
    }

    if (*(x_pts + middle) < *(x_pts + top) && *(y_pts + middle) > *(y_pts + top)) {
        V2X = *(x_pts + top);
        V2Y = *(y_pts + middle);
    } else if (*(x_pts + middle) < *(x_pts + top) && *(y_pts + middle) < *(y_pts + top)) {
        V2X = *(x_pts + middle);
        V2Y = *(y_pts + top);
    } else if (*(x_pts + middle) > *(x_pts + top) && *(y_pts + middle) > *(y_pts + top)) {
        V2X = *(x_pts + top);
        V2Y = *(y_pts + middle);
    } else if (*(x_pts + middle) > *(x_pts + top) && *(y_pts + middle) < *(y_pts + top)) {
        V2X = *(x_pts + middle);
        V2Y = *(y_pts + top);
    }

    cout << V1X << V1Y << V2X << V2Y << *(x_pts + middle) << *(y_pts + middle) << endl;

    if (V1Y == V2Y) {
        luOverlap = true;
    }

    //------------------------------------------------------------------------------------------

    //------------------------------------------------------------------------------------------
    // UPPER, LOWER

    if (*(x_pts + bottom) < *(x_pts + middle) && *(y_pts + bottom) > *(y_pts + middle)) {
        V1X = *(x_pts + middle);
        V1Y = *(y_pts + bottom);
    } else if (*(x_pts + bottom) < *(x_pts + middle) && *(y_pts + bottom) < *(y_pts + middle)) {
        V1X = *(x_pts + bottom);
        V1Y = *(y_pts + middle);
    } else if (*(x_pts + bottom) > *(x_pts + middle) && *(y_pts + bottom) > *(y_pts + middle)) {
        V1X = *(x_pts + middle);
        V1Y = *(y_pts + bottom);
    } else if (*(x_pts + bottom) > *(x_pts + middle) && *(y_pts + bottom) < *(y_pts + middle)) {
        V1X = *(x_pts + bottom);
        V1Y = *(y_pts + middle);
    }

    if (*(x_pts + middle) < *(x_pts + top) && *(y_pts + middle) > *(y_pts + top)) {
        V2X = *(x_pts + middle);
        V2Y = *(y_pts + top);
    } else if (*(x_pts + middle) < *(x_pts + top) && *(y_pts + middle) < *(y_pts + top)) {
        V2X = *(x_pts + top);
        V2Y = *(y_pts + middle);
    } else if (*(x_pts + middle) > *(x_pts + top) && *(y_pts + middle) > *(y_pts + top)) {
        V2X = *(x_pts + middle);
        V2Y = *(y_pts + top);
    } else if (*(x_pts + middle) > *(x_pts + top) && *(y_pts + middle) < *(y_pts + top)) {
        V2X = *(x_pts + top);
        V2Y = *(y_pts + middle);
    }

    cout << V1X << V1Y << V2X << V2Y << *(x_pts + middle) << *(y_pts + middle) << endl;

    if (V1Y == V2Y) {
        ulOverlap = true;
    }

    //------------------------------------------------------------------------------------------

    //------------------------------------------------------------------------------------------
    // UPPER, UPPER

    if (*(x_pts + bottom) < *(x_pts + middle) && *(y_pts + bottom) > *(y_pts + middle)) {
        V1X = *(x_pts + middle);
        V1Y = *(y_pts + bottom);
    } else if (*(x_pts + bottom) < *(x_pts + middle) && *(y_pts + bottom) < *(y_pts + middle)) {
        V1X = *(x_pts + bottom);
        V1Y = *(y_pts + middle);
    } else if (*(x_pts + bottom) > *(x_pts + middle) && *(y_pts + bottom) > *(y_pts + middle)) {
        V1X = *(x_pts + middle);
        V1Y = *(y_pts + bottom);
    } else if (*(x_pts + bottom) > *(x_pts + middle) && *(y_pts + bottom) < *(y_pts + middle)) {
        V1X = *(x_pts + bottom);
        V1Y = *(y_pts + middle);
    }

    if (*(x_pts + middle) < *(x_pts + top) && *(y_pts + middle) > *(y_pts + top)) {
        V2X = *(x_pts + top);
        V2Y = *(y_pts + middle);
    } else if (*(x_pts + middle) < *(x_pts + top) && *(y_pts + middle) < *(y_pts + top)) {
        V2X = *(x_pts + middle);
        V2Y = *(y_pts + top);
    } else if (*(x_pts + middle) > *(x_pts + top) && *(y_pts + middle) > *(y_pts + top)) {
        V2X = *(x_pts + top);
        V2Y = *(y_pts + middle);
    } else if (*(x_pts + middle) > *(x_pts + top) && *(y_pts + middle) < *(y_pts + top)) {
        V2X = *(x_pts + middle);
        V2Y = *(y_pts + top);
    }

    cout << V1X << V1Y << V2X << V2Y << *(x_pts + middle) << *(y_pts + middle) << endl;

    if (V1X == V2X) {
        cout << "yo" << endl;
        uuOverlap = true;
    }
    
    if (llOverlap) {
        value = 1;
        totalOverlap++;
    } else if (luOverlap) {
        value = 2;
        totalOverlap++;
    } else if (ulOverlap) {
        value = 3;
        totalOverlap++;
    } else if (uuOverlap) {
        value = 4;
        totalOverlap++;
    } else {
        value = 5;
    }

    cout << "VAAAALLLLLLUUUUUEEEEEEEEEE: " << value << endl;
    cout << llOverlap << luOverlap << ulOverlap << uuOverlap << endl;

    return value;
}

void runSteiner(int *graph_D, int *graph_y, int *graph_x, int *x_pts, int *y_pts, int parent_nodes[], int child_nodes[], int N) {
    // First, we need to identify the nodes with more than one child and their children
    int orientation = 0;

    bool parents[N]; // will keep parents
    bool seen[N];

    for (int l = 0; l < N; l++) {
        parents[l] = false;
        seen[N] = false;
    }

    // parent_nodes has a list of all potential parents

    int potentialParent;
    int parentCounter = 0; // if this is greater than 1, then that node has a parent with atleast 2 children
    int counter = 0;
    int frequency[N];
    for (int i = 0; i < N; i++) {

        //first, get a potential parent
        potentialParent = i;

        //cout << "potential: " << potentialParent << endl;

        // now, we need to see how any times this parent has appeared

        for (int j = 0; j < N; j++) {
            if (potentialParent == parent_nodes[j]) {
                parentCounter++;
            }
        }

        //cout << "counter: " << parentCounter << endl;
        frequency[i] = parentCounter;

        if (parentCounter > 1) {
            parents[i] =  true;
        }
        parentCounter = 0;
    }

    for (int y = 0; y < N; y++) {
        cout << "counter: " << frequency[y] << endl;
    }

    int sortedChildren[N];
    int sortedParents[N];

    for (int i = 0; i < N; i++) {
        // i is the index of sorted children
        for (int j = 0; j < N; j++) {
            if (i == child_nodes[j]) {
                sortedParents[i] = parent_nodes[j];

            }
        }
    }

    int edgesDone = 0;
    int middleNode;
    int topNode;
    //int direction[N - 1];
    counter = 0;

    for (int i = 0; i < N; i++) {
        sortedChildren[i] = 0;
    }

    int maxChildren = 0;

    for (int i = 0; i < N; i++) {
        if (frequency[i] > maxChildren) {
            maxChildren = frequency[i];
        }
    }

    // begin bottom up traversal

    // we need a starting point

    int startingNode;
    int parentOfStartingNode;

    int numberConnected[N];
    bool parentFound = false;
    int nextTop;
    bool rootConnected = false;
    int direction[N];
    int sortedDirections[N];
    int potentialStart;
    int potentialMiddle;
    bool failedFirstTest;
    bool failedSecondTest;

    for (int i = 0; i < N; i++) {
        numberConnected[i] = 0;
        direction[i] = 0;
        sortedDirections[i] = 0;
    }
    while (!rootConnected) {
        for (int i = 0; i < N; i++) {
            if (numberConnected[0] == frequency[0]) {
                // root node calculations
                rootConnected = true;
                cout << "ROOT" << endl;
                break;
            }

            parentFound = false;
            

            failedFirstTest = false;
            failedSecondTest = false;

            cout << "This is a potential starting point: " << potentialStart << endl;
            cout << "This is a potential middle point: " << potentialMiddle << endl;


            if (frequency[potentialStart] != numberConnected[potentialStart]) { // parent's children arent all connected
                //cout << "Potential Start: " << potentialStart << endl;
                //cout << "Potential Middle: " << potentialMiddle << endl;
                failedFirstTest = true;
                cout << "Failed 1st test" << endl;
            } else {
                cout << "Passed 1st test" << endl;
            }


            // need to see if there is or isnt a connection above the starting node
            // if there is, then the test fails
            // if there isn't, then you pass

            if (sortedDirections[potentialStart] == 0) {
                cout << "Passed 2nd test" << endl;
            } else {
                failedSecondTest = true;
                cout << "Failed 2nd test" << endl;
            }

            if (failedFirstTest || failedSecondTest) {
                // the potential starting point is not valid
                // look for another one
                continue;
            } else {
                // we now a real starting point and end point
                startingNode = potentialStart;
                // starting node = all the children are connected and no connection right above
                middleNode = potentialMiddle;
                cout << "This is an actual starting point: " << startingNode << endl;
                cout << "This is an actual middle point: " << middleNode << endl;
            }

            //------------------------------------------------------------

            // we now have a starting point
            // we have to go up this branch until we've reached a parent

            // Remember: BOTTOM UP TRAVERSAL

            // check to see if the middle node is a parent
            if (frequency[middleNode] > 1) {
                // we know middle node is a parent
                cout << "This is case h and i: " << startingNode << endl;
                cout << "This is case h and i: " << middleNode << endl;
                //function call
                //needs to return the direction

                numberConnected[middleNode]++;
                sortedDirections[startingNode] = 1;
            } else {
                // the middle node just has 1 child
                topNode = sortedParents[middleNode];

                if (frequency[topNode] > 1) {
                    // middle: 1 child
                    // top: > 1 child
                    // case d, e, f, and g
                    cout << "This is case d, e, f, and g: " << startingNode << endl;
                    cout << "This is case d, e, f, and g: " << middleNode << endl;
                    cout << "This is case d, e, f, and g: " << topNode << endl;
                    //function call
                    orientation = steinerCalculation(graph_D, graph_y, graph_x, x_pts, y_pts, startingNode, middleNode, topNode);
                    if (orientation == 1) {
                        sortedDirections[startingNode] = 1;
                        sortedDirections[middleNode] = 1;
                    } else if (orientation == 2) {
                        sortedDirections[startingNode] = 1;
                        sortedDirections[middleNode] = 2;
                    } else if (orientation == 3) {
                        sortedDirections[startingNode] = 2;
                        sortedDirections[middleNode] = 1;
                    } else if (orientation == 4) {
                        sortedDirections[startingNode] = 2;
                        sortedDirections[middleNode] = 2;
                    } else {
                        sortedDirections[startingNode] = 3;
                        sortedDirections[middleNode] = 3;
                    }
                    numberConnected[middleNode]++;
                    numberConnected[topNode]++;
                    //sortedDirections[startingNode] = 1;
                    //sortedDirections[middleNode] = 1;
                } else {
                    // middle: 1 child
                    // top: 1 child
                    cout << "This is case a, b, and c: " << startingNode << endl;
                    cout << "This is case a, b, and c: " << middleNode << endl;
                    cout << "This is case a, b, and c: " << topNode << endl;
                    //function call

                    if (frequency[startingNode] == 0) {
                        orientation = steinerCalculation(graph_D, graph_y, graph_x, x_pts, y_pts, startingNode, middleNode, topNode);
                    } else {
                        // do this calculation (not in the slides)
                        //orientation = steiner2(graph_D, graph_y, graph_x, x_pts, y_pts, startingNode, middleNode, topNode, parent_nodes, child_nodes, sortedDirections, N);
                    }

                    orientation = steinerCalculation(graph_D, graph_y, graph_x, x_pts, y_pts, startingNode, middleNode, topNode);

                    if (orientation == 1) {
                        sortedDirections[startingNode] = 1;
                        sortedDirections[middleNode] = 1;
                    } else if (orientation == 2) {
                        sortedDirections[startingNode] = 1;
                        sortedDirections[middleNode] = 2;
                    } else if (orientation == 3) {
                        sortedDirections[startingNode] = 2;
                        sortedDirections[middleNode] = 1;
                    } else if (orientation == 4) {
                        sortedDirections[startingNode] = 2;
                        sortedDirections[middleNode] = 2;
                    } else {
                        sortedDirections[startingNode] = 3;
                        sortedDirections[middleNode] = 3;
                    }

                    numberConnected[middleNode] = 1;
                    numberConnected[topNode] = 1;
                    //sortedDirections[startingNode] = 1;
                    //sortedDirections[middleNode] = 1;
                }
            }

            for (int k = 0; k < N; k++) {
                cout << "Number Connected: " << numberConnected[k] << endl;
            }

            for (int k = 0; k < N; k++) {
                cout << "Sorted Directions: " << sortedDirections[k] << endl;
            }



        }

    }

    // Need to put sortedDirections in correct place of directions

    int directions[N];

    for (int i = 0; i < N; i++) {
        directions[i] = sortedDirections[child_nodes[i]]; 
    }

    cout << "DIRECTIONS" << endl;

    for (int p = 0; p < N; p++) {
        cout << directions[p] << endl;
    }
        
}

int steiner2(int *graph_D, int *graph_y, int *graph_x, int *x_pts, int *y_pts, int bottom, int middle, int top, int parent_nodes[], int child_nodes[], int sortedDirections[], int N) {
    // with this function, the starting node has atleast 1 child
    // bottom: atleast 1 child
    // middle: 1 child
    // top: 1 child

    int bottomDirection = 0;
    for (int i = 0; i < N; i++) {
        if (child_nodes[i] == bottom) {
            bottomDirection = sortedDirections[child_nodes[i]];
        }
    }
    return 0;
}

int steiner3(int *graph_D, int *graph_y, int *graph_x, int *x_pts, int *y_pts, int bottom, int middle) {
    // need to find connection between bottom and middle

    // first, fix the connecton between bottom and middle to lower

    int connection  = 0;

    // get the vertex between the 2 points for an upper connection
    int VX = *(x_pts + bottom);
    int VY = *(y_pts + middle);

    // we now have to solve for Z(T df)
    return 0;


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
    runPrim(*graph_D, *graph_y, *graph_x, parent_nodes, child_nodes, N);

    // Run L-RST Algorithm
    runLRST(*graph_D, *graph_y, *graph_x, x, y, parent_nodes, child_nodes, N);
    //used for testing
    //steinerCalculation(*graph_D, *graph_y, *graph_x, x, y, 5, 6, 7);

    // Run Steiner Algorithm
    //runSteiner(*graph_D, *graph_y, *graph_x, x, y, parent_nodes, child_nodes, N);

    calculateWL(parent_nodes, child_nodes, *graph_D, *graph_y, *graph_x, x, y, N);

    cout << "Overlap Total: " << totalOverlap << endl;
    
    return 0;
}
