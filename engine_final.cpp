#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <cstdio>
#include <vector>
#include <stack>
#include <math.h>
#include <algorithm> //for min/max
#include <ctime>


using namespace std;

// set to 0 to disable debugging prints
#define DEBUG 0
#define VERBOSE_DEBUG 0


// Global variables for KR
// ---------------------------------------------------------
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
// ---------------------------------------------------------


// update graph 2D arrays with weights of Manhattan Distance, y difference, x maximum
// row and column dimensions should be the same for all arrays
void createWeightGraphs(int *graph_D, int *graph_y, int *graph_x, int *x_pts, int *y_pts, int N)
{
    int x1, x2, y1, y2; // x and y coordinates of points to compare
    int ptr_val = 0;
    for (int i=0; i < N; ++i) {
        x1 = *(x_pts + i);
        y1 = *(y_pts + i);
        for (int j=0; j < N; ++j) {
            ptr_val = i*N + j;  // used to access next array element
            x2 = *(x_pts + j);
            y2 = *(y_pts + j);
            *(graph_D + ptr_val) = abs(x1 - x2) + abs(y1 - y2); // Manhattan distance
            *(graph_y + ptr_val) = -abs(y1 - y2); // difference in y points
            // x max
            if (x1 > x2) {
                *(graph_x + ptr_val) = -1 * x1;
            } else {
                *(graph_x + ptr_val) = -1 * x2;
            }
            
            // Print for debug 
            if (VERBOSE_DEBUG) {
                cout << "i: " << i << ", j: " << j << endl;
                cout << "<" << *(graph_D + ptr_val) << ", " << *(graph_y + ptr_val) << ", " << *(graph_x + ptr_val) << ">" << endl;
            }
        } 
    }
    return;
}
            

// write results of Prim's Algorithm to file
// each line of format: parent child 
// -1 is parent of root
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

// class to store important parameters for running Prim
class PrimParams {
    public:
        PrimParams(int N);
        // basic getters/setters
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

// constructor for initialization
PrimParams::PrimParams(int N)
{
    for(int i = 0; i<N; ++i) {
        connectedPoints.push_back(false);
    }

    connectedPoints[0] = true;
}

// reset parameters for new run
void PrimParams::refresh_all()
{
    dist = 999999;
    y_diff = 999999;
    x_max = 999999;
    parent = 0;
    child = 0;
}


// set all major parameters at once
void PrimParams::set_all(int d, int y, int x, int p, int c) 
{
    dist = d;
    y_diff = y;
    x_max = x;
    parent = p;
    child = c;
}

// check weights and update values if necessary
// used to determine the proper parent/child relationships
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

// add node to MST with correct parent/child relationship
void PrimParams::add_node(int *graph_D, int *graph_y, int *graph_x, int N)
{
    int d, y, x, reference, endPoint; // temp parameters

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

            check_params(d, y, x, reference, endPoint); // check against current best weights and update if necessary
        }
    }
}

// generate Minimum Spanning Tree
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

        parent_nodes[numberConnected] = params.get_parent();
        child_nodes[numberConnected] = params.get_child();
        numberConnected++;

        params.refresh_all(); // reset params for next node
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

// class to hold information for each tree node
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
        int get_lower_Z(void) {return lower_Z;}
        int get_upper_Z(void) {return upper_Z;}
        int get_upper(void) {return upper;}
        int get_lower(void) {return lower;}
        int get_child_pointer(void) {return child_pointer;}
        
        // funcs for tree traversal 
        int pick_next_child(void);
        void reset_child_pointer(void) {child_pointer = 0;}
        void partialLRST(vector<Node> &tree);
        void update_upper_Z(int overlap, int orientation);
        void update_lower_Z(int overlap, int orientation);
        void calculate_upperL_overlap(vector<Node> &tree);
        void calculate_lowerL_overlap(vector<Node> &tree);
        int processRoot(vector<Node> &tree);

        //debugging
        void print_children();
    private:
        int tag;
        int x;
        int y;
        int lower_Z;
        int upper_Z;
        int upper; // binary representation of orientations for children if node is upper. 0 = lower, 1 = upper. LSB first index of children
        int lower; // binary representation of orientations for children if node is lower. 0 = lower, 1 = upper. LSB first index of children
        std::vector<int> children; // tags of all children
        std::vector<int> children_x;
        std::vector<int> children_y;
        int parent; // tag of parent
        int parent_x;
        int parent_y;
        int child_pointer; // keep track of which branches have been traversed
};

// for debugging
void Node::print_children()
{
    cout << "printing children: " << endl;
    for(auto it = begin(children); it != end(children); ++it)
        cout << *it << endl;
}
  
// constructor for parameter initialization
Node::Node(int val, int x_coord, int y_coord)
{
    tag = val;
    upper = 0;
    lower = 0;
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


// get the best gain for upper L
void Node::update_upper_Z(int overlap, int orientation)
{
    if (overlap > upper_Z) {
        upper_Z = overlap;
        upper = orientation;
    }
}

// get the best gain for lower L
void Node::update_lower_Z(int overlap, int orientation)
{
    if (overlap > lower_Z) {
        lower_Z = overlap;
        lower = orientation;
    }
}

struct HorizLine {
    int y;
    int x1;
    int x2;
};

struct VertLine {
    int x;
    int y1;
    int y2;
};


// add break lower L into horizontal and vertical components
// add components to collection of line segments
void add_LowerL_to_Line_vectors(vector<VertLine> &vert_lines, vector<HorizLine> &horiz_lines, int x1, int x2, int y1, int y2)
{
    HorizLine Lbottom;
    VertLine Lside;
    Lbottom.x1 = min(x1, x2);
    Lbottom.x2 = max(x1, x2);
    Lbottom.y = (y1 < y2) ? y1 : y2;

    Lside.y1 = min(y1, y2);
    Lside.y2 = max(y1, y2);
    Lside.x = (y1 < y2) ? x2 : x1;
    
    vert_lines.push_back(Lside);
    horiz_lines.push_back(Lbottom);
   

    // cout << "adding Lower L: bottom - {" << Lbottom.y << ", " << Lbottom.x1 << ", " << Lbottom.x2 << "}; side - {" << Lside.x << ", " << Lside.y1 << ", " << Lside.y2 << "}" << endl;
}


// add break upper L into horizontal and vertical components
// add components to collection of line segments
void add_UpperL_to_Line_vectors(vector<VertLine> &vert_lines, vector<HorizLine> &horiz_lines, int x1, int x2, int y1, int y2)
{
    HorizLine Ltop;
    VertLine Lside;
    Ltop.x1 = min(x1, x2);
    Ltop.x2 = max(x1, x2);
    Ltop.y = (y1 > y2) ? y1 : y2;

    Lside.y1 = min(y1, y2);
    Lside.y2 = max(y1, y2);
    Lside.x = (y1 > y2) ? x2 : x1;
    
    vert_lines.push_back(Lside);
    horiz_lines.push_back(Ltop);
   
    /*
    if (DEBUG) { 
        cout << "adding Upper L: top - {" << Ltop.y << ", " << Ltop.x1 << ", " << Ltop.x2 << "}; side - {" << Lside.x << ", " << Lside.y1 << ", " << Lside.y2 << "}" << endl;
    }
    */
}



// calculate amount of overlap given collection of vertical line segments
int calculate_vertical_overlap(vector<VertLine> &vert_lines)
{
    int overlap = 0;
    int size_c;
    vector<VertLine> combined_vert_lines;
    for (auto vert = vert_lines.begin(); vert != vert_lines.end(); ++vert) {
        if (combined_vert_lines.size() == 0) {
            combined_vert_lines.push_back((*vert));
            continue;
        }
        size_c = combined_vert_lines.size();
        for (int i = 0; i < size_c; ++i) {
            // overlap only if x coordinates are the same
            if ((*vert).x == combined_vert_lines[i].x) {
                overlap += max(0, min((*vert).y2, combined_vert_lines[i].y2) - max((*vert).y1, combined_vert_lines[i].y1));
                if ((*vert).y1 < combined_vert_lines[i].y1) {
                    combined_vert_lines[i].y1 = (*vert).y1;
                }
                if ((*vert).y2 > combined_vert_lines[i].y2) {
                    combined_vert_lines[i].y2 = (*vert).y2;
                }
            } else {
                combined_vert_lines.push_back((*vert));
            }
        }
    }
    if(VERBOSE_DEBUG) {cout << "vertical overlap: " << overlap << endl;}
    return overlap;
}


// calculate amount of overlap given collection of horizontal line segments
int calculate_horizontal_overlap(vector<HorizLine> &horiz_lines)
{
    int overlap = 0;
    int size_c;
    vector<HorizLine> combined_horiz_lines;
    for (auto horiz = horiz_lines.begin(); horiz != horiz_lines.end(); ++horiz) {
        if (combined_horiz_lines.size() == 0) {
            combined_horiz_lines.push_back((*horiz));
            continue;
        }
        size_c = combined_horiz_lines.size();
        for (int i = 0; i < size_c; ++i) {
            // overlap only if y coordinates are the same
            if ((*horiz).y == combined_horiz_lines[i].y) {
                overlap += max(0, min((*horiz).x2, combined_horiz_lines[i].x2) - max((*horiz).x1, combined_horiz_lines[i].x1));
                if ((*horiz).x1 < combined_horiz_lines[i].x1) {
                    combined_horiz_lines[i].x1 = (*horiz).x1;
                }
                if ((*horiz).x2 > combined_horiz_lines[i].x2) {
                    combined_horiz_lines[i].x2 = (*horiz).x2;
                }
            } else {
                combined_horiz_lines.push_back((*horiz));
            }
        }
    }
    if(VERBOSE_DEBUG) {cout << "horizontal overlap: " << overlap << endl;}
    return overlap;
}


// calculate Z1 overlap for collection of vertical and horizontal line segments
int calculate_Z1(vector<VertLine> &vert_lines, vector<HorizLine> &horiz_lines)
{
    int overlap = 0;
    overlap += calculate_vertical_overlap(vert_lines);
    overlap += calculate_horizontal_overlap(horiz_lines);
    
    //cout << "z1 Overlap: " << overlap << endl;
    return overlap;
}

// determine best orientation for producing most overlap for a Node
// assign Node/parent edge to upper L
void Node::calculate_upperL_overlap(vector<Node> &tree)
{
    int c = children.size();
    int mask, child, child_x, child_y, overlap, Z1, Z2; 
    // every possible orientation of children
    for (int i = 0; i < pow(2,c); ++i) {
        vector<VertLine> vert_lines;
        vector<HorizLine> horiz_lines;
        //cout << "parent (x,y) for tag " << tag << " {" << parent_x << ", " << parent_y << "}" << endl;
        add_UpperL_to_Line_vectors(vert_lines, horiz_lines, x, parent_x, y, parent_y);
        Z2 = 0;

        // add all lines
        for (int j = 0; j < c; ++j) {
            child = children[j];
            child_x = children_x[j];
            child_y = children_y[j];
            // fancy masking to get orientation of child
            mask = i & (1<<j);
            if(mask == 0) { // lower L
                if(VERBOSE_DEBUG) {cout << "Lower L " << i << ", " << j << endl;}
                Z2 += tree[child].get_lower_Z();
                add_LowerL_to_Line_vectors(vert_lines, horiz_lines, child_x, x, child_y, y);
            } else { // upper L
                if(VERBOSE_DEBUG) {cout << "Upper L " << i << ", " << j << endl;}
                Z2 += tree[child].get_upper_Z();
                add_UpperL_to_Line_vectors(vert_lines, horiz_lines, child_x, x, child_y, y);
            }
            Z1 = calculate_Z1(vert_lines, horiz_lines);
            overlap = Z1 + Z2;
            if(VERBOSE_DEBUG) {cout << "update upper: " << overlap << ", " << i << endl;}
            update_upper_Z(overlap, i);
        }
        if(VERBOSE_DEBUG) {cout << "UPPER iteration: " << i << ", best orientation: " << get_upper() << ", best Z: " << get_upper_Z() <<  endl << "--------------------------------" << endl;}
    }
}


// determine best orientation for producing most overlap for a Node
// assign Node/parent edge to lower L
void Node::calculate_lowerL_overlap(vector<Node> &tree)
{
    int c = children.size();
    int mask, child, child_x, child_y, overlap, Z1, Z2;
    // every possible orientation of children
    for (int i = 0; i < pow(2,c); ++i) {
        vector<VertLine> vert_lines;
        vector<HorizLine> horiz_lines;
        //cout << "parent (x,y) for tag " << tag << " {" << parent_x << ", " << parent_y << "}" << endl;
        add_LowerL_to_Line_vectors(vert_lines, horiz_lines, x, parent_x, y, parent_y);
        Z2 = 0;

        // add all lines
        for (int j = 0; j < c; ++j) {
            child = children[j];
            child_x = children_x[j];
            child_y = children_y[j];
            // fancy masking to get all orientations
            mask = i & (1<<j);
            if(mask == 0) { // lower L
                if(VERBOSE_DEBUG) {cout << "Lower L " << i << ", " << j << endl;}
                Z2 += tree[child].get_lower_Z();
                add_LowerL_to_Line_vectors(vert_lines, horiz_lines, child_x, x, child_y, y);
            } else { // upper L
                if(VERBOSE_DEBUG) {cout << "Upper L " << i << ", " << j << endl;}
                Z2 += tree[child].get_upper_Z();
                add_UpperL_to_Line_vectors(vert_lines, horiz_lines, child_x, x, child_y, y);
            }
            Z1 = calculate_Z1(vert_lines, horiz_lines);
            overlap = Z1 + Z2;
            if(VERBOSE_DEBUG) {cout << "update lower: " << overlap << ", " << i << endl;}
            update_lower_Z(overlap, i);
        }
        //DEBUG
        if(VERBOSE_DEBUG) {cout << "LOWER iteration: " << i << ", best orientation: " << get_lower() << ", best Z: " << get_lower_Z() <<  endl << "--------------------------------" << endl;}
    }
}

// determine best orientation of children to produce most possible overlap in the tree
// return reduction in wirelength
int Node::processRoot(vector<Node> &tree)
{
    int c = children.size();
    int mask, child, child_x, child_y, overlap, Z1, Z2; 
    // for every possible orientation of children
    for (int i = 0; i < pow(2,c); ++i) {
        vector<VertLine> vert_lines;
        vector<HorizLine> horiz_lines;
        Z2 = 0;

        // add all lines
        for (int j = 0; j < c; ++j) {
            child = children[j];
            child_x = children_x[j];
            child_y = children_y[j];
            // fancy masking to get all orientations
            mask = i & (1<<j);
            if(mask == 0) { // lower L
                if(VERBOSE_DEBUG) {cout << "Lower L" << i << ", " << j << endl;}
                Z2 += tree[child].get_lower_Z();
                add_LowerL_to_Line_vectors(vert_lines, horiz_lines, child_x, x, child_y, y);
            } else { // upper L
                if(VERBOSE_DEBUG) {cout << "Upper L" << i << ", " << j << endl;}
                Z2 += tree[child].get_upper_Z();
                add_UpperL_to_Line_vectors(vert_lines, horiz_lines, child_x, x, child_y, y);
            }
            Z1 = calculate_Z1(vert_lines, horiz_lines);
            overlap = Z1 + Z2;
            update_upper_Z(overlap, i); //using upper Z even though there is no parent 
        }
        if(VERBOSE_DEBUG) {cout << "ROOT iteration: " << i << ", best orientation: " << get_upper() << ", best Z: " << get_upper_Z() << endl << "--------------------------------" << endl;}
    }
    if (DEBUG) {
        cout << "--------------------------------" << endl;
        cout << "ROOT: orientation = " << upper << ", Z = " << upper_Z << endl;
        cout << "--------------------------------" << endl;
    }
    return upper_Z;
}

// run partial L-RST for a node
void Node::partialLRST(vector<Node> &tree)
{
    // skip bottom nodes
    if (children.size() != 0) {
        // upper L
        calculate_upperL_overlap(tree);

        // lower L
        calculate_lowerL_overlap(tree);
        if (DEBUG) {
            cout << "--------------------------------" << endl;
            cout << "UPPER: orientation = " << upper << ", Z = " << upper_Z << endl;
            cout << "LOWER: orientation = " << lower << ", Z = " << lower_Z << endl;
            cout << "--------------------------------" << endl;
        }
    }

}

// generate a vector of nodes that represent the minimum spanning tree
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

// traverse the tree bottom up to determine best gain possible from overlap
// return total overlap
int calculateOverlapBottomUp(vector<Node> &tree, int N)
{
    stack<int> node_stack;
    int nodesProcessed = 0;
    int current_node = 0;
    int next_node = 0;
    node_stack.push(current_node);
      
    //process all except root node
    while(nodesProcessed < (N-1)) {
        current_node = node_stack.top();
        next_node = tree[current_node].pick_next_child();
        if(next_node == -1) {
            node_stack.pop(); 
            nodesProcessed++;
            if(DEBUG) {
                cout << "-------------------------------" << endl << "popping: " << current_node << "...Nodes processed: " << nodesProcessed << endl << "-------------------------------" << endl;
            }
            tree[current_node].partialLRST(tree);
        } else {
            //cout << "next node is: " << next_node << endl;
            node_stack.push(next_node);
        }
    }
    cout << "processing root node..." << endl;
    int wirelength_reduction = tree[0].processRoot(tree);
    if(VERBOSE_DEBUG) {
        for (int a = 0; a < N; ++a) {
            cout << "INDEX: " << a << " ----------upper Z: " << tree[a].get_upper_Z() << "...lower Z: " << tree[a].get_lower_Z() << "----------" << endl;
            cout << "upper orientation: " << tree[a].get_upper() << "...lower orientation: " << tree[a].get_lower() << endl;
        }
    }
    return wirelength_reduction;
}

// set upper or lower orientation for L formed by parent and child node
int setChildOrientation(vector<Node> &tree, int parent_node, int child_index, int isUpper)
{
    int orientation = (isUpper) ? tree[parent_node].get_upper() : tree[parent_node].get_lower();
    int masked = orientation & (1<<child_index);
    if(VERBOSE_DEBUG) {cout << "orientation: " << orientation << "...mask: " << masked << endl;}
    if (masked == 0) {
        return 0; // lower
    } else { 
        return 1; // upper
    }

}

// traverse tree top --> down to assign upper/lower to all edges in tree
void assign_Ls_top_down(vector<Node> &tree, int L_assignments[], int N)
{
    int numberAssigned = 1;
    int current_node = 0;
    int next_node = 0;
    int overall_Z = 0;
    int isUpper = 1; // 1 for upper, 0 for lower
    int child_index = 0;
    stack<int> top_down_stack;
    top_down_stack.push(current_node);

    for(int i = 0; i<N; ++i) 
        tree[i].reset_child_pointer();

    while (numberAssigned < N) {
        current_node = top_down_stack.top();    
        next_node = tree[current_node].pick_next_child();
        //cout << "next child: " << next_node << endl;
        if(next_node == -1) {
            top_down_stack.pop();
            //cout << "popping: " << current_node << endl;
            isUpper = L_assignments[top_down_stack.top()];
        } else {
            numberAssigned++;
            if(VERBOSE_DEBUG) {cout << "nodes assigned: " << numberAssigned << endl;}
            child_index = (tree[current_node].get_child_pointer()) - 1;
            //cout << "Child indx: " << child_index << endl;
            int temp = setChildOrientation(tree, current_node, child_index, isUpper);
            L_assignments[next_node] = temp;
            isUpper = temp;
            top_down_stack.push(next_node);
        }
    }
}

// write results of LRST to file
// line format: parent child orientation
// orientation is 1 for upper, 0 for lower
// root parent is -1
void writeLRSTResults(int L[], int parent_nodes[], int child_nodes[], int N)
{
    ofstream results;
    results.open("LRSTResults.txt");
    for(int i = 0; i < N; ++i)
        results << parent_nodes[i] << " " << child_nodes[i] << " " <<  L[child_nodes[i]] << endl;
    results.close();
    return;

}

// run L-RST Algorithm
// return overlap
int runLRST(int *graph_D, int *graph_y, int *graph_x, int x[], int y[], int parent_nodes[], int child_nodes[], int N)
{
    // initialize tree
    vector<Node> tree;
    if(DEBUG) {cout << "--- creating tree / nodes ---" << endl;}
    createTree(tree, parent_nodes, child_nodes, x, y, N); 

    if(DEBUG) {cout << "--- calculating overlap ---" << endl;}
    int overlap = calculateOverlapBottomUp(tree, N);

    int L_assignments[N]; // for each index, 0 is lower, 1 is upper
    L_assignments[0] = 1;
    if(DEBUG) {cout << endl << "--- assigning Ls top --> down ---" << endl << endl;}
    assign_Ls_top_down(tree, L_assignments, N);

    if(DEBUG) {
        for (int i = 0; i < N; ++i)
            cout << "i: " << i << ", L assigned: " << L_assignments[i] << endl;
    }
    writeLRSTResults(L_assignments, parent_nodes, child_nodes, N);
    return overlap;
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
    results.open("wirelenghtResults.txt");
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


int main(int argc, char** argv)
{
    const char *filename = argv[1];
    std::string temp_filename = argv[1];
   
    int N = getNumNodes(temp_filename);
    exportFilename(filename);

    cout << "NUM NODES: " << N << endl;

    // N is number of nodes in graph
    // x and y coordinates for nodes
    int x[N];
    int y[N];
    
    parseFile(x, y, N, filename);
    
    /*
    if (DEBUG) {
        for (int i = 0; i < N; ++i) {
            cout << x[i] << ", " << y[i] << endl;
        }
    } 
    */

    // init graphs for Manhattan Distance, y difference, and x maximum
    int graph_D[N][N];
    int graph_y[N][N];
    int graph_x[N][N];

    int parent_nodes[N];
    int child_nodes[N];

    
    clock_t c_start = clock();

    cout << endl << "============== Create Weight Graph =================" << endl << endl;
    // create 2D arrays with values for Prim
    createWeightGraphs(*graph_D, *graph_y, *graph_x, x, y, N); 

    // Run Prim's Algorithm
    cout << endl << "============== Generating MST =================" << endl << endl;
    runPrim(*graph_D, *graph_y, *graph_x, parent_nodes, child_nodes, N);

    clock_t c_endPrim = clock();
    // Run L-RST Algorithm
    cout << endl << "============== Running L-RST =================" << endl << endl;
    int overlap_LRST = runLRST(*graph_D, *graph_y, *graph_x, x, y, parent_nodes, child_nodes, N);

    clock_t c_endLRST = clock();

    // Run Kahng / Robins Algorithm 
    finalSteinerParents.push_back(0);
    finalSteinerChildren.push_back(0);
    int gridSize = getGridSize(temp_filename);
    int initialWL = calculateWL(parent_nodes, child_nodes, *graph_D, *graph_y, *graph_x, x, y, N);
    hananGrid(parent_nodes, child_nodes, *graph_D, *graph_y, *graph_x, x, y, N, gridSize, initialWL);
    writeKRResults();


    clock_t c_endKR = clock();

    // calculate Wirelength
    int LRST_WL = initialWL - overlap_LRST;
    int KR_WL = 0; // TODO - update this 

    writeWirelengths(initialWL, LRST_WL, KR_WL);

    clock_t c_end = clock();
    double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC; 
    double time_elapsed_MST_ms = 1000.0 * (c_endPrim - c_start) / CLOCKS_PER_SEC; 
    double time_elapsed_LRST_ms = 1000.0 * (c_endLRST - c_endPrim) / CLOCKS_PER_SEC; 
    double time_elapsed_KR_ms = 1000.0 * (c_endKR - c_endLRST) / CLOCKS_PER_SEC; 

    cout << "===================================" << endl;
    cout << "Original Wirelength: " << initialWL << endl;
    cout << "Wirelength after L-RST: " << LRST_WL << endl;
    cout << "Wirelength after Kahng / Robins: " << KR_WL << endl;
    cout << "CPU time used for generating MST: " << time_elapsed_MST_ms / 1000.0 << " s\n";
    cout << "CPU time used for running L-RST: " << time_elapsed_LRST_ms / 1000.0 << " s\n";
    cout << "CPU time used for running Kahng / Robins: " << time_elapsed_KR_ms / 1000.0 << " s\n";
    cout << "===================================" << endl;



    return 0;
}
