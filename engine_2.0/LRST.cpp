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

#include "HorizLine.h"
#include "VertLine.h"
#include "LRST.h"

using namespace std;
#define DEBUG 0
#define VERBOSE_DEBUG 0

// Node functions


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

// for debugging
void Node::print_children()
{
    cout << "printing children: " << endl;
    for(auto it = begin(children); it != end(children); ++it)
        cout << *it << endl;
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




// ==============================================================


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




