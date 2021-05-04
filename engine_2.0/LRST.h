#ifndef LRST_H
#define LRST_H

#include <vector>

#include "VertLine.h"
#include "HorizLine.h"

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
        void partialLRST(std::vector<Node> &tree);
        void update_upper_Z(int overlap, int orientation);
        void update_lower_Z(int overlap, int orientation);
        void calculate_upperL_overlap(std::vector<Node> &tree);
        void calculate_lowerL_overlap(std::vector<Node> &tree);
        int processRoot(std::vector<Node> &tree);

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


void add_LowerL_to_Line_vectors(std::vector<VertLine> &vert_lines, std::vector<HorizLine> &horiz_lines, int x1, int x2, int y1, int y2);
void add_UpperL_to_Line_vectors(std::vector<VertLine> &vert_lines, std::vector<HorizLine> &horiz_lines, int x1, int x2, int y1, int y2);
int calculate_vertical_overlap(std::vector<VertLine> &vert_lines);
int calculate_horizontal_overlap(std::vector<HorizLine> &horiz_lines);
int calculate_Z1(std::vector<VertLine> &vert_lines, std::vector<HorizLine> &horiz_lines);
void createTree(std::vector<Node> &tree, int parent_nodes[], int child_nodes[], int x[], int y[], int N);
int calculateOverlapBottomUp(std::vector<Node> &tree, int N);
int setChildOrientation(std::vector<Node> &tree, int parent_node, int child_index, int isUpper);
void assign_Ls_top_down(std::vector<Node> &tree, int L_assignments[], int N);
void writeLRSTResults(int L[], int parent_nodes[], int child_nodes[], int N);
int runLRST(int *graph_D, int *graph_y, int *graph_x, int x[], int y[], int parent_nodes[], int child_nodes[], int N);


#endif
