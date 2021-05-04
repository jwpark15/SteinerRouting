#ifndef PRIMPARAMS_H
#define PRIMPARAMS_H

#include <vector>

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
        std::vector<bool> connectedPoints;
};
#endif
