import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import argparse

def get_points(filename):
    x_pts = []
    y_pts = []
    with open(filename, 'r') as f:
        line = f.readline()
        while line != '':
            x_pts.append(int(line.split()[0]))
            y_pts.append(int(line.split()[1]))
            
            line = f.readline()
            
    return x_pts, y_pts

def plot_MST(file_pts, file_mst):
    x, y = get_points(file_pts) 
    parents = []
    children = []
    with open(file_mst, 'r') as f:
        line = f.readline()
        while line != '':
            parents.append(int(line.split()[0]))
            children.append(int(line.split()[1]))

            line = f.readline()
    
    plot_lines = [] 
    colors = []
    for i in range(1,len(parents)):
        i1 = parents[i]
        i2 = children[i]
        plot_lines.append([(x[i1], y[i1]), (x[i2], y[i2])])
        colors.append((0,0,0,1))

    print(plot_lines)
    lc = LineCollection(plot_lines, colors=colors)
    fig, ax = plt.subplots()
    ax.add_collection(lc)
    ax.autoscale()
    
    plt.show()
    
    return



if __name__ == "__main__":
    plot_MST("Points/points_10_5.pts", "primResults.txt")
