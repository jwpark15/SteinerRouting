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


def get_MST_lines(filename, x, y):
    parents = []
    children = []
    with open(filename, 'r') as f:
        line = f.readline()
        while line != '':
            parents.append(int(line.split()[0]))
            children.append(int(line.split()[1]))

            line = f.readline()

    plot_lines = [] 
    for i in range(1,len(parents)):
        i1 = parents[i]
        i2 = children[i]
        plot_lines.append([(x[i1], y[i1]), (x[i2], y[i2])])

    return plot_lines


def plot_MST(file_pts, file_mst, ax):
    x, y = get_points(file_pts) 
    mst_lines = get_MST_lines(file_mst, x, y)    

    # plot
    lc = LineCollection(mst_lines)
    ax.scatter(x,y)
    ax.add_collection(lc)
    
    return ax


def get_grid_size(file_pts):
    end_ind = file_pts.rfind('_')
    first_ind = file_pts[:end_ind].rfind('_')
    grid_size = int(file_pts[first_ind+1:end_ind])
    return grid_size

# TODO add plots for L-RST and Kahng/Robins
def generate_plots(file_pts, file_mst):
    # initialize plots
    fig, ax = plt.subplots(1,3)
    grid_size = get_grid_size(file_pts)
    plt.setp(ax, xlim=(0, grid_size), ylim=(0, grid_size))
    ax[0].set_title("MST")
    ax[1].set_title("Final L-RST")
    ax[2].set_title("Final Kahng/Robins")
    
    # MST
    ax[0] = plot_MST(file_pts, file_mst, ax[0])


    # L-RST

    # Kahng/Robins


    plt.show()
    return

if __name__ == "__main__":
    generate_plots("Points/points_10_5.pts", "primResults.txt")
