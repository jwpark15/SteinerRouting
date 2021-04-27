import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import argparse
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import imageio
from PIL import Image 
import PIL 

images = []
imageCounter = 0
krWirelengths = []

# REFERENCE
# number = (y*size + y) + x
# y = number / (size + 1)
# x = number - (y*size+y)

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

def get_kr_pts(filename, gridSize):
    global imageCounter
    global images
    global krWirelengths
    parentX = []
    parentY = []
    childX = []
    childY = []
    currentPX = 0
    currentPY = 0
    currentCY = 0
    currentCX = 0
    parent = 0
    child = 0
    with open(filename, 'r') as f:
        line = f.readline()
        line = f.readline()
        while line != '':
            #print(line)
            # y = number / (size + 1)
            # x = number - (y*size+y)
            parent = int(line.split()[0])
            child = int(line.split()[1])
            #print(parent)
            #print(child)

            if parent == 0:
                # this is a new iteration
                # right now, we have all the x and y data we need
                # plot the data
                lines = get_KR_lines(parentX, parentY, childX, childY)

                fig = plt.figure()
                ax = plt.axes(xlim=(0, 31), ylim=(0, 31))

                lc = LineCollection(lines)
                ax.scatter(parentX, parentY)
                ax.scatter(childX, childY)
                plt.grid()
                plt.title("Wirelength: " + krWirelengths[imageCounter])
                plt.xticks(np.arange(0, 31))
                plt.yticks(np.arange(0, 31))
                ax.add_collection(lc)
                name = str(imageCounter) + ".png"
                print(imageCounter)
                imageCounter = imageCounter + 1
                images.append(plt.savefig(name))
                plt.savefig(name)

                line = f.readline()
                if (line != ''):
                    parentY.clear()
                    parentX.clear()
                    childY.clear()
                    childX.clear()
            else:
                currentPY = int(parent / (gridSize + 1))
                currentCY = int(child / (gridSize + 1))
                currentPX = int(parent - (currentPY*gridSize + currentPY))
                currentCX = int(child - (currentCY*gridSize + currentCY))
                parentY.append(currentPY)
                parentX.append(currentPX)
                childY.append(currentCY)
                childX.append(currentCX)
                line = f.readline()
            
    return parentX, parentY, childX, childY


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

def get_KR_lines(parentX, parentY, childX, childY):

    plot_lines = [] 
    for i in range(0,len(parentX)):
        # y = number / (size + 1)
        # x = number - (y*size+y)
        plot_lines.append([(parentX[i], parentY[i]), (childX[i], childY[i])])
    
    #print(plot_lines)

    return plot_lines



def plot_MST(file_pts, file_mst, ax):
    x, y = get_points(file_pts) 
    mst_lines = get_MST_lines(file_mst, x, y)    

    # plot
    lc = LineCollection(mst_lines)
    ax.scatter(x,y)
    ax.add_collection(lc)
    
    return ax

def plot_KR(file_kr, grid_size, ax):
    global images
    parentX, parentY, childX, childY = get_kr_pts(file_kr, grid_size)
    print(images[0])
    #imageio.mimsave('./test2.gif', images)
    kr_lines = get_KR_lines(parentX, parentY, childX, childY)
    #print(kr_lines)

    # plot
    lc = LineCollection(kr_lines)
    print(lc)
    ax.scatter(parentX, parentY)
    ax.scatter(childX, childY)
    ax.add_collection(lc)

    return ax


def get_grid_size(file_pts):
    end_ind = file_pts.rfind('_')
    first_ind = file_pts[:end_ind].rfind('_')
    grid_size = int(file_pts[first_ind+1:end_ind])
    return grid_size

# TODO add plots for L-RST and Kahng/Robins
def generate_plots(file_pts, file_mst, file_kr, file_kr_lengths):
    global krWirelengths
    with open(file_kr_lengths, 'r') as f:
        line = f.readline()
        while line != '':
            krWirelengths.append(line)
            line = f.readline()
    # initialize plots
    fig, ax = plt.subplots(1,2)
    grid_size = get_grid_size(file_pts)
    plt.setp(ax, xlim=(0, grid_size + 1), ylim=(0, grid_size + 1))
    plt.grid()
    ax[0].set_title("MST")
    ax[1].set_title("Final L-RST")
    plt.xticks(np.arange(0, grid_size+1))
    plt.yticks(np.arange(0, grid_size+1))
    #ax[2].set_title("Final Kahng/Robins")
    
    # MST
    ax[0] = plot_MST(file_pts, file_mst, ax[0])


    # L-RST

    # Kahng/Robins
    ax[1] = plot_KR(file_kr, grid_size, ax[1])

    plt.show()
    return



if __name__ == "__main__":
    generate_plots("Points/points_30_100.pts", "primResults.txt", "KRResults.txt", "KRWirelengths.txt")
