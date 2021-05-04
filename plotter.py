import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from matplotlib.collections import LineCollection
import random

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
    ax.scatter(x[0], y[0], color="red") # source node
    ax.scatter(x[1:],y[1:], color="tab:orange")
    lc = LineCollection(mst_lines)
    ax.add_collection(lc)
    ax.grid()
    
    return ax



def get_LRST_lines(filename, x, y):
    parents = []
    children = []
    orientations = []
    with open(filename, 'r') as f:
        line = f.readline()
        while line != '':
            parents.append(int(line.split()[0]))
            children.append(int(line.split()[1]))
            orientations.append(int(line.split()[2]))

            line = f.readline()

    plot_lines = [] 
    corner_pts = []
    for i in range(1,len(parents)):
        i1 = parents[i]
        i2 = children[i]
        if (orientations[i] == 0):
            # lower L
            min_y = min(y[i1], y[i2])
            plot_lines.append([(min(x[i1], x[i2]), min_y), (max(x[i1], x[i2]), min_y)])
            if (y[i1] > y[i2]):
                plot_lines.append([(x[i1], y[i2]), (x[i1], y[i1])])
                corner_pts.append([x[i1], min_y])
            else:
                plot_lines.append([(x[i2], y[i1]), (x[i2], y[i2])])
                corner_pts.append([x[i2], min_y])
        else:
            # upper L
            max_y = max(y[i1], y[i2])
            plot_lines.append([(min(x[i1], x[i2]), max_y), (max(x[i1], x[i2]), max_y)])
            if (y[i1] < y[i2]):
                plot_lines.append([(x[i1], y[i1]), (x[i1], y[i2])])
                corner_pts.append([x[i1], max_y])
            else:
                plot_lines.append([(x[i2], y[i2]), (x[i2], y[i1])])
                corner_pts.append([x[i2], max_y])

    return plot_lines, corner_pts

def get_LRST_Steiner(lrst_lines, corner_pts, nodes):
    steiner_x = []
    steiner_y = []
    for pt in corner_pts:
        x = pt[0]
        y = pt[1]
        if pt in nodes: 
            continue
        counter = 0
        for line in lrst_lines:
            x_range = range(line[0][0], line[1][0] + 1)
            y_range = range(line[0][1], line[1][1] + 1)
            if ((len(x_range)==1) and (len(y_range) == 1)):
                continue
            if (x in x_range) and (y in y_range):
                counter += 1
        if counter >= 3:
            steiner_x.append(x)
            steiner_y.append(y)
    
    return steiner_x, steiner_y


def plot_LRST(file_pts, file_lrst, ax):
    x, y = get_points(file_pts) 
    nodes = [[x[i], y[i]] for i in range(len(x))]
    lrst_lines, corner_pts = get_LRST_lines(file_lrst, x, y)    
    steiner_x, steiner_y = get_LRST_Steiner(lrst_lines, corner_pts, nodes)

    # plot
    lc = LineCollection(lrst_lines)
    ax.scatter(x[0], y[0], color="red") # source node
    ax.scatter(x[1:],y[1:], color="tab:orange")
    ax.scatter(steiner_x, steiner_y, color="purple")
    ax.add_collection(lc)
    ax.grid()
    
    return ax

# Creates the recti-linearlized lines and points for last iteration's plot
# Last iteration refers to the set of coordinates that occur after the last line that is "0, 0" in KRResults.txt
def get_KR_lines(krConnections, numberOfIterations, gridSize):
    kr_lines = [] # contains the lines
    currentIteration = 0
    with open(krConnections, 'r') as f:
        while (currentIteration < numberOfIterations): # keep looping until the last iteration is found
            line = f.readline()
            if (int(line.split()[0]) == 0 and int(line.split()[1]) == 0):
                currentIteration = currentIteration + 1 # found a new iteration
            if (currentIteration == numberOfIterations): # reached the last iteration
                line = f.readline() # read line to get non-zero numbers
                break
        while(line != ''): # Line number for last iteration obtained
            currentPY = int(int(line.split()[0]) / (gridSize + 1))
            currentPX = int(int(line.split()[0]) - (currentPY * gridSize + currentPY))
            currentCY = int(int(line.split()[1]) / (gridSize + 1))
            currentCX = int(int(line.split()[1]) - (currentCY * gridSize + currentCY))
            if (currentPX != currentCX and currentPY != currentCY): # diagonal
                kr_lines.append([(currentPX, currentPY), (currentCX, currentPY)])
                kr_lines.append([(currentCX, currentPY), (currentCX, currentCY)])
            else:
                kr_lines.append([(currentPX, currentPY), (currentCX, currentCY)]) # verticle or horizontal
            line = f.readline()
    return kr_lines

# Function calls get_KR_lines to get the points and connection
# Points and connections are then plotted
def plot_KR(file_pts, krConnections, file_kr_lengths, file_kr_steiners, gridSize, ax):
    x, y = get_points(file_pts) 
    numberOfIterations = 0
    with open(file_kr_lengths, 'r') as f:
        line = f.readline()
        while line != '':
            numberOfIterations += 1
            line = f.readline()
    
    kr_lines = get_KR_lines(krConnections, numberOfIterations, gridSize)

    # plot
    ax.scatter(x[0], y[0], color="red") # source node
    ax.scatter(x[1:],y[1:], color="tab:orange")
    lc = LineCollection(kr_lines)
    #ax.scatter(parentX, parentY)
    #ax.scatter(childX, childY)
    with open(file_kr_steiners, 'r') as f:
        line = f.readline()
        while line != '':
            steinerY = int(int(line) / (gridSize + 1))
            steinerX = int(int(line) - (steinerY*gridSize + steinerY))
            ax.scatter(steinerX, steinerY, c = 'Purple')
            line = f.readline()
    ax.add_collection(lc)
    ax.grid()
    
    return ax

# get grid size from file name
def get_grid_size(file_pts):
    end_ind = file_pts.rfind('_')
    first_ind = file_pts[:end_ind].rfind('_')
    grid_size = int(file_pts[first_ind+1:end_ind])
    return grid_size

def get_wirelengths(file_WLs):
    with open(file_WLs, 'r') as f:
        line1 = f.readline()
        line2 = f.readline()
        line3 = f.readline()

    mst_wl = line1.strip().split()[-1]
    lrst_wl = line2.strip().split()[-1]
    kr_wl = line3.strip().split()[-1]
   
    return mst_wl, lrst_wl, kr_wl

def generate_plots(file_pts, file_mst, file_lrst, file_kr, file_kr_lengths, file_kr_steiners, file_WLs):
    # initialize plots
    fig, ax = plt.subplots(1,3)
    grid_size = get_grid_size(file_pts)
    plt.setp(ax, xlim=(-1, grid_size+1), ylim=(-1, grid_size+1))
    
    # get wirelengths for each graph
    MST_WL, LRST_WL, KR_WL = get_wirelengths(file_WLs)

    # MST
    ax[0].set_title("Minimum Spanning Tree | WL = " + str(MST_WL))
    ax[0] = plot_MST(file_pts, file_mst, ax[0])
    
    # L-RST
    ax[1].set_title("Final L-RST | WL = " + str(LRST_WL))
    ax[1] = plot_LRST(file_pts, file_lrst, ax[1])

    # Kahng/Robins
    ax[2].set_title("Final Kahng/Robins | WL = " + str(KR_WL))
    ax[2] = plot_KR(file_pts, file_kr, file_kr_lengths, file_kr_steiners, grid_size, ax[2])

    # display
    fig = plt.gcf()
    fig.set_size_inches(16, 5)
    try: 
        plt.show()
    except: 
        print("cannot show plot in Matplotlib")

    benchmark = file_pts[:-4]
    plt.savefig(f'{benchmark}.png')
    return

if __name__ == "__main__":
    with open("benchmark.txt", 'r') as f:
        points_file = f.readline()

    print(points_file)
    generate_plots(points_file, "primResults.txt", "LRSTResults.txt", "KRResults.txt", "KRWirelengths.txt", "KRSteiners.txt", "wirelengthResults.txt")
