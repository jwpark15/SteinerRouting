import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from celluloid import Camera
from math import floor

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
        #plot_lines.append([(x[i1], y[i1]), (x[i2], y[i2])])
        plot_lines.append(([x[i1], x[i2]], [y[i1], y[i2]]))

    return plot_lines


def plot_MST(file_pts, file_mst, ax, camera):
    x, y = get_points(file_pts) 
    plot_lines = get_MST_lines(file_mst, x, y)    
    for i in range(len(plot_lines)):
        ax.scatter(x,y, color='orange')
        for line in plot_lines[:i]:
            plt.plot(line[0], line[1], color='blue')
        camera.snap()

    # fix bug where last line doesn't plot
    ax.scatter(x,y, color='orange')
    for line in plot_lines:
        plt.plot(line[0], line[1], color='blue')
    camera.snap()

    return plot_lines



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
    for i in range(1,len(parents)):
        i1 = parents[i]
        i2 = children[i]
        if (orientations[i] == 0):
            # lower L
            min_y = min(y[i1], y[i2])
            plot_lines.append(( [min(x[i1], x[i2]), max(x[i1], x[i2])] , [min_y, min_y]))
            if (y[i1] > y[i2]):
                plot_lines.append(([x[i1], x[i1]], [y[i2], y[i1]]))
            else:
                plot_lines.append(([x[i2], x[i2]], [y[i1], y[i2]]))
        else:
            # upper L
            max_y = max(y[i1], y[i2])
            plot_lines.append(( [min(x[i1], x[i2]), max(x[i1], x[i2])] , [max_y, max_y]))
            if (y[i1] < y[i2]):
                plot_lines.append(([x[i1], x[i1]], [y[i2], y[i1]]))
            else:
                plot_lines.append(([x[i2], x[i2]], [y[i1], y[i2]]))

    return plot_lines


def plot_LRST(file_pts, file_lrst, MST_lines, ax, camera):
    x, y = get_points(file_pts) 
    plot_lines = get_LRST_lines(file_lrst, x, y)    
    for i in range(len(plot_lines)):
        ax.scatter(x,y, color='orange')
        for line in plot_lines[:i]:
            plt.plot(line[0], line[1], color='red')
        for m_line in MST_lines[floor(i/2):]:
            plt.plot(m_line[0], m_line[1], color='blue')
        camera.snap()

    # fix bug where last line doesn't plot
    ax.scatter(x,y, color='orange')
    for line in plot_lines:
        plt.plot(line[0], line[1], color='red')
    camera.snap()

    return plot_lines
    

def get_grid_size(file_pts):
    end_ind = file_pts.rfind('_')
    first_ind = file_pts[:end_ind].rfind('_')
    grid_size = int(file_pts[first_ind+1:end_ind]) + 1
    return grid_size

# TODO add plots for L-RST and Kahng/Robins
def generate_plots(file_pts, file_mst, file_lrst, file_wirelengths):
    # initialize plots
    fig, ax = plt.subplots(1)
    camera = Camera(fig)
    grid_size = get_grid_size(file_pts)
    plt.setp(ax, xlim=(0, grid_size), ylim=(0, grid_size))
    ax.set_title("L-RST Animation")
    ax.grid()
    MST_lines = plot_MST(file_pts, file_mst, ax, camera)

    ax = plot_LRST(file_pts, file_lrst, MST_lines, ax, camera)
    camera.snap()

    animation = camera.animate()
    animation.save("test.gif", writer = 'imagemagick')

    plt.show()
    return

if __name__ == "__main__":
    with open("benchmark.txt", 'r') as f:
        points_file = f.readline()

    print(points_file)
    generate_plots(points_file, "primResults.txt", "LRSTResults.txt", "wirelengthResults.txt")
