"""!
@file compute_granudrum_aor.py
@brief This file contains functions to compute the angle of repose in a similar way to the GranuDrum from GranuTools.

It includes functions to compute the angle of repose, write the result to a file, and a main function to parse command
line arguments.
"""
import argparse
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.legend_handler import HandlerPatch
from matplotlib.patches import Circle

from Tools.read_mercury_cg import mercury_cg


def compute_granudrum_aor(filename):
    """!
    @brief This function computes the angle of repose in a similar way to the GranuDrum from GranuTools.

    @param[in] filename The name of the coarse-grained .stat file from MercuryCG.
    """

    # load the data from a coarse grained stat file
    data = mercury_cg(filename + ".stat")

    # Create the figure and axes
    fig, ax = plt.subplots(figsize=(10, 8))
    dim = len(data["x"].shape)
    # this line of code is needed for this code to work for timesmoothed and timeaveraged CG data
    if dim == 2:
        # If the shape of the data is 2D, i.e. no time dimension we add an empty singleton dimension
        for key in data:
            data[key] = np.expand_dims(data[key], axis=-1)
    # Plot the contour plot
    contour = ax.contourf(data["x"][:, :, -1], data["z"][:, :, -1], data["VolumeFraction"][:, :, -1], cmap='viridis',
                          levels=20)
    # Add a colorbar with label and formatting
    fig.colorbar(contour, label='VolumeFraction', format='%1.2f')
    # Add labels and title
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Z-axis')
    ax.set_title('Contour Plot of volume fraction')
    # Add a grid for better readability
    ax.grid(True, linestyle='--', alpha=0.7)

    # Plot the 0.1 volumen fraction contour line
    contour_line_01 = ax.contour(data["x"][:, :, -1], data["z"][:, :, -1], data["VolumeFraction"][:, :, -1],
                                 colors="red", levels=[np.max(data.VolumeFraction[:, :, -1])*0.5])
    ax.clabel(contour_line_01, inline=True, fmt='%1.1f', fontsize=16)

    # Calculate the center and radius of the circle
    center_x = (np.max(data["x"]) + np.min(data["x"])) / 2
    center_z = (np.max(data["z"]) + np.min(data["z"])) / 2
    radius = (np.max(data["x"]) - np.min(data["x"])) / 10  # D_drum/5

    # Create a circle of d = D_drum/5
    circle = Circle((center_x, center_z), radius, fill=False, label='$D_{Drum}/5$ circle', color='b')

    # Find the intersection points with a tolerance of 1% of the circle radius
    intersection_points = []
    for line in contour_line_01.collections[0].get_paths():
        for vert in line.vertices:
            x, z = vert
            if np.isclose((x - center_x) ** 2 + (z - center_z) ** 2, radius ** 2, atol=radius * 1e-2):
                intersection_points.append((x, z))
    # Save the x and z values of the intersection points
    intersection_x = [point[0] for point in intersection_points]
    intersection_z = [point[1] for point in intersection_points]
    # Draw the circle and the line connecting the intersection points
    ax.add_patch(circle)
    ax.plot(intersection_x, intersection_z, 'ro-', label='Intersection points')

    # Fit a 1st degree polynomial (a line) to the intersection points
    slope, intercept = np.polyfit(intersection_x, intersection_z, 1)
    # Generate x values
    line_x = np.linspace(min(intersection_x), max(intersection_x), 100)
    # Calculate the corresponding y values
    line_z = slope * line_x + intercept
    # Plot the line
    ax.plot(line_x, line_z, 'g-', label='angle of repose line', linewidth=2)

    # print the function of the line
    print(f"f(x) = {slope} * x + {intercept}", )
    # calculate the absolute angle of repose
    angle_drum = np.abs(np.arctan(slope)) * 180 / np.pi
    print(f"absolute angle of repose = {angle_drum} degrees")

    # matplotlib inline function which creates a legend entry in form a circle for our d/5 circle
    class HandlerCircle(HandlerPatch):
        def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
            center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
            p = mpatches.Circle(xy=center, radius=height / 2)
            self.update_prop(p, orig_handle, legend)
            p.set_transform(trans)
            return [p]

    # add legend and specify the custom handler for the circle
    ax.legend(handler_map={mpatches.Circle: HandlerCircle()})
    # show the plot. set block to True to see it
    plt.show(block=False)
    # save the created figure
    plt.savefig(f'{filename}_AOR.png', dpi=300)
    # write the angle of repose to a file
    write_to_file(f'{filename}.txt', f"{angle_drum:.2f}")


# function to write a string to a file
def write_to_file(filename, filecontent):
    """!
    @brief This function writes a string to a file.

    @param[in] filename The name of the file to write to.
    @param[in] filecontent The content to write to the file.

    @return True if the file was successfully written, False otherwise.
    """

    try:
        with open(filename, 'w') as file:
            file.write(filecontent)
        return True
    except OSError:
        print('\033[91m' + 'Error in writeToFile: file could not be opened' + '\033[0m')
        return False


"""!
@brief This is the main function. It parses command line arguments and calls the compute_granudrum_aor function.
"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute the angle of repose in a similar way to the GranuDrum from '
                                                 'GranuTools', usage='%(prog)s <filename>')
    parser.add_argument('filename', type=str, help='filename which should be the name of a coarse-grained '
                                                   '.stat file from MercuryCG')
    args = parser.parse_args()

    compute_granudrum_aor(args.filename)
