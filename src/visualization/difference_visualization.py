import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from pathlib import Path
import argparse
import os

def isBaseMatrix(filename):
    endings = ["Est.csv", "Nord.csv", "Ouest.csv", "Sud.csv"]
    for ending in endings:
        if str(filename)[-len(ending):] == ending:
            return True
    return False

def findBaseMatrices(matrix_list):
    base_matrices = []
    for matrix in matrix_list:
        if isBaseMatrix(matrix):
            base_matrices.append(matrix)
    return base_matrices

def plotShadingMask(base_matrix, matrix_files):
    # Load base matrix
    base_data = np.genfromtxt(base_matrix, delimiter=',')
    comparison_matrices = []

    # Find the comparison matrices
    for matrix in matrix_files:
        # compare the firs part of the name and check if it corresponds to the name of the base matrix appended with something
        if str(matrix.name)[:len(base_matrix.stem)] == base_matrix.stem:
            comparison_matrices.append(matrix)

    for matrix in comparison_matrices:
        # Load comparison matrix
        comparison_data = np.genfromtxt(matrix, delimiter=',')
        
        # Compute difference
        diff_data = base_data - comparison_data
        
        # Plotting code...
        fig = plt.figure()
        
        # azimuth-latitude coordinates
        altitude = np.linspace(0, 90, diff_data.shape[1]+1)
        azimuth = np.linspace(0, 2*np.pi, diff_data.shape[0]+1)

        # Coordinate grid for the pcolormesh
        r, th = np.meshgrid(altitude, azimuth)

        ax1 = plt.subplot(projection="polar")
        plt.grid(False)

        im = plt.pcolormesh(th, r, diff_data, cmap=cm.gray_r , vmin=0, vmax=1)

        # Setting colorbar ticks
        v1 = np.linspace(0, 1, 11)
        cbar = plt.colorbar(im,ticks=v1)
        cbar.ax.set_yticklabels(["{:4.2f}".format(i) for i in v1]) # add the labels


        # Grid on azimuth: 0 is in the North position, 
        # the axis is reverted with respect to usual polar coordinates
        # steps of 15° in the plot
        plt.thetagrids([i*15 for i in range(0,24)])
        ax1.set_theta_direction(-1)
        ax1.set_theta_zero_location("N")

        # Grid on altitude: 0 is in the outer circle, 
        # steps inwards of 10° in the plot
        plt.rgrids([i*10 for i in range(0,10)])
        ax1.set_rlim(bottom=90, top=0)

        # Save the plot...
        output_file = os.path.join(p.destination, 'Difference_Shading_mask_'+Path(matrix).stem+'.png')
        plt.savefig(output_file)
        plt.close(fig)

# Read the command line
parser = argparse.ArgumentParser()
parser.add_argument("--dir_path", type=Path, help='Insert the name of the directory containing the csv files with the shading mask matrices')
parser.add_argument("--destination", type=Path, help='Destination of the shading mask images. Default: location of the python script')

p = parser.parse_args()

if p.dir_path and p.dir_path.exists():
    matrix_files = list(p.dir_path.glob("*.csv"))
    base_matrices = findBaseMatrices(matrix_files)

    for base_matrix in base_matrices:
        plotShadingMask(base_matrix, matrix_files)
