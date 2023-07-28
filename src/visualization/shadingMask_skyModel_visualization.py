import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from pathlib import Path
import argparse
import os
import re
from scipy import interpolate

def interpolate_to_fit(small, large):
    x = np.arange(0, small.shape[1], 1)
    y = np.arange(0, small.shape[0], 1)
    f = interpolate.interp2d(x, y, small, kind='linear')

    xnew = np.linspace(0, small.shape[1], large.shape[1])
    ynew = np.linspace(0, small.shape[0], large.shape[0])
    return f(xnew, ynew)

def plotShadingMask(csv_filename, destination):

    shading_test_values = np.genfromtxt(csv_filename, delimiter=',')
    
    # Extract the key part from the shading mask name
    shading_name = Path(csv_filename).stem
    key_part = re.search('SM_Matrix_(.*)_.*', shading_name).group(1)

    sky_models_directory = Path(csv_filename).parent / ".." / "skyModels"

    for hour in range(5, 20):
        sky_model_filename = sky_models_directory / f"SkyModel_Matrix_{key_part}_{hour}H.csv"
        if sky_model_filename.exists():
            sky_model_values = np.genfromtxt(sky_model_filename, delimiter=',')
            if shading_test_values.shape != sky_model_values.shape:
                if np.prod(shading_test_values.shape) < np.prod(sky_model_values.shape):
                    shading_test_values = interpolate_to_fit(shading_test_values, sky_model_values)
                else:
                    sky_model_values = interpolate_to_fit(sky_model_values, shading_test_values)

            # Transposing the sky_model_values if it is not in the same shape as shading_test_values
            if shading_test_values.shape != sky_model_values.shape:
                sky_model_values = sky_model_values.T

            result = (1 - shading_test_values ) * sky_model_values
            
            # calculate the maximum, but avoid zero
            max_val = np.max(result)
            if max_val < 1e-10:
                max_val = 1e-10
            result /= max_val

            result = np.nan_to_num(result)  # replace any NaN or Inf values

            fig = plt.figure()
            altitude = np.linspace(0, 90, result.shape[1]+1)
            azimuth = np.linspace(0, 2*np.pi, result.shape[0]+1)
            r, th = np.meshgrid(altitude, azimuth)

            ax1 = plt.subplot(projection="polar")
            plt.grid(False)
            im = plt.pcolormesh(th, r, result, cmap=cm.viridis , vmin=0, vmax=1) #cmap=cm.gray_r

            v1 = np.linspace(0, 1, 11)
            cbar = plt.colorbar(im,ticks=v1)
            cbar.ax.set_yticklabels(["{:4.2f}".format(i) for i in v1])

            plt.thetagrids([i*15 for i in range(0,24)])
            ax1.set_theta_direction(-1)
            ax1.set_theta_zero_location("N")

            plt.rgrids([i*10 for i in range(0,10)])
            ax1.set_rlim(bottom=90, top=0)

            output_file = destination / f'Shading_mask_{shading_name}_{hour}H.png'
            plt.savefig(output_file)
            plt.close(fig)

def plotShadingMaskDir(directory_path, destination):
    for root, dirs, files in os.walk(directory_path):
        for file in files:
            if file.endswith(".csv"):
                plotShadingMask(os.path.join(root,file), destination)

parser = argparse.ArgumentParser()
parser.add_argument("--dir_path", type=Path, help='Directory containing the shading mask matrices')
parser.add_argument("--destination",type=Path, required=True, help='Destination of the shading mask images')

p = parser.parse_args()

if p.dir_path and p.dir_path.exists():
    plotShadingMaskDir(p.dir_path, p.destination)
