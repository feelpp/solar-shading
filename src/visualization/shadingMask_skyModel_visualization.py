import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from pathlib import Path
import argparse
import os
import re
from scipy.ndimage import zoom

# spline interpolation => utiliser fonciton scipy vue avec vigon
def resample_array_to_match_shape(array, target_shape):
    zoom_ratio = [t/s for t, s in zip(target_shape, array.shape)]
    return zoom(array, zoom_ratio, order=5)

def calculate_and_store_max(csv_filename):
    shading_test_values = np.genfromtxt(csv_filename, delimiter=',')
    shading_name = Path(csv_filename).stem
    key_part = re.search('SM_Matrix_(.*)_.*', shading_name).group(1)

    sky_models_directory = Path(csv_filename).parent / ".." / "skyModels"
    max_value = 1e-10
    for hour in range(5, 20):
        sky_model_filename = sky_models_directory / f"SkyModel_Matrix_{key_part}_{hour}H.csv"
        if sky_model_filename.exists():
            sky_model_values = np.genfromtxt(sky_model_filename, delimiter=',')
            if shading_test_values.shape != sky_model_values.shape:
                shading_test_values = resample_array_to_match_shape(shading_test_values, sky_model_values.shape)
            result = (1 - shading_test_values ) * sky_model_values
            max_val = np.max(result)
            if max_val > max_value:
                max_value = max_val
    return max_value

def plotShadingMask(csv_filename, destination, max_value):
    shading_test_values = np.genfromtxt(csv_filename, delimiter=',')
    shading_name = Path(csv_filename).stem
    key_part = re.search('SM_Matrix_(.*)_.*', shading_name).group(1)

    sky_models_directory = Path(csv_filename).parent / ".." / "skyModels"
    for hour in range(5, 20):
        sky_model_filename = sky_models_directory / f"SkyModel_Matrix_{key_part}_{hour}H.csv"
        if sky_model_filename.exists():
            sky_model_values = np.genfromtxt(sky_model_filename, delimiter=',')
            if shading_test_values.shape != sky_model_values.shape:
                shading_test_values = resample_array_to_match_shape(shading_test_values, sky_model_values.shape)
            result = (1 - shading_test_values ) * sky_model_values
            result = np.nan_to_num(result)  # replace any NaN or Inf values

            fig = plt.figure()
            altitude = np.linspace(0, 90, result.shape[1]+1)
            azimuth = np.linspace(0, 2*np.pi, result.shape[0]+1)
            r, th = np.meshgrid(altitude, azimuth)

            ax1 = plt.subplot(projection="polar")
            plt.grid(False)
            im = plt.pcolormesh(th, r, result, cmap=cm.viridis , vmin=0, vmax=max_value) 

            v1 = np.linspace(0, max_value, 11)
            cbar = plt.colorbar(im,ticks=v1)
            cbar.ax.set_yticklabels(["{:d}".format(int(i)) for i in v1])
            cbar.set_label('Illuminance [lux]', rotation=270, labelpad=20)

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
                csv_file_path = os.path.join(root,file)
                max_value = calculate_and_store_max(csv_file_path)
                plotShadingMask(csv_file_path, destination, max_value)

parser = argparse.ArgumentParser()
parser.add_argument("--dir_path", type=Path, required=True, help='Directory containing the shading mask matrices')
parser.add_argument("--destination",type=Path, default= Path("results/images") , help='Destination of the shading mask images')
# parser.add_argument("--RNG", type=str, default="_STD", help='RNG name to be added to the output')

p = parser.parse_args()

if p.dir_path and p.dir_path.exists():
    plotShadingMaskDir(p.dir_path, p.destination)
