import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pathlib import Path
import pandas as pd
import pvlib
from pvlib import location

# This script is currently not used in the project since the path isn't correct, but it is kept for future reference and improvements

def plotShadingMask(csv_filename, latitude, longitude, tz):

    # Define a location for pvlib
    site = location.Location(latitude, longitude, tz=tz)

    # Calculate the solar position for every day of the year at sunrise and sunset
    times = pd.date_range('2023-01-01 00:00:00', '2023-12-31 23:59:59', freq='1D', tz=site.tz)
    
    # calculate sunrise, sunset and transit times
    sunrise, sunset, _ = pvlib.solarposition.sun_rise_set_transit_spa(times, latitude, longitude)
    
    fig = plt.figure()

    # The CSV file must contain a comma-separated matrix of size
    # azimuth x altitude
    shading_test_values = np.genfromtxt(csv_filename, delimiter=',')

    # azimuth-latitude coordinates
    altitude = np.linspace(0, 90, shading_test_values.shape[1]+1)
    azimuth = np.linspace(0, 2*np.pi, shading_test_values.shape[0]+1)

    # Coordinate grid for the pcolormesh
    r, th = np.meshgrid(altitude, azimuth)

    ax1 = plt.subplot(projection="polar")
    plt.grid(False)

    im = plt.pcolormesh(th, r, shading_test_values, cmap=cm.gray_r , vmin=0, vmax=1)

    # Plot solar path
    for day in times:
        # Convert fractional days to timestamps
        sunrise_time = pd.Timestamp(day.date()) + pd.Timedelta(days=sunrise[day.dayofyear - 1])
        sunset_time = pd.Timestamp(day.date()) + pd.Timedelta(days=sunset[day.dayofyear - 1])
        
        # Generate time range for each day
        times_day = pd.date_range(start=sunrise_time, end=sunset_time, freq='5min')
        solpos_day = site.get_solarposition(times_day)
        az_rad = np.radians(360 - solpos_day.azimuth)
        alt_rad = 90 - solpos_day.apparent_elevation
        ax1.plot(az_rad, alt_rad, color='r')

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

    # assuming `destination` is a variable containing your destination path
    output_file = os.path.join(p.destination, 'Shading_mask'+Path(csv_filename).stem+'.png')
    plt.savefig(output_file)

    print("Shading mask image saved in: ", output_file)

    plt.close(fig)

def plotShadingMaskDir(directory_path, latitude, longitude, tz):
    for root, dirs, files in os.walk(directory_path):
        for file in files:
            if file.endswith(".csv"):
                plotShadingMask(os.path.join(root,file), latitude, longitude, tz)

# Read the command line
parser = argparse.ArgumentParser()
parser.add_argument("--file_path", type=Path, help='Insert the name of the csv file containing the shading mask matrix')
parser.add_argument("--dir_path", type=Path, help='Insert the name of the directory containing the csv file with the shading mask matrix')
parser.add_argument("--destination",type=Path, help='Destination of the shading mask images. Default: location of the python script')
parser.add_argument("--latitude", type=float, required=True, help='Latitude of the location')
parser.add_argument("--longitude", type=float, required=True, help='Longitude of the location')
parser.add_argument("--timezone", type=str, required=True, help='Timezone of the location')

p = parser.parse_args()

if(p.dir_path and p.dir_path.exists()):
    plotShadingMaskDir(p.dir_path, p.latitude, p.longitude, p.timezone)

if(p.file_path and p.file_path.exists()):
    plotShadingMask(p.file_path, p.latitude, p.longitude, p.timezone)
