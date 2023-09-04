# Test for shading mask diagram plotting

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from pathlib import Path
import argparse
import os
import plotly.graph_objs as go

def plotShadingMask(csv_filename,save=True):

    # The CSV file must contain a comma-separated matrix of size
    # azimuth x altitude
    shading_test_values = np.genfromtxt(csv_filename, delimiter=',')
    
    if save:
        fig = plt.figure()

        # azimuth-latitude coordinates
        altitude = np.linspace(0, 90, shading_test_values.shape[1]+1)
        azimuth = np.linspace(0, 2*np.pi, shading_test_values.shape[0]+1)
        
        # Coordinate grid for the pcolormesh
        r, th = np.meshgrid(altitude, azimuth)

        ax1 = plt.subplot(projection="polar")
        plt.grid(False)

        im = plt.pcolormesh(th, r, shading_test_values, cmap=cm.gray_r , vmin=0, vmax=1)

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

        plt.savefig('Shading_mask'+Path(csv_filename).stem+'.png')
        plt.close()
    else:

        hovertemplate = ('Altitude: %{r}<br>'
                    'Azimuth: %{theta}<br>'
                    'Solar mask: %{customdata}<br>'
                    '<extra></extra>')
        
        # azimuth-altitude coordinates
        r, theta = np.mgrid[0:90:shading_test_values.shape[1]+1j, 0:360:shading_test_values.shape[0]+1j]
        
        # Plotly polar representation of the shading mask
        fig = go.Figure(go.Barpolar(
                r=r.ravel(),
                theta=theta.ravel(),
                hovertemplate=hovertemplate,
                customdata=shading_test_values.ravel(order='F'),
                marker=dict(
                    colorscale='gray',
                    showscale=True,
                    color=shading_test_values.ravel(order='F'),
                    reversescale=True)
                )
            )
        
        fig.update_layout(
                title='',
                polar=dict(
                    angularaxis=dict(tickvals=np.arange(0, 360, 15),
                                    direction='clockwise'),
                    radialaxis=dict(angle=60,
                                    tickvals=[],
                                    autorange="reversed"),
                    bargap=0
                ),
                
            )
        fig.show()


def plotShadingMaskDir(directory_path):
    for root, dirs, files in os.walk(directory_path):
        for file in files:
            if file.endswith(".csv"):
                plotShadingMask(os.path.join(root,file))

if __name__ == '__main__':
    # Read the command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--file_path", type=Path, help='Insert the name of the csv file containing the shading mask matrix')
    parser.add_argument("--dir_path", type=Path, help='Insert the name of the directory containing the csv file with the shading mask matrix')
    parser.add_argument("--destination",type=Path, help='Destination of the shading mask images. Default: location of the python script')

    p, unknown = parser.parse_known_args()

    if(p.dir_path and p.dir_path.exists()):
        plotShadingMaskDir(p.dir_path)

    if( p.file_path  and p.file_path.exists() ):
        plotShadingMask(p.file_path)