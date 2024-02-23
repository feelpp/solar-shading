# Test for shading mask diagram plotting

from ktirio.shadingmask import visualization as visu
from pathlib import Path
import argparse

if __name__ == '__main__':
    # Read the command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--file_path", type=Path, help='Insert the name of the csv file containing the shading mask matrix')
    parser.add_argument("--dir_path", type=Path, help='Insert the name of the directory containing the csv file with the shading mask matrix')
    parser.add_argument("--destination",type=Path, help='Destination of the shading mask images. Default: location of the python script')

    p, unknown = parser.parse_known_args()

    if(p.dir_path and p.dir_path.exists()):
        visu.plotShadingMaskDir(p.dir_path)

    if( p.file_path  and p.file_path.exists() ):
        visu.plotShadingMask(p.file_path)