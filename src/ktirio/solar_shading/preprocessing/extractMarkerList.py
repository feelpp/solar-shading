# Extract the list of surface and volume markers from a .msh file

import os
from pathlib import Path
import argparse

def extractMSHmarkers(input_file):
    """
    Extract the list of surface and volume markers from a .msh file
    """
    outputfile_surf = os.path.splitext(input_file)[0]+"_surf.txt" # file with all surface markers
    outputfile_vol = os.path.splitext(input_file)[0]+"_vol.txt" # file with all volume markers
    outputfile_surf_id = os.path.splitext(input_file)[0]+"_surfid.txt" # file with all building-only surface markers (not face by face nor soil)
    outputfile_face_id = os.path.splitext(input_file)[0]+"_faceid.txt" # file with all building-only surface markers (not face by face nor soil)
    with open(input_file,"r") as inputfile, open(outputfile_surf,"w") as outputfilesurf, \
         open(outputfile_vol,"w") as outputfilevol, open(outputfile_surf_id,"w") as outputfilesurfid, \
         open(outputfile_face_id,"w") as outputfilefaceid :
        copyLine = False
        for line in inputfile:
            if line.strip() == "$PhysicalNames":
                copyLine=True
                continue
            elif line.strip() == "$EndPhysicalNames":
                copyLine=False
                break
            elif copyLine:
                split_line=line.split()
                if(split_line[0] == "2"):
                    outputfilesurf.write(split_line[2].replace('"', '')+"\n")
                    if "face" in split_line[2]:
                        outputfilefaceid.write(split_line[2].replace('"', '')+"\n")
                    elif "face" not in split_line[2] and "terrain" not in split_line[2] :
                        outputfilesurfid.write(split_line[2].replace('"', '')+"\n")
                elif(split_line[0] == "3"):
                    outputfilevol.write(split_line[2].replace('"', '')+"\n")                    

        

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh_file_path", type=Path, help='Insert the name of the mesh file from which surface and volume markers need to be extracted')

    p = parser.parse_args()

    if( p.mesh_file_path  and p.mesh_file_path.exists() ):
        extractMSHmarkers(p.mesh_file_path)