from shadingmask import *
import json
import feelpp as fpp
import sys
import os

# def test_Init():
#     e = fpp.Environment(
#         sys.argv,
#         config=fpp.globalRepository("surfaceMask"))

def test_surfaceMask(init_feelpp):
    
    fpp.Environment.changeRepository(
        directory="shading-tests/surface")
    
    mesh_file_path = "../../src/cases/example_shading_mask/skin_mesh_example.geo"
    json_file_path = "../../src/cases/example_shading_mask/example_shading_mask_surface.json"

    mesh_file_path=os.path.join(os.path.dirname(__file__),mesh_file_path)

    # Load the mesh
    mesh=fpp.load( fpp.mesh(dim=2,realdim=3) , str(mesh_file_path) )

    # Load the JSON file configuring the mask computation
    json_file_path = os.path.join(os.path.dirname(__file__),json_file_path)
    with open(json_file_path,'r') as f:
        jsonfile = json.load(f)
        print(jsonfile)

    # Launch the underlying C++ application on 1 thread only
    jsonfile["Nthreads"]=1

    # Read the file with the list of buildings
    fileSurfaces = os.path.basename(jsonfile["Buildings"]["fileSurfaces"])
    jsonfile["Buildings"]["fileSurfaces"] = os.path.join(os.path.dirname(json_file_path),fileSurfaces)


    SM = shadingmask(mesh,jsonfile)

    SM.computeMasks()
