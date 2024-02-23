from ktirio.shadingmask import *
import json
import feelpp as fpp
import sys
import os

def test_volumeMask(init_feelpp):

    fpp.Environment.changeRepository(
        directory="shading-tests/volume")

    mesh_file_path = "../../src/cases/example_shading_mask/example_shading_mask.geo"
    json_file_path = "../../src/cases/example_shading_mask/example_shading_mask.json"

    mesh_file_path=os.path.join(os.path.dirname(__file__),mesh_file_path)

    # Load the mesh
    mesh=fpp.load( fpp.mesh(dim=3,realdim=3) , str(mesh_file_path) )

    # Load the JSON file configuring the mask computation
    json_file_path = os.path.join(os.path.dirname(__file__),json_file_path)
    with open(json_file_path,'r') as f:
        jsonfile = json.load(f)
        print(jsonfile)

    # Launch the underlying C++ application on 1 thread only
    jsonfile["Nthreads"]=1

    SM = shadingmask(mesh,jsonfile)

    SM.computeMasks()

    assert(True)
