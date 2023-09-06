try:
    from ._shadingmask import *

    _smclass={
            '2-3': ShadingMask_2D3D,
            '3-3': ShadingMask_3D3D        
        }    
except ImportError as e:
    print('Import feelpp.toolboxes.fluid failed: Feel++ Toolbox Fluid is not available')
    pass


def shadingmask(mesh,jsonfile,intervalsAzimuth=72,intervalsAltitude=10):
    realdim = mesh.realDimension()
    dim = mesh.dimension()

    key = str(dim) + '-' +  str(realdim) 


    return _smclass[key](mesh=mesh,json=jsonfile,intervalsAzimuth=intervalsAzimuth,intervalsAltitude=intervalsAltitude)

