import json
import pandas as pd
import numpy as np
#import os
import glob

#from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
#from matplotlib import cm

import math


#Directory="C:/YOO/TMP/001/SRC/np9"
Directory="/data/scratch/lemoinep/mon_home/feelppdb/str0/np_1"
PathSave=Directory

FilesList=glob.glob(Directory+"/shadingmask_metadata*.json")
NbFiles=len(FilesList)
print("NbFiles="+str(NbFiles))

TabNthreads = []
TabNraysPerElement = []
TabMaskComputation = []
TabDeltaTime = []


TabData=[]  
#for i in data['Nthreads']:
#    print(i)

NumOption=1

if (NumOption==1):

    for k in range(NbFiles):
        FICH = open(FilesList[k])
        data = json.load(FICH)
        FICH.close()
 
        Value01=data['shadingMask']['Nthreads']
        Value02=data['shadingMask']['NraysPerElement']
        Value03=data['shadingMask']['nBuildingFaces']
        Value04=data['shadingMask']['nMarkers']
        
        Value05=data['shadingMask']['Timer']['BVH_building_time']
        Value06=data['shadingMask']['Timer']['DataStructures_building_time']
        Value07=data['shadingMask']['Timer']['MaskCSVsaving']
        Value08=data['shadingMask']['Timer']['MaskComputation']
        Value09=data['shadingMask']['Timestamp']['Delta Time']
    
        Line=[Value01,Value02,Value03,Value04,Value05,
              Value06,Value07,Value08,Value09]
    
        if (k==0):
            TabData=np.array(Line)       
        else:
            TabData=np.vstack((TabData,Line))

    print(TabData.shape)

    df = pd.DataFrame(TabData, 
    index =  [str(n) for n in range(0,NbFiles)], 
    columns = 
    ['Nthreads'
    ,'NraysPerElement'
    ,'nBuildingFaces'
    ,'nMarkers'
    ,'BVH_building_time'
    ,'DataStructures_building_time'
    ,'MaskCSVsaving'
    ,'MaskComputation'
    ,'DeltaTime'])
    df.to_csv(PathSave+'/DATA.csv')  




if (NumOption==0):

    for k in range(NbFiles):
        FICH = open(FilesList[k])
        data = json.load(FICH)
        FICH.close()
    
        Value01=data['shadingMask']['Nthreads']
        Value02=data['shadingMask']['NraysPerElement']
        Value03=data['shadingMask']['nBuildings']

        
        Value04=data['shadingMask']['Timer']['BVHs_total_building_time']
        Value05=data['shadingMask']['Timer']['MaskComputation']
    
        Line=[Value01,Value02,Value03,Value04,Value05]
    
        if (k==0):
            TabData=np.array(Line)       
        else:
            TabData=np.vstack((TabData,Line))

    print(TabData.shape)

    df = pd.DataFrame(TabData, 
    index =  [str(n) for n in range(0,NbFiles)], 
    columns = 
    ['Nthreads'
    ,'NraysPerElement'
    ,'nBuildingFaces'
    ,'BVH_building_time'
    ,'MaskComputation'])
    df.to_csv(PathSave+'/DATA.csv')  


if (1==1):
    df = pd.read_csv(PathSave+'/DATA.csv')
    print(df.to_string()) 
    T1=df.loc[ (df['NraysPerElement']==1000) & (df['Nthreads']==1) ]["MaskComputation"]
    T1=df.loc[ (df['NraysPerElement']==1000) & (df['Nthreads']==1) ]["DeltaTime"]
    print("T1="+str(T1))
    print("1/T1="+str(1/T1))
    
    Nn=df["Nthreads"]
    Tn=df["MaskComputation"]
    Tn=df["DeltaTime"]
    NbTn=Tn.shape[0]
    print("nb Tn="+str(NbTn))
    
    print("Tn=")
    print(str(Tn))
    print("")
    
    SpeedUP=[]
    Efficiency=[]
    Amdalh=[]
    Amdalh=[]
    LogSpeedUP=[]
    
    for k in range(0,NbTn):
        #print(k)
        Index=Tn.index[k]
        SpeedUP.append(T1/Tn[Index])
        Efficiency.append(T1/(Tn[Index]*Nn[Index]))
        Amdalh.append(1/(T1+Tn[Index]/Nn[Index]))
        LogSpeedUP.append(math.log(T1)-math.log(Tn[Index]))
    
    print()
    print(str(SpeedUP))
    print("")
    
    fig1 = plt.figure("Figure 1")
    x=df["Nthreads"]
    y=df["DeltaTime"]
    plt.plot(x,y, 'ro')
    plt.xlabel('Numbers Of Treads')
    plt.ylabel('Time')
    plt.title('Solar Shading Performance')
    plt.grid()
    plt.legend()
    plt.savefig(PathSave+'/GR_Performance_Scale.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    fig2 = plt.figure("Figure 2")
    x=df["Nthreads"]
    y=SpeedUP
    plt.plot(x,y, 'ro')
    plt.xlabel('Numbers Of Treads')
    plt.ylabel('SpeedUp')
    plt.title('Solar Shading Performance')
    plt.grid()
    plt.legend()
    plt.savefig(PathSave+'/GR_SpeedUP_Scale.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    fig3 = plt.figure("Figure 3")
    x=df["Nthreads"]
    y=Amdalh
    plt.plot(x,y, 'ro')
    plt.xlabel('Numbers Of Treads')
    plt.ylabel('Amdalh Law')
    plt.title('Solar Shading Performance')
    plt.grid()
    plt.legend()
    plt.savefig(PathSave+'/GR_AmdalhLaw_Scale.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    fig4 = plt.figure("Figure 4")
    x=df["Nthreads"]
    y=Efficiency
    plt.plot(x,y, 'ro')
    plt.xlabel('Numbers Of Treads')
    plt.ylabel('Efficiency')
    plt.title('Solar Shading Performance')
    plt.grid()
    plt.legend()
    plt.savefig(PathSave+'/GR_Efficiency_Scale.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    
    
    
    
    
    

    