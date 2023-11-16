import json
import pandas as pd
import numpy as np
import os
import glob

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import math

from scipy.interpolate import make_interp_spline
from scipy.interpolate import interp1d

#Directory="C:/YOO/TMP/001/SRC/np11"

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

for k in range(NbFiles):
    FICH = open(FilesList[k])
    data = json.load(FICH)
    FICH.close()

    NraysPerElement              =data['shadingMask']['NraysPerElement']
    Nthreads                     =data['shadingMask']['Nthreads']
    nBuildingFaces               =data['shadingMask']['nBuildingFaces']
    nMarkers                     =data['shadingMask']['nMarkers']

    BVH_building_time            =data['shadingMask']['Timer']['BVH_building_time']
    DataStructures_building_time =data['shadingMask']['Timer']['DataStructures_building_time']
    MaskCSVsaving                =data['shadingMask']['Timer']['MaskCSVsaving']
    MaskComputation              =data['shadingMask']['Timer']['MaskComputation']

    MaskComputationGenerate      =data['shadingMask']['Timer Specx']['MaskComputation Generate']
    MaskSavingGenerate           =data['shadingMask']['Timer Specx']['MaskSaving Generate']
    MaskSavingTask1              =data['shadingMask']['Timer Specx']['Compute Task 1']
    #TimeCTRLMasks                =data['shadingMask']['Timer Specx']['Time CTRL Masks']
    TimeSaveParallelMode         =data['shadingMask']['Timer Specx']['Time save parallel mode']
    #TimeSaveSequentialMode       =data['shadingMask']['Timer Specx']['Time save sequential mode']
    
    DeltaTime                    =data['shadingMask']['Timestamp']['Delta Time']
    

    TabNthreads.append(Nthreads)
    TabNraysPerElement.append(NraysPerElement)
    TabMaskComputation.append(MaskComputation)
    TabDeltaTime.append(DeltaTime)

    Value01=data['shadingMask']['Nthreads']
    Value02=data['shadingMask']['NraysPerElement']
    Value03=data['shadingMask']['nBuildingFaces']
    Value04=data['shadingMask']['nMarkers']
    
    Value05=data['shadingMask']['Timer']['BVH_building_time']
    Value06=data['shadingMask']['Timer']['DataStructures_building_time']
    Value07=data['shadingMask']['Timer']['MaskCSVsaving']
    Value08=data['shadingMask']['Timer']['MaskComputation']

    Value09=data['shadingMask']['Timer Specx']['MaskComputation Generate']
    Value10=data['shadingMask']['Timer Specx']['MaskSaving Generate']
    Value11=data['shadingMask']['Timer Specx']['Compute Task 1']
    #Value12=data['shadingMask']['Timer Specx']['Time CTRL Masks']
    Value13=data['shadingMask']['Timer Specx']['Time save parallel mode']
    #Value14=data['shadingMask']['Timer Specx']['Time save sequential mode']
    
    Value15=data['shadingMask']['Timestamp']['Delta Time']

    Line=[Value01,Value02,Value03,Value04,Value05,
          Value06,Value07,Value08,Value09,Value10,
          Value11,Value13,Value15]

    if (k==0):
        TabData=np.array(Line)       
    else:
        #print("Append")
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
,'MaskComputation Generate'
,'MaskSaving Generate'
,'Compute Task 1'
#,'Time CTRL Masks'
,'Time save parallel mode'
#,'Time save sequential mode'])
,'DeltaTime'])

df.to_csv(PathSave+'/DATA.csv')  


df = pd.read_csv(PathSave+'/DATA.csv')

print(df.to_string()) 




if (1==1):
    fig1 = plt.figure("Figure 1")
    
    subdef=df.loc[df['Nthreads'] == 1]
    x=subdef["NraysPerElement"]
    y=subdef["MaskComputation"]
    plt.plot(x,y, 'bo-',label='1 Thread')
    
    subdef=df.loc[df['Nthreads'] == 10]
    x=subdef["NraysPerElement"]
    y=subdef["MaskComputation"]
    plt.plot(x,y, 'go-',label='10 Thread')
    
    subdef=df.loc[df['Nthreads'] == 20]
    x=subdef["NraysPerElement"]
    y=subdef["MaskComputation"]
    plt.plot(x,y, 'ro-',label='20 Thread')
    
    subdef=df.loc[df['Nthreads'] == 30]
    x=subdef["NraysPerElement"]
    y=subdef["MaskComputation"]
    plt.plot(x,y, 'co-',label='30 Thread')
    
    subdef=df.loc[df['Nthreads'] == 40]
    x=subdef["NraysPerElement"]
    y=subdef["MaskComputation"]
    plt.plot(x,y, 'mo-',label='40 Thread')
    
    subdef=df.loc[df['Nthreads'] == 50]
    x=subdef["NraysPerElement"]
    y=subdef["MaskComputation"]
    plt.plot(x,y, 'yo-',label='50 Thread')
    
    subdef=df.loc[df['Nthreads'] == 60]
    x=subdef["NraysPerElement"]
    y=subdef["MaskComputation"]
    plt.plot(x,y, 'ko-',label='60 Thread')
    
    subdef=df.loc[df['Nthreads'] == 70]
    x=subdef["NraysPerElement"]
    y=subdef["MaskComputation"]
    plt.plot(x,y, 'gx-',label='70 Thread')
    
    subdef=df.loc[df['Nthreads'] == 80]
    x=subdef["NraysPerElement"]
    y=subdef["MaskComputation"]
    plt.plot(x,y, 'yx-',label='80 Thread')
    
    
    #plt.hist(x, bins=4, edgecolor='black')
    plt.xlabel('Nb Rays')
    plt.ylabel('Time Computation')
    plt.title('Solar Shading Performance')
    plt.legend()
    plt.savefig(PathSave+'/GR_NbRayPerf.png', dpi=300, bbox_inches='tight')
    plt.show()



if (1==1):
    fig0 = plt.figure("Figure 0")
    

    subdef=df.loc[df['NraysPerElement'] ==30000]
    x=subdef['Nthreads']
    y=subdef["MaskComputation"]
    plt.plot(x,y, 'ko',label='30000 Rays') 
    

    subdef=df.loc[df['NraysPerElement'] ==50000]
    x=subdef['Nthreads']
    y=subdef["MaskComputation"]
    plt.plot(x,y, 'yx',label='50000 Rays') 
    #X_Y_Spline = make_interp_spline(sorted(x), y)
    #X_ = np.linspace(x.min(), x.max(), 500)
    #Y_ = X_Y_Spline(X_)   
    #plt.plot(X_, Y_)
    
    #cubic_interpolation_model =interp1d(sorted(x), y, kind = "cubic")
    #Y_=cubic_interpolation_model(X_)
    #plt.plot(X_, Y_)
    
    #plt.hist(x, bins=4, edgecolor='black')
    plt.xlabel('Nb Threads')
    plt.ylabel('Time Computation')
    plt.title('Solar Shading Performance')
    plt.legend()
    plt.savefig(PathSave+'/GR_NbRayPerf0.png', dpi=300, bbox_inches='tight')
    plt.show()




#==============================================================================


for k1 in range(1,11):
    #SpeedUP Part
    NbRay=10000*k1
    NbRay=10000
    NumThreadT1=1
    subdef=df.loc[df['NraysPerElement'] == NbRay]
    print(subdef)
    T1=df.loc[ (df['NraysPerElement']==NbRay) & (df['Nthreads']==NumThreadT1) ]["MaskComputation"]
    print("T1="+str(T1))
    print("1/T1="+str(1/T1))
    
    Nn=subdef["Nthreads"]
    Tn=subdef["MaskComputation"]
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
    
    
    
    fig2 = plt.figure("Figure 2")
    x=subdef["Nthreads"]
    y=SpeedUP
    plt.plot(x,y, 'ro',label='NbRay='+str(NbRay))
    plt.plot([100,100])
    plt.xlabel('Numbers Of Treads')
    plt.ylabel('SpeedUp')
    plt.title('Solar Shading Performance')
    plt.grid()
    plt.legend()
    #X_Y_Spline = make_interp_spline(sorted(x), y)
    #X_ = np.linspace(x.min(), x.max(), 500)
    #Y_ = X_Y_Spline(X_)
    #plt.plot(X_, Y_)
    plt.savefig(PathSave+'/GR_SpeedUP_'+str(NbRay)+'.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    
    fig3 = plt.figure("Figure 3")
    x=subdef["Nthreads"]
    y=Amdalh
    plt.plot(x,y, 'ro',label='NbRay='+str(NbRay))
    plt.xlabel('Numbers Of Treads')
    plt.ylabel('Amdalh Law')
    plt.title('Solar Shading Performance')
    plt.grid()
    plt.legend()
    #X_Y_Spline = make_interp_spline(sorted(x), y)
    #X_ = np.linspace(x.min(), x.max(), 500)
    #Y_ = X_Y_Spline(X_)
    #plt.plot(X_, Y_)
    plt.savefig(PathSave+'/GR_AmdalhLaw_'+str(NbRay)+'.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    fig4 = plt.figure("Figure 4")
    x=subdef["Nthreads"]
    y=Efficiency
    plt.plot(x,y, 'ro',label='NbRay='+str(NbRay))
    plt.xlabel('Numbers Of Treads')
    plt.ylabel('Efficiency')
    plt.title('Solar Shading Performance')
    plt.grid()
    plt.legend()
    #X_Y_Spline = make_interp_spline(sorted(x), y)
    #X_ = np.linspace(x.min(), x.max(), 500)
    #Y_ = X_Y_Spline(X_)
    #plt.plot(X_, Y_)
    plt.savefig(PathSave+'/GR_Efficiency_'+str(NbRay)+'.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    
    fig5 = plt.figure("Figure 5")
    x=subdef["Nthreads"]
    y=Amdalh
    plt.xscale("log")
    #plt.yscale("log")
    plt.plot(x,y, 'mo',label='NbRay='+str(NbRay))
    plt.xlabel('Numbers Of Treads')
    plt.ylabel('Log SpeedUP')
    plt.title('Solar Shading Performance')
    plt.grid()
    plt.legend()
    
    #X_Y_Spline = make_interp_spline(sorted(x), y)
    #X_ = np.linspace(x.min(), x.max(), 500)
    #Y_ = X_Y_Spline(X_)
    #plt.plot(X_, Y_)
    
    
    plt.savefig(PathSave+'/GR_LogSpeedUP_'+str(NbRay)+'.png', dpi=300, bbox_inches='tight')
    plt.show()

#==============================================================================





#==============================================================================


for k1 in range(1,11):
    #SpeedUP Part
    NbRay=5000*k1
    NbRay=30000
    NbRay=50000
    NbRay=10000
    NumThreadT1=1
    subdef=df.loc[df['NraysPerElement'] == NbRay]
    print(subdef)
    T1=df.loc[ (df['NraysPerElement']==NbRay) & (df['Nthreads']==NumThreadT1) ]["MaskComputation"]
    print("T1="+str(T1))
    print("1/T1="+str(1/T1))
    
    Nn=subdef["Nthreads"]
    Tn=subdef["MaskComputation"]
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
    plt.subplot(311)
    x=subdef["Nthreads"]
    y=SpeedUP
    plt.plot(x,y, 'ro',label='NbRay='+str(NbRay))
    #plt.plot([100,100])
    plt.xlabel('Numbers Of Treads')
    plt.ylabel('SpeedUp')
    plt.title('Solar Shading Performance')
    plt.grid()
    plt.legend()
    #X_Y_Spline = make_interp_spline(sorted(x), y)
    #X_ = np.linspace(x.min(), x.max(), 500)
    #Y_ = X_Y_Spline(X_)
    #plt.plot(X_, Y_)

    
    plt.subplot(312)
    x=subdef["Nthreads"]
    y=Amdalh
    plt.plot(x,y, 'ro',label='NbRay='+str(NbRay))
    plt.xlabel('Numbers Of Treads')
    plt.ylabel('Amdalh Law')
    #plt.title('Solar Shading Performance')
    plt.grid()
    plt.legend()
    #X_Y_Spline = make_interp_spline(sorted(x), y)
    #X_ = np.linspace(x.min(), x.max(), 500)
    #Y_ = X_Y_Spline(X_)
    #plt.plot(X_, Y_)

    plt.subplot(313)
    x=subdef["Nthreads"]
    y=Efficiency
    plt.plot(x,y, 'ro',label='NbRay='+str(NbRay))
    plt.xlabel('Numbers Of Treads')
    plt.ylabel('Efficiency')
    #plt.title('Solar Shading Performance')
    plt.grid()
    plt.legend()
    #X_Y_Spline = make_interp_spline(sorted(x), y)
    #X_ = np.linspace(x.min(), x.max(), 500)
    #Y_ = X_Y_Spline(X_)
    #plt.plot(X_, Y_)


    
    
    plt.savefig(PathSave+'/GR_SYNTHESE_'+str(NbRay)+'.png', dpi=300, bbox_inches='tight')
    plt.show()

#==============================================================================










