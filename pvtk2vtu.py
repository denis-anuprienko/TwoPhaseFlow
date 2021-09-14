# trace generated using paraview version 5.9.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

dir = ''
N = 1

for i in range(N):
    # create a new 'Partitioned Legacy VTK Reader'
    pvtk = PartitionedLegacyVTKReader(registrationName='pvtk', FileName=dir+'/sol'+str(i) +'.pvtk')
    # save data
    SaveData(dir+'/tmp.vtu', proxy=pvtk, PointDataArrays=[],CellDataArrays=['Liquid_Saturation'],DataMode='Binary',
    CompressorType='LZMA',
    CompressionLevel='9')
    Delete(pvtk)
    del pvtk
    
    tmp = XMLUnstructuredGridReader(registrationName='psol0.vtu', FileName=[dir+'/tmp.vtu'])
    
    tmp.CellArrayStatus = ['Liquid_Saturation']
    tmp.PointArrayStatus = []

    # save data
    SaveData(dir+'/psol'+str(i)+'.vtu', proxy=tmp, CellDataArrays=['Liquid_Saturation'])
    Delete(tmp)
    del tmp
    
