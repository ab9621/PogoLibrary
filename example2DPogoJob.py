# -*- coding: utf-8 -*-
'''
Author: Alexander Ballisat

This is an example of generating a 2D model, running it, analysing the
results and animating the field data. It is also an example of how to set up
Pogo for batch simulations.

For a more detailed instructions on how to set up Pogo please refer to the
example3DPogoJob script.

For 2D problems, the triangular meshing tool Triangle is regularly used. This
is needed for this script but can be downloaded from:
https://www.cs.cmu.edu/~quake/triangle.html

Hopefully this script is reasonably self explanatory. If you do have any 
questions or find any bugs please email:
alexander.ballisat@bristol.ac.uk
'''

# The necessary imports of common scripts
import numpy as np
import matplotlib.pyplot as plt
import subprocess # this is a standard library of functions that allows 
                  # system calls to be made as if in a command prompt

# The Pogo specific imports
import savePogoInputFile as sPIF
import GeometryFunctions as gF
import loadPogoHistFile as lPH
import pogoFunctions as pF
import pogoJobSubmission as pJS
import SingleTraceAnalysis as sTA

'''
Define a general function for generating a 2D job
'''
def pogo2DJob(path,
              jobName='example2DPogoJob',
              xSize = 50E-3,
              ySize = 50E-3,
              f0=5E6,
              nElementsPerWavelength=10,
              transducerWidth=5E-3,
              inputSignal='gauss',
              E = 70E9,
              nu = 0.34,
              rho = 2700.,
              transducerCentre=None,
              modelTime = None,
              modelTimeStep = None):
    '''
    General function for generating and running a 2D model for an isotropic
    solid.
    
    NEEDS MORE COMMENTING
    '''

    gF.write2DRectanglePolyFile(xSize, ySize, r'{}\{}'.format(path, jobName))
    
    vL, vP = pF.waveVelocity(E, nu, rho)
    print 'Longitudinal wave velocity = {}'.format(vL)
    # Step 2 - Set some properties of the model
    if modelTime == None:
        modelTime = 3.0*ySize/vL
    
    targetMeshSize = (vL/(f0*nElementsPerWavelength))**2/2. #convert the time to a distance
    
    subprocess.call('triangle -a{:.12f} {}.poly'.format(targetMeshSize, jobName),
                    cwd=path)
    
    nodes = pF.loadNodeFile('{}\{}.1.node'.format(path, jobName), nDims=2) # the .1 is added by tetgen
    
    elements = pF.loadElementFile('{}\{}.1.ele'.format(path, jobName)) # the .1 is added by tetgen
    
    elementType = 'CPE3' # look at the Pogo docs for more information
    
    
    if modelTimeStep == None:
        safetyFactor = 13.0
        minEdgeLength = pF.minEdgeLength(nodes, elements)
        print 'minimum edge length = {}'.format(minEdgeLength)
        modelTimeStep = minEdgeLength/vL/safetyFactor

    print 'modelTimeStep = {}s'.format(modelTimeStep)
    
    nTimeSteps = int(modelTime/modelTimeStep)
    print 'nTimeSteps = {}'.format(nTimeSteps)
    
    if transducerCentre == None:
        transducerCentre = np.array([xSize/2., ySize])
    
    transducerNodes = pF.findNodesInTransducer(nodes,
                                               transducerCentre,
                                               transducerWidth)
    
    timeBase = np.linspace(0, nTimeSteps-1, nTimeSteps) * modelTimeStep
    
    inputSignal = pF.gaussTone(timeBase, 3.0/f0, f0) #centre of pulse at 2.5/f0
    # but allow a bit of spare room at the start.
    
    signalAmplitude = 100.0
    signalType = 0 # 0 corresponds to a force
    signalDoFAppliedTo = 2 # 2 is the y axis
    
    inSig = [transducerNodes, signalAmplitude, signalDoFAppliedTo, signalType, inputSignal]
    
    history = [transducerNodes, signalDoFAppliedTo]
    
    materialType = 0 # Pogo definition, 0 is an isotropic material
    material1 = [materialType, E, nu, rho]
    
    sPIF.writePogoInputFile(jobName,
                            nodes,
                            elements,
                            (elementType,),
                            (inSig,),
                            (history,),
                            nDims = 2,
                            nt = nTimeSteps,
                            dt = modelTimeStep,
                            materials = (material1,),
                            historyMeasurementFrequency = 1,
                            folderIn = path)
    
    pJS.SubmitRun(jobName, path, cleanup=False, dimensions=2)
    
    result = lPH.loadPogoHistory('{}\{}.pogo-hist'.format(path, jobName))
    
    finalTrace = np.sum(result.histTraces, axis=0)
    time = np.linspace(0, result.nt-1, result.nt)*result.dt
    
    final = np.vstack((time, finalTrace))
    
    return final
    
if __name__ == '__main__':
    path = r'D:\example2DPogoJob'
    plt.close('all')
    trace = pogo2DJob(path)
    sTA.Analyse(trace, trace=True, hilbert=True)
