'''
Author: Alexander Ballisat

PLEASE READ THROUGH THIS SECTION AS IT EXPLAINS WHAT IS REQURIED TO RUN THIS 
SCRIPT. The 'path' variable in Step 0 must be set, otherwise everything 
should run if you have it set up correctly.

All units are SI for consistency.

This is an example script to show how a 3D Pogo model can be built, written 
into the Pogo input file, submitted for running and the results accessed.

The Pogo programs needed to run this script are:
- pogoBlock3DGreedy
- pogoSolve3d

A meshing tool is also required. In this script the tetrehedral meshing tool
Tetgen is used. This can be downloaded from:
http://wias-berlin.de/software/tetgen/

For 2D problems, the triangular meshing tool Triangle is regularly used. This
is not needed for this script but can be downloaded from:
https://www.cs.cmu.edu/~quake/triangle.html

Ideally pogoBlock3DGreedy, pogoSolve3d and Tetgen should be placed in a 
folder on the system path and this script assumes that this has been done. To
test if this has been done correctly, open a command prompt and call 
'pogoSolve3d --test' and 'tetgen' (without the apostrophes). The first should
run a simple test to check that Pogo will run on your machine from any 
directory and  the second that Tetgen is running and can be called from any 
directory.

The shapes of the models are written into the .poly format, more information 
about which can be found in the Tetgen manual. A useful tool for visualising 
.poly files and meshes with small numbers of elements is Tetview:
http://wias-berlin.de/software/tetgen/tetview.html

This script is written in Python 2.7.12. If you are new to Python, Chris 
Woods has an excellent tutorial on his website:
http://chryswoods.com/beginning_python/index.html

Note that there are some subtle differences between Python 2 and 3, this 
script will not run on a Python 3 interpreter.

I would recommend using the Anaconda Python distribution as it has all the 
common libraries needed to run this script. It can be downloaded from:
https://www.continuum.io/downloads

The Spyder development environment is a very good IDE, similar in design to 
MATLAB. It is assumed that if you have downloaded this that you are familiar 
with GitHub and using Git in general. If not, TortoiseGit is a good piece of 
software to use:
https://tortoisegit.org/

There are two sets of libraries needed to run this script, common ones that 
are used in many Python scripts and Pogo specific ones. The common ones are:
- scipy (needed for numpy and matplotlib to work but not explicitly imported)
- numpy
- matplotlib
If you have installed Anaconda Python these should be working. The Pogo 
libraries needed are:
- geometryFunctions
- loadPogoHistFile
- pogoFunctons
- savePogoInputFile
- pogoJobSubmission
- singleTraceAnalysis

It is suggested that a folder is created in which you clone the GitHub repo 
PogoLibrary so that you can easily pull updates. In Spyder, this folder can 
be added to the path by going to Tools -> PYTHONPATH manager.

The GPU requires (based on the meshing parameters used in this script) 
approximately 350 MB of RAM for a 10 mm cube (approximately 960K elements),
2500 MB for a 20 mm cube (approximately 7.5M elements) and 8100 MB for a 
30 mm cube (approximately 25.5M elements).

The 10 mm cube should run on most (if not all) modern nVidia graphics cards.

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

# Step 0 - Set the properties of the job
'''
This script is designed to call Pogo directly from Python. In order to 
achieve this, a jobname, which all the files created by this script will 
have, and a path, the directory in which this script is being run, must be 
set, e.g.
path = r'C:\Pogo\Example'
The 'r' delimits the '\' in the string.
A warning on file names: Pogo discards everything after the first full stop,
e.g. 'abc.def' will become 'abc'.
'''
jobName = 'ExamplePogoJob'
path =  # MAKE SURE THIS IS SET

# Test if path is set
if path == '':
    raise ValueError('path must be set.\n')

# Step 1 - Create the geometry
'''
In this case a simple 10 mm cube of aluminium. 
'''
xSize = 0.03
ySize = 0.03
zSize = 0.03

gF.write3DBlockPolyFile(xSize,
                        ySize,
                        zSize,
                        '{}\{}'.format(path, jobName))

# Step 2 - Set some properties of the model 
'''
This is important as the highest frequency of interest and the wavespeed in 
the material determines the mesh size needed. How long to run the model for
and the time step need to be chosen is most easily calculated by
working out the travel time for the furthest reflection on interest. In this
model, the transducer is placed on the top on the block, at z=0.01 and we are
only interested in the first reflection off the backwall.
'''
f0 = 5E6
v = 6400.
modelTime = 4.0 * 2*zSize/v #travel time to back surface and back plus a bit for comfort

'''
The model time step should be approximately less than the travel time between
two nodes. The mesh will use 5 elements per wavelength so can use this in the
calculation. A Gaussian input pulse will be used so the highest frequency is
roughly 2*f0. The time step = v/(2*f0) / 5 / v = 1/(5*2*f0)
'''
modelTimeStep = 1./(5.*2.*f0) # this is just used for the meshing

# Step 3 - mesh the block and load in the nodes and elements
'''
Having of the order of 5 elements per wavelength of the highest  frequency of
interest is necessary to get a good result. The -a switch sets a maximum
volume for any element.
'''
targetMeshSize = modelTimeStep*v #convert the time to a distance

subprocess.call('tetgen -a{:.12f}F {}.poly'.format(targetMeshSize**3, jobName),
                cwd=path)

# With the standard options this should produce approximately 960000 elements

nodes = pF.loadNodeFile('{}\{}.1.node'.format(path, jobName)) # the .1 is added by tetgen

elements = pF.loadElementFile('{}\{}.1.ele'.format(path, jobName)) # the .1 is added by tetgen

elementType = 'C3D4' # look at the Pogo docs for more information

# Step 4 - Set the transducer properties, find the nodes in the transducer
# and generate the input pulse
'''
This needs to include the size of the transducer and the input signal applied
to each node in the footprint of the transducer. In this example, all are 
applied with the same input however applying delays allows the beam to be 
focussed. In this case a square transducer is used. Will also want the
output to be taken at the transducer. The safety factor is to ensure the
stability of the calculation.
'''
safetyFactor = 8.0
modelTimeStep *= 1.0/safetyFactor # reduce the time step to ensure FE calculation is stable
print 'modelTimeStep = {}s'.format(modelTimeStep)

nTimeSteps = int(modelTime/modelTimeStep)
print 'nTimeSteps = {}'.format(nTimeSteps)

transducerWidth = xSize/5.
transducerCentre = np.array([xSize/2., ySize/2., zSize])

transducerNodes = pF.findNodesInTransducer(nodes,
                                           transducerCentre,
                                           transducerWidth)

timeBase = np.linspace(0, nTimeSteps-1, nTimeSteps) * modelTimeStep

inputSignal = pF.gauss_tone(timeBase, 3.0/f0, f0) #centre of pulse at 2.5/f0
# but allow a bit of spare room at the start.

signalAmplitude = 10.0
signalType = 0 # 0 corresponds to a force
signalDoFAppliedTo = 3 # 3 is the z axis

inSig = [transducerNodes, signalAmplitude, signalDoFAppliedTo, signalType, inputSignal]

history = [transducerNodes, signalDoFAppliedTo]

# Step 5 - Set the material properties of the model
'''
These include things like material properties, the length of time of the 
simulation and the time step.
'''
materialType = 0 # Pogo definition, 0 is an isotropic material
E = 7E10
nu = 0.34
rho = 2700.
material1 = [materialType, E, nu, rho]

# Step 6 - Write the Pogo input file
sPIF.writePogoInputFile(jobName,
                        nodes,
                        elements,
                        (elementType,),
                        (inSig,),
                        (history,),
                        nDims = 3,
                        nt = nTimeSteps,
                        dt = modelTimeStep,
                        materials = (material1,),
                        historyMeasurementFrequency=1,
                        folderIn = path)

# Step 6 - Submit the job for blocking and running
'''
The library SubmitRun was written with running large jobs in mind. The 
standard function to call from it SubmitRun has several options which allows 
intermediary files, such as the input file and the block file to be deleted. 
These files can become very large, especially for large models so this can 
save significant hard disk space.
'''
pJS.SubmitRun(jobName, path, cleanup=False, dimensions=3)

# Step 7 - Read in the output of the model
'''
The output file is read into an instance of a class created to store the
outputs of Pogo models. This allows mutliple output files to be loaded and 
compared without having to do to much coding.
'''
result = lPH.loadPogoHistory('{}\{}.pogo-hist'.format(path, jobName))

# Step 8 - Sum the result to a single trace for the whole transducer
# and plot it
finalTrace = np.sum(result.histTraces, axis=0)
time = np.linspace(0, result.nt-1, result.nt)*result.dt

final = np.vstack((time, finalTrace))

sTA.Analyse(final, trace=True, hilbert=True)

'''
This concludes the tutorial script. Contained in the PogoLibrary are many
other functions that may be useful for generating, running and analysing Pogo
models.

If you have many small jobs to run simultaneously, it is worth looking at
the Python Multiprocessing module. A good tutorial can be found here:
http://chryswoods.com/parallel_python/index.html

Modern graphics cards allow Multiple Input Multiple Output (MIMO), thus a
single job can be run per CPU thread. In practice this means that, if the
graphics card you are using has enough RAM, you can run as many jobs in
parallel as you have CPU threads. This does not result in a linear decrease
in the time for each job, i.e. running eight jobs in parallel will not be
eight times faster, but it does provide a significant speed up which makes
it worth while. An example script demonstrating this is currently under
development.
'''
