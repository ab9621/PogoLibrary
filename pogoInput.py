
import numpy as np
import warnings
import subprocess
import pogoFunctions as pF
import pdb
from PolyInterface import poly
class PogoInput:
    def __init__(self,
                 fileName,
                 elementTypes,
                 signals,
                 historyMeasurement,
                 nodes = None,
                 elements = None,
                 geometryFile = None,
                 precision=8,
                 targetMeshSize = 5e-5,
                 nDims=2,
                 nDofPerNode = None,
                 notes = None,
                 runName = 'pogoJob',
                 nt = 100,
                 dt = 1e-8,
                 elementTypeRefs = None,
                 materialTypeRefs = None,
                 orientationRefs = None,
                 elementParameters = None,
                 materials = [[0,7e10,0.34,2700],],
                 orientations = None,
                 boundaryConditions = None,
                 historyMeasurementFrequency = 20,
                 fieldStoreIncrements = None,
                 folderIn = None,
                 totalForce = False,
                 version = 1.03,
                 writeFile = True):
        self.fileName = fileName
        
        ### Header
        self.header = np.array(['']*20, dtype='str')
        if version not in [1.03, 1.04]:
            raise ValueError('Input file version must be 1.03 or 1.04.')
        headerString = '%pogo-inp{}'.format(version)
        for c1 in range(0, len(headerString)):
            self.header[c1] = headerString[c1]
            
        ### Precision
        if precision not in [4,8]:
            raise ValueError('Precision must be 4 or 8.')
        self.precision = np.array([precision,],dtype='int32')
        self.nDims = np.array([nDims,],dtype='int32')
        
        ### Number of degrees of freedom per node
        if nDofPerNode == None:
            nDofPerNode = self.nDims
            
        if nDofPerNode not in [1,2,3]:
            raise ValueError('Number of degrees of freedom must be 1, 2 or 3')
        self.nDofPerNode = np.array([nDofPerNode,],dtype='int32')
        
        ### Set notes
        self.notes = np.array(['']*1024, dtype='str')
        if notes != None:
            if len(notes) > 1024:
                notes = notes[:1024]
                for character in range(len(notes)):
                    self.notes[character] = notes[character]
        
        ### Set runname
        self.runName = np.array(['']*80, dtype='str')
        
        if len(runName) > 80:
            runName = runName[:80]
            
        for character in range(0, len(runName)):
            self.runName[character] = runName[character]
        
        ### Set time step and run time
        self.nt = np.array([nt,],dtype='int32')
        self.dt = np.array([dt,],dtype=self.getPrecString())
        
        ### Node generation if necessary
        if not np.any(nodes) and not geometryFile:
            raise ValueError('Either a poly file or node/element definitions are required')
        elif geometryFile and targetMeshSize and not np.any(elements) and not np.any(nodes):
            if geometryFile.split('.')[-1] == 'dxf':
                print 'Creating poly file from {}'.format(geometryFile)
                poly.poly(geometryFile,elementSize = targetMeshSize,writeFile=True)
            
            if geometryFile.split('.')[-1] == 'poly':
                geometryFile = geometryFile[:-5]
            if self.nDims == 2:
                targetMeshArea = targetMeshSize*targetMeshSize
                subprocess.call('triangle -q -j -a{:.12}F {}.poly'.format(targetMeshArea,geometryFile))
            elif self.nDims == 3:
                targetMeshVolume = targetMeshSize*targetMeshSize*targetMeshSize
                ### Add cwd
                subprocess.call('tetgen {:.12}F {}.poly'.format(targetMeshVolume,geometryFile))
            nodes = pF.loadNodeFile(geometryFile+'.1.node')
            elements = pF.loadElementFile(geometryFile+'.1.ele')
        
        ### Number of nodes and node positions
        if np.shape(nodes)[0] != nDims:
            raise ValueError('nodes must be in shape (nDims, nNodes).')
        self.nNodes = np.array([np.shape(nodes)[1],],dtype = 'int32')
        self.nodes = nodes.astype(self.getPrecString()).T
        
        ### Number of elements and nodes per element
        self.nElements = np.array([np.shape(elements)[1],],dtype='int32')
        self.nNodesPerElement = np.array([np.shape(elements)[0],],dtype='int32')
        
        ### Element type refs
        if elementTypeRefs == None:
            elementTypeRefs = np.zeros(self.nElements)
        if len(elementTypeRefs) != self.nElements:
            raise ValueError('elementTypeRefs must be of length nElements.')
        #if min(elementTypeRefs) != 0:
        #    raise ValueError('elementTypeRefs must be 1 indexed.')
        
        self.elementTypeRefs = elementTypeRefs.astype('int32')# - 1
        
        ### Material type refs
        if materialTypeRefs == None:
            materialTypeRefs = np.zeros(self.nElements)
            
        if len(materialTypeRefs) != self.nElements:
            raise ValueError('materialTypeRefs must be of length nElements.')
            
        #if min(materialTypeRefs) != 1:
        #    raise ValueError('materialTypeRefs must be 1 indexed.')
            
        self.materialTypeRefs = materialTypeRefs.astype('int32') #- 1
        
        ### Element orientations
        if orientationRefs == None:
            orientationRefs = np.zeros(self.nElements,dtype = 'int32')
            
        if len(orientationRefs)!= self.nElements:
            raise ValueError('orientationRefs must be of length nElements.')
            
        if min(elementTypeRefs) < 0: #unused values are set to 0 so -1 in zero indexing
            raise ValueError('orientationRefs must be 1 indexed.')
            
        self.orientationRefs = orientationRefs.astype('int32')# - 1
        
        ### Elements
        if np.max(elements) > self.nNodes:
            raise ValueError('elements points to nodes which are greater than nNodes.')
        
        if np.min(elements) < 0:
            raise ValueError('elements must be 1 indexed.')
            
        self.elements = elements.astype('int32') - 1 #convert to zero indexing
        self.elements = self.elements.T
        
        ### PML sets
        self.nPmlSets = np.array([0,],dtype = 'int32')
        self.pmlParams = np.array([0,],dtype = 'int32')
        
        ### Element types
        self.nElementTypes = np.array([len(elementTypes),],dtype = 'int32')
        if elementParameters == None:
            elementParameters = np.array([0,]*len(elementTypes), dtype = 'int32')
        if np.max(self.elementTypeRefs) > self.nElementTypes - 1:
            raise ValueError('elementTypeRefs points to element types greater than the number of types of element.')    
        self.elementTypes = []
        for ii,elementType in enumerate(elementTypes):
            self.elementTypes.append(ElementType(elementType,elementParameters[ii],self.getPrecString()))
        
        ### Material types
        self.nMaterials = np.array([len(materials),], dtype = 'int32')
        self.materials = []
        for material in materials:
            self.materials.append(Material(material,self.getPrecString()))
        
        ### Orientations
        if orientations == None:
            self.nOr = np.array([0,],dtype ='int32')
            self.orientations = None
        else:
            self.orientations = []
            self.nOr = np.array([len(orientations),],dtype = 'int32')
            for orientation in orientations:
                self.orientations.append(Orientation(orientation,self.getPrecString()))
        
        ### Boundary conditions
        if boundaryConditions == None:
            self.nFixDof =  np.array([0,],dtype ='int32')
            self.boundaryConditions = None
        else:
            nSets = len(boundaryConditions) / 2
            self.nFixDof = np.array([sum([len(boundaryConditions[c1*2]) for c1 in range(nSets)]),],dtype = 'int32')
            self.boundaryConditions = []
            for c1 in range(0,nSets):
                #self.boundaryConditions.append(BoundaryCondition(boundaryConditions[c1]))
                self.boundaryConditions.append(np.array([(boundaryConditions[c1*2]-1)*4 + boundaryConditions[c1*2+1]-1,],dtype='int32'))
        ### Input signals
        self.nInputSignals = np.array([len(signals),],dtype = 'int32')
        self.signals = []
        for signal in signals:
            self.signals.append(Signal(signal,totalForce,self.getPrecString(),dt))
            
        ### History measurements
        if historyMeasurement == None:
            warnings.warn('Warning : No history measurements requested.')
            self.nMeas = 0
            self.historyMeasurement = 0
        else:
            self.nMeas = np.array([len(historyMeasurement),],dtype = 'int32')
            self.historyMeasurement = HistoryMeasurement(historyMeasurement,historyMeasurementFrequency)
            
        ### Field measurements
        if fieldStoreIncrements == None:
            self.nFieldStore = np.array([0,],dtype='int32')
            self.fieldStoreIncrements = np.array([0,],dtype ='int32')
        else:
            self.nFieldStore = np.array([len(fieldStoreIncrements),],dtype = 'int32')
            if np.max(fieldStoreIncrements) > nt or np.min(fieldStoreIncrements) < 1:
                raise ValueError('fieldStoreIncrements out of range [1, nt].')
            self.fieldStoreIncrements = np.array([fieldStoreIncrements-1,],dtype = 'int32')
        
        ### Write to file
        if writeFile:
            self.writeFile()
        
    def getPrecString(self):
        precString = 'float64'
        if self.precision == 4:
            self.precString = 'float32'
        return precString
    
        
    def writeFile(self):
        with open(self.fileName + '.pogo-inp','wb') as f:
            self.header.tofile(f)
            self.precision.tofile(f)
            self.nDims.tofile(f)
            self.nDofPerNode.tofile(f)
            self.notes.tofile(f)
            self.runName.tofile(f)
            self.nt.tofile(f)
            self.dt.tofile(f)
            self.nNodes.tofile(f)
            self.nodes.tofile(f)
            self.nElements.tofile(f)
            self.nNodesPerElement.tofile(f)
            self.elementTypeRefs.tofile(f)
            self.materialTypeRefs.tofile(f)
            self.orientationRefs.tofile(f)
            self.elements.tofile(f)
            self.nPmlSets.tofile(f)
            self.pmlParams.tofile(f)
            self.nElementTypes.tofile(f)
            for elementType in self.elementTypes:
                elementType.writeElementType(f)
            self.nMaterials.tofile(f)
            for material in self.materials:
                material.writeMaterial(f)
            self.nOr.tofile(f)
            if not self.orientations == None:
                for orientation in self.orientations:
                    orientation.writeOrientation(f)
            self.nFixDof.tofile(f)
            if not self.boundaryConditions == None:
                for bc in self.boundaryConditions:
                    bc.tofile(f)
            self.nInputSignals.tofile(f)
            self.signals[0].nt.tofile(f)
            self.signals[0].dt.tofile(f)
            for signal in self.signals:
                signal.writeSignal(f)
            if self.nMeas>0:
                self.historyMeasurement.writeHistory(f)
            else:
                
                np.array([0,], dtype='int32').tofile(f)
                np.array([0,], dtype='int32').tofile(f)
            
            self.nFieldStore.tofile(f)
            self.fieldStoreIncrements.tofile(f)


class Material:
    def __init__(self,materialInfo,precString):
        self.matType = np.array([materialInfo[0],],dtype='int32')
        self.matProps = np.array([materialInfo[1:],],dtype=precString)
        self.nMatParams = np.array([len(materialInfo[1:]),],dtype='int32')
    
    def writeMaterial(self,fileId):
        self.matType.tofile(fileId)
        self.nMatParams.tofile(fileId)
        self.matProps.tofile(fileId)
    
class ElementType:
    def __init__(self,elementType,elementParams,precString):
        self.elTypeSave = np.array(['']*20,dtype='str')
        for character in range(len(elementType)):
            self.elTypeSave[character] = elementType[character]
        
        if elementParams:
            self.nParams = np.array([len(elementParameters),],dtype='int32')
            self.params = np.array(elementParams,dtype = precString)
        else:
            self.params = np.array([0,],dtype='int32')
            self.nParams = np.array([0,],dtype='int32')
            
    def writeElementType(self,fileId):
        self.elTypeSave.tofile(fileId)
        self.nParams.tofile(fileId)
        self.params.tofile(fileId)

class Orientation:
    def __init__(self,orInfo,precString):
        self.paramType = np.array([orInfo[0],], dtype='int32')
        self.nOrParams = np.array([len(orInfo[1:]),],dtype='int32')
        self.paramValues = np.array([orInfo[1:],],dtype = precString)
        
    def writeOrientation(self,fileId):
        self.paramType.tofile(fileId)
        self.nOrParams.tofile(fileId)
        self.paramValues.tofile(fileId)
        
class BoundaryCondition:
    def __init__(self,BCs):
        self.nodes = np.array(BCs[0])
        self.dof = np.array(BCs[1])
    def writeBoundaryCondition(self,fileId):
        dofOut = np.array([(self.nodes-1)*4 + self.dof-1,],dtype='int32')
        dofOut.tofile(fileId)
        
class HistoryMeasurement:
    nodes = np.array([],dtype='int32')
    dofs = np.array([],dtype='int32')
    def __init__(self,histInfo,frequency):
        ###Add Input checking
        
        for history in histInfo:
            self.nodes = np.hstack((self.nodes,history[0]))
            self.dofs = np.hstack((self.dofs,history[1]))
        
        self.frequency = np.array([frequency,],dtype = 'int32')
        self.nMeas = np.array([len(self.nodes),],dtype = 'int32')
        ###Must Add Version 1.04 support
    
    def writeHistory(self,fileId):
        self.nMeas.tofile(fileId)
        self.frequency.tofile(fileId)
        pdb.set_trace()
        outHist = self.nodes*4 + self.dofs - 1
        outHist.tofile(fileId)
class FieldMeasurement:
    def __init__(self,increments=0):
        ###Add input checking
        self.increments = np.array([increments - 1],dtype='int32')
                
class Signal:
    def __init__(self, signalInfo, totalForce, precString,dt):
        if signalInfo:
            nNodes = len(signalInfo[0])
            self.type = np.array([signalInfo[3],],dtype = 'int32')
            # if len(np.unique(signalInfo[0])) != nNodes:
            #     errStr = 'Duplicate nodes cannot be specified for a signal'
            #     raise ValueError(errStr) 
            
            if np.size(signalInfo[1]) != 1 and len(signalInfo[1]) != nNodes:
                raise ValueError('Signal amplitude must be a scalar or a vector of amplitudes for each node signal applied to.')
            
            if signalInfo[3] not in [0,1]:
                raise ValueError('Signal type for signal {} must be 0 or 1.'.format(ii))
            self.nNodes = np.array([len(signalInfo[0]),],dtype='int32')
            self.nodes = np.array(signalInfo[0],dtype = 'int32')
            if type(signalInfo[1]) is float:
                if totalForce == True:
                    if sigType == 1:
                        raise ValueError('totalForce not supported for displacement load.')
                    else:
                        ampVal = signalInfo[1]/nNodes

                else:
                    ampVal = signalInfo[1]
                    
                amp = np.array(np.ones(nNodes)*ampVal, dtype=precString)
            elif type(signalInfo[1]) is np.ndarray:
                if len(signalInfo[1]) != self.nNodes:
                    raise ValueError('If signal amplitude is an array, a value must be specified for each node in the transducer.')
                
                if totalForce == True:
                    raise Warning('totalForce is not supported for loads specified for individual nodes.')

                amp = np.array([signalInfo[1],], dtype=precString)
                    
            else:
                raise ValueError('Signal amplitude not recognised')
            self.amplitude = amp
            self.dof = np.array(signalInfo[2],dtype ='int32')
            
            self.shape = np.array(signalInfo[4],dtype = precString)
            self.dt = np.array(dt,dtype=precString)
            self.nt = np.array(len(signalInfo[4]),dtype = 'int32')
            
    
    def writeSignal(self,fileId):
        self.nNodes.tofile(fileId)
        self.type.tofile(fileId)
        dof = self.nodes*4 + self.dof-1
        dof.tofile(fileId)
        self.amplitude.tofile(fileId)
        self.shape.tofile(fileId)