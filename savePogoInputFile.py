# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 15:53:28 2016

@author: Alexander Ballisat

Script to generate a Pogo input file directly from Python.

14/10/2016 Version 1.0

Please email bugs to alexander.ballisat@bristol.ac.uk
"""

import numpy as np

def writePogoInputFile(fileName,
                       nodePositions,
                       elements,
                       elementTypes,
                       signals,
                       historyMeasurement,
                       precision = 8,
                       nDims = 2,
                       nDofPerNode = None,
                       notes = None,
                       runName = 'pogoJob',
                       nt = 100,
                       dt = 1E-8,
                       elementTypeRefs = None,
                       materialTypeRefs = None,
                       orientationRefs = None,
                       elementParameters = None,
                       materials = [[0, 7E10, 0.34, 2700.],], #Aluminium
                       orientations = None,
                       boundaryConditions = None,
                       historyMeasurementFrequency = 20,
                       historyMeasurementStart = 1,
                       fieldStoreIncrements = None,
                       folderIn = None,
                       totalForce = False,
                       version=1.04):
    '''
    Function to write a Pogo input file directly to the Pogo format. Defaults
    are supplied for all values but most should be set for a model to be run
    correctly. This is the best first thing to check, that the values set are
    correct, if the model does not provide a sensible answer.
    
    Arguments which are by default None calculate themselves in the script.
    
    As is currently the case in Pogo (as of version 1.04), PMLs are not
    supported in this script.
    
    Parameters
    ----------
    fileName : string
        The name of the Pogo input file to be saved.
        
    nodePositions : array, float
        The positions of the nodes, size (nDims, nNodes).
        
    elements : array, int
        The nodes for each element, size (nNodesPerElementMax, nElements).
        1 indexing must be used. Any set to zero are unused.
        
    elementTypes : iterable, string
        Iterable of the element types used, codes must match supported Abaqus
        definitions.
        
    signal : iterable
        Iterable of [nodes signal applied to, signal amplitude, DOF applied
        to, signal type, signal]. Nodes signal applied to must be a sequence,
        Signal amplitude should be a float or of length of the number of 
        nodes the signal is appled to, DOF applied to must be in range 
        [1, nDofPerNode], signal type is either 0 for a force or 1 for a 
        displacement, signal must be of length nt so must be sampled at dt.

    historyMeasurement : iterable
        The nodes and degree of freedom at which to record history data
        specified in (node numbers, degree of freedom). Default is None but
        normally this will be set unless only field measurements are desired.
        Duplicate degrees of freedom for the same node are not allowed.
        
    precision : int, optional
        The precision of the output. Default is 8 (float64) but can be set to
        4 for float32.
        
    nDims : int, optional
        The dumber of dimensions of the model. Must be in the range [1,3], 
        default is 2.
        
    nDofPerNode : float, optional
        The number of degrees of freedom per node. Default is None in which
        case it defautls to the number of dimensions of the model.
        
    notes : string, optional
        Any notes that want to be stored with the model. Default is None but
        up to 1024 characters can be supplied. A cut off is applied to
        anything greater than 1024.
        
    runName : string, optional
        The name of the Pogo job used within Pogo, does not effect any file
        names.
       
    nt : int, optional
        The number of time steps that will be computed. Default is 100.
        
    dt : float, optional
        The time step. Default is 1E-8
        
    elementTypeRefs : array, int, optional
        The references to the different element types used in the model. Must
        be of size (nElements). Default is None which corresponds to a single
        element type which all elements have and is generated in this script.
        
    materialTypeRefs : array, int, optional
        The references to the different material types used in the model. 
        Must be of size (nElements). Default is None which corresponds to a
        single material type which all elements have and is generated in this
        script. Must be 1 indexed.
        
    orientationRefs : array, int, optional
        The references to the orientations for each element. Must be of size
        (nElements). Default is None which corresponds to a single
        orientation type which all elements have and is generated in this
        script.
        
    elementParameters : iterable, optional
        Iterable of the parameters needed for the element types. If none are
        needed set as zero. Default is None which generates a list of zeros
        in the script.
        
    materials : iterable, optional
        List of material properties, each entry should be the material type 
        (0 is isotropic, 1 is anisotropic) followed by a list of material
        properties for that material. for isotropic materials, parameters are
        E, nu, rho, alpha (optional). The order of materials should match 
        that desired in materialTypeRefs. Default is a single material with 
        the properties of Aluminium.
        
    orientations : iterable, optional
        Can specify orientations. Each orientation type should be a list with
        with the parameter type followed by the parameters associated with
        that orientation. Default is None which is orientations not used.
        
    boundaryConditions : iterable, optional
        The boundary conditions applied to the model. The iterable should be
        in the form [nodeSet1, nodeSet1FixedDOF, nodeSet2, nodeSet2FixedDOF,
        etc.]. Each pair should be a set of nodes (1 indexed) and a degree of
        freedom. Default is None in which case no boundary conditions are
        applied.
        
    historyMeasurementFrequency : int, optional
        The number of time increments between history measurements being
        taken. Default is 20.
        
    historyMeasurementStart : int, optional
        The starting increment of the history measurements (1 indexed).
        Default is 1.
        
    fieldStoreIncrements : array, int, optional
        The increments (1 indexed) at which to output the field at. Default
        is None which does not record any field data.
        
    folderIn : string, optional
        Use this if the file is being written in a different location to the
        one that this function is being called from. Default is None.
        
    totalForce : boolean, optional
        This is to be used in conjunction with a force load. Set to True to
        apply a total load which is then divided by the number of nodes that
        a load is applied to (True) or have the supplied force applied to all
        nodes (False). Default is False.
        
    version : float, optional
        The version of the input file, this is used with older versions of
        Pogo. Default is most current, 1.04, but can be changed to 1.03. It
        will fail if anything else is passed.
        
    Returns
    -------
    None
    
    '''
    
    if folderIn == None:
        dFileName = fileName
        
    else:
        dFileName = r'{}\{}'.format(folderIn, fileName)
    
    with open(dFileName + '.pogo-inp', 'wb') as f:

        ##### Generate the file header
        header = np.array(['']*20, dtype='str')
        if version not in [1.03, 1.04]:
            raise ValueError('Input file version must be 1.03 or 1.04.')
        headerString = '%pogo-inp{}'.format(version)
        for c1 in range(0, len(headerString)):
            header[c1] = headerString[c1]
        
        header.tofile(f)
        
        ##### Model precision
        if precision not in [4,8]:
            raise ValueError('Precision must be 4 or 8.')
            
        precString = 'float32'
        if precision == 8:
            precString = 'float64'
            
        precision = np.array([precision,], dtype='int32')
        precision.tofile(f)
        
        ##### Number of dimensions
        if nDims not in [1,2,3]:
            raise ValueError('Number of dimensions must be 1, 2 or 3.')
        nDims = np.array([nDims,], dtype='int32')
        nDims.tofile(f)
        
        ##### Number of degrees of freedom per node
        if nDofPerNode == None:
            nDofPerNode = nDims
            
        if nDofPerNode not in [1,2,3]:
            raise ValueError('Number of degrees of freedom must be 1, 2 or 3')
            
        nDofPerNode = np.array([nDofPerNode,], dtype='int32')
        nDofPerNode.tofile(f)
        
        ##### Notes
        notesA = np.array(['']*1024, dtype='str')
        
        if notes != None:
            if len(notes) > 1024:
                notes = notes[:1024]
            for c1 in range(0, len(notes)):
                notesA[c1] = notes[c1]
        
        notesA.tofile(f)
        
        ##### Run name
        runNameA = np.array(['']*80, dtype='str')
        
        if len(runName) > 80:
            runName = runName[:80]
            
        for c1 in range(0, len(runName)):
            runNameA[c1] = runName[c1]
            
        runNameA.tofile(f)
        
        ##### Number of time steps
        nt = np.array([nt,], dtype='int32')
        nt.tofile(f)
        
        ##### Time step
        dt = np.array([dt,], dtype=precString)
        dt.tofile(f)
        
        ##### Number of nodes
        nNodes = np.shape(nodePositions)[1]
        if np.shape(nodePositions)[0] != nDims:
            raise ValueError('nodePositions must be in shape (nDims, nNodes).')
        nNodesA = np.array([nNodes,], dtype='int32')
        nNodesA.tofile(f)
        
        ##### Node positions
        nodePositions = nodePositions.astype(precString)
        nodePositions = nodePositions.T
        nodePositions.tofile(f)
        
        ##### Number of elements
        nElements = np.shape(elements)[1]
        nElementsA = np.array([nElements,], dtype='int32')
        nElementsA.tofile(f)
        
        ##### Nodes per element
        nNodesPerElement = np.shape(elements)[0]
        nNodesPerElementA = np.array([nNodesPerElement,], dtype='int32')
        nNodesPerElementA.tofile(f)
        
        ##### Element type references
        if elementTypeRefs == None:
            elementTypeRefs = np.ones(nElements)
            
        if len(elementTypeRefs) != nElements:
            raise ValueError('elementTypeRefs must be of length nElements.')
            
        if min(elementTypeRefs) != 1:
            raise ValueError('elementTypeRefs must be 1 indexed.')
            
        elementTypeRefs = elementTypeRefs.astype('int32') - 1
        elementTypeRefs.tofile(f)
        
        ##### Material type refs
        if materialTypeRefs == None:
            materialTypeRefs = np.ones(nElements)
            
        if len(materialTypeRefs) != nElements:
            raise ValueError('materialTypeRefs must be of length nElements.')
            
        if min(materialTypeRefs) != 1:
            raise ValueError('materialTypeRefs must be 1 indexed.')
            
        materialTypeRefs = elementTypeRefs.astype('int32') #- 1
        materialTypeRefs.tofile(f)
        
        ##### Element orientations
        if orientationRefs == None:
            orientationRefs = np.zeros(nElements)
            
        if len(orientationRefs)!= nElements:
            raise ValueError('orientationRefs must be of length nElements.')
            
        if min(elementTypeRefs) < 0: #unused values are set to 0 so -1 in zero indexing
            raise ValueError('orientationRefs must be 1 indexed.')
            
        orientationRefs = orientationRefs.astype('int32') - 1
        orientationRefs.tofile(f)
        
        ##### Elements
        if np.max(elements) > nNodes:
            raise ValueError('elements points to nodes which are greater than nNodes.')
        
        if np.min(elements) < 0:
            raise ValueError('elements must be 1 indexed.')
            
        elements = elements.astype('int32') - 1 #convert to zero indexing
        elements = elements.T
        elements.tofile(f)
        
        ##### PML sets
        zero = np.array([0,], dtype='int32')
        zero.tofile(f) #for nPML sets
        zero.tofile(f) #for nPML parameters
        
        ##### Element types
        nElementTypes = len(elementTypes)
        if elementParameters == None:
            elementParameters = np.array([0,]*nElementTypes, dtype='int32')
        
        if np.max(elementTypeRefs) > nElementTypes - 1:
            raise ValueError('elementTypeRefs points to element types greater than the number of types of element.')
        
        nElementTypesA = np.array([nElementTypes,], dtype='int32')
        nElementTypesA.tofile(f)
        
        for c1 in range(0, nElementTypes):
            if elementTypes[c1] == '':
                raise ValueError('Element code for element type {} not defined'.format(c1+1))
        
            elTypeSave = np.array(['']*20, dtype='str')
            for c2 in range(0, len(elementTypes[c1])):
                elTypeSave[c2] = elementTypes[c1][c2]
            elTypeSave.tofile(f)
            
            if elementParameters[c1] == 0:
                zero.tofile(f)
                zero.tofile(f)
            
            else:
                nParams = len(elementParameters[c1])
                nParams = np.array([nParams,], dtype='int32')
                params = np.array(elementParameters[c1], dtype=precString)
                nParams.tofile(f)
                params.tofile(f)
                
        ##### Material types
        nMaterials = len(materials)
        nMaterialsA = np.array([nMaterials,], dtype='int32')
        nMaterialsA.tofile(f)
        
        for c1 in range(0, nMaterials):
            matType = np.array([materials[c1][0],], dtype='int32')
            matType.tofile(f)
            matProps = np.array(materials[c1][1:], dtype=precString)
            nMatParams = np.array([len(matProps),], dtype='int32')
            
            nMatParams.tofile(f)
            matProps.tofile(f)
            
        ##### Orientations
        if orientations == None:
            nOr = 0
            zero.tofile(f)
            
        else:
            nOr = len(orientations)
            
            for c1 in range(0, nOr):
                paramType = np.array([orientations[c1][0],], dtype='int32')
                paramType.tofile(f)
                nOrParams = np.array([len(orientations[c1][1:]),], dtype='int32')
                nOrParams.tofile(f)
                paramValues = np.array([orientations[c1][1:],], dtype=precString)
                paramValues.tofile(f)
                
        ##### Boundary Conditions
        if boundaryConditions == None:
            nFixDof = 0
            
        else:
            nFixDof = len(boundaryConditions)/2
            
        if nFixDof == 0:
            zero.tofile(f)
        
        else:
            nFixDofA = np.array([nFixDof,], dtype='int32')
            nFixDofA.tofile(f)
            
            for c1 in range(0, nFixDof):
                dof = (boundaryConditions[c1*2]-1)*4 + boundaryConditions[c1*2+1]-1
                dof = np.array([dof,], dtype='int32')
                
        ##### Input signals
        nInputSignals = len(signals)
        nInputSignalsA = np.array([nInputSignals,], dtype='int32')
        nInputSignalsA.tofile(f)
        
        nt.tofile(f)
        dt.tofile(f)
        
        for c1 in range(0, nInputSignals):
            if len(signals[c1]) != 5:
                raise ValueError('signal {} is not properly defined, only {} parts set.'.format(c1, len(signals[c1])))
                
            nNodes = len(signals[c1][0])
            nNodesA = np.array([nNodes,], dtype='int32')
            
            if len(np.unique(signals[c1][0])) != nNodes:
                raise ValueError('Duplicate nodes cannot be specified for a signal.')
                
            nNodesA.tofile(f)
            
            if signals[c1][3] not in [0,1]:
                raise ValueError('Signal type for signal {} must be 0 or 1.'.format(c1))
                
            sigType = np.array([signals[c1][3],], dtype='int32')
            sigType.tofile(f)
            
            if signals[c1][2] < 1 or signals[c1][2] > nDofPerNode:
                raise ValueError('DoF for signal {} not valid.'.format(c1))
                
            dof = (np.array(signals[c1][0]))*4 + signals[c1][2]-1
            dof = np.array([dof,], dtype='int32')
            dof.tofile(f)
            
            if len(signals[c1][4]) != nt:
                raise ValueError('Signal {} must be length nt'.format(c1))
            
            if np.size(signals[c1][1]) != 1 and len(signals[c1][1]) != nNodes:
                raise ValueError('signal {} amplitude must be a scalar or a vector of amplitudes for each node signal applied to.')
            
            ##### this needs fixing - need to generate an array if only a float is passed.
            if type(signals[c1][1]) is float or type(signals[c1][1]) is np.float64:
                if totalForce == True:
                    if sigType == 1:
                        raise ValueError('totalForce not supported for displacement load.')
                    else:
                        ampVal = signals[c1][1]/nNodes

                else:
                    ampVal = signals[c1][1]
                    
                amp = np.array(np.ones(nNodes)*ampVal, dtype=precString)

            elif type(signals[c1][1]) is np.ndarray:
                if len(signals[c1][1]) != nNodes:
                    raise ValueError('If signal amplitude is an array, a value must be specified for each node in the transducer.')
                
                if totalForce == True:
                    raise Warning('totalForce is not supported for loads specified for individual nodes.')

                amp = np.array([signals[c1][1],], dtype=precString)
                    
            else:
                raise ValueError('Signal amplitude not recognised')
                
            amp.tofile(f)
            
            sig = np.array(signals[c1][4], dtype=precString)
            sig.tofile(f)
            

            
        ##### History measurements
        if historyMeasurement == None:
            print '\nWarning : No history measurements requested.\n'
            zero.tofile(f)
            zero.tofile(f)
            zero.tofile(f)
        
        else:
            hist = np.array([])
            for c1 in range(0, len(historyMeasurement)):
                r = (np.array(historyMeasurement[c1][0]))*4 + historyMeasurement[c1][1]-1
                hist = np.hstack((hist,r))
                
            if len(np.unique(hist)) != len(hist):
                raise ValueError('Duplicate degrees of freedom for the same node(s) found.')
                
            nMeas = np.array([len(hist),], dtype='int32')
            nMeas.tofile(f)
            
            historyMeasurementFrequency = np.array([historyMeasurementFrequency,], dtype='int32')
            historyMeasurementFrequency.tofile(f)
            
            if version == 1.04:
                historyMeasurementStart = np.array([historyMeasurementStart-1,], dtype='int32')
                historyMeasurementStart.tofile(f)
            
            hist = hist.astype('int32')
            hist.tofile(f)
            
        ##### Field measurements
        if fieldStoreIncrements == None:
            zero.tofile(f)
            
        else:
            nFieldStore = np.array([len(fieldStoreIncrements),], dtype='int32')
            nFieldStore.tofile(f)
            
            if np.max(fieldStoreIncrements) > nt or np.min(fieldStoreIncrements) < 1:
                raise ValueError('fieldStoreIncrements out of range [1, nt].')
                
            fieldStoreIncrements = np.array([fieldStoreIncrements,], dtype='int32')
            fieldStoreIncrements.tofile(f)
        
        print '\nPogo input file written.\n'
        
        return
        
        
#if __name__ == '__main__':
#    name = 'TestPogoInputFile'
#    w = np.linspace(0, 50, 101)
#    nodes = np.vstack((w,w))
#    p = np.linspace(0,100,101)
#    elements = np.vstack((np.roll(np.copy(p), 1),np.roll(np.copy(p), 2), np.roll(np.copy(p), 3), np.roll(np.copy(p), 4)))
#    elTypes = ['C3D4',]
#    signal = [[(3,), 11.5, 2, 0, np.linspace(0,1.0,100)],]
#    hist = [[[1,2,3,4,5,6], 2],]
#    
#    writePogoInputFile(name, nodes, elements, elTypes, signal, historyMeasurement=hist)