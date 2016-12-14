# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 10:04:35 2016

@author: Alexander Ballisat

This is a library of functions that have been used for various simulations
performed using pogo. This idea of this script is to bring them all together
into one place to reduce duplication and provide a useful library for
anyone writing pogo input files from Python.

The functions included are:
loadNodeFile
loadElementFile

10/10/2016 Version 1.0
"""
import numpy as np

def loadNodeFile(fileName):
    '''
    Function to load in a .node file as used with tetgen and triangle, store
    the nodes in an array of shape (nDimensions, nNodes) and return it.
    
    Parameters
    ----------
    fileName : string
        The name of the .node file to read. The extension '.node' must be
        passed in fileName.
        
    Returns
    -------
    nodes : array, float
        The array of shape (nDimensions, nNodes) of the nodes loaded in.
    '''
    if fileName[-5:] != '.node':
        raise ValueError('File supplied is not .node file')
    
    nodes = (np.loadtxt(fileName, skiprows=1)[:,1:]).T
    return nodes
    
def loadElementFile(fileName):
    '''
    Function to load in a .ele file as used with tetgen and triangle, store
    the nodes in an array of shape (nNodesPerElement, nElements) and return
    it.
    
    Parameters
    ----------
    fileName : string
        The name of the .ele file to read.
        
    Returns
    -------
    elements : array, float
        The array of shape (nNodesPereElement, nElements) of the nodes loaded
        in.
    '''
    if fileName[-4:] != '.ele':
        raise ValueError('File supplied is not .node file')
    
    elements = (np.loadtxt(fileName, skiprows=1)[:,1:]).T
    return elements.astype('int32')
    
def gauss_tone(t, tau, f0):
    '''
    Function to generate a 5 cycle Gaussian tone burst.
    
    Parameters
    ----------
    t : array, float
        The time vector in which the pulse is to be generated.
        
    tau : float
        The time of the centre of the pulse.
        
    f0 : float
        The centre frequency of the pulse.
        
    Returns
    -------
    gaussian5CycleToneBurst : array, float
        The amplitudes at all times of the tone burst.
    '''
    
    return np.exp(-1.*np.power((t-tau)/(1/(1.0*f0)), 2)) * np.sin(2*np.pi*f0*(t-tau))

def findNodesInTransducer(nodes, transducerCentre, xWidth, yWidth=None, angle=None, shape='rectangular'):
    '''
    Function to find the nodes in a mesh that are on the surface of a
    transducer. It is assumed that the transducer is flat on a surface with
    constant z and has a rectangular profile. In 2D only the x and z axes are
    used, in 3D it is assumed the transducer is fixed on the z axis and has a
    profile in the x-y plane. It will return the inidces of the nodes that 
    are on the transducer.
    
    Parameters
    ----------
    nodes : array, float
        The locations of the nodes in the shape (nDims, nNodes).
        
    transducerCentre : array, float
        The centre of the array in global coordinates. Should be of length
        nDims.
        
    xWidth : float
        The width of the transducer in the x dimension.
        
    yWidth : float, optional
        The width of the transducer in the y dimension. If None, it assumes
        it is square and is set to xWidth
        
    angle : float, optional
        This allows for a rotation of the transducer on the surface of the
        specimen. Angle must be in degrees and is anti-clockwise around the
        positive z axis. 0 degrees points in the direction of the positive
        x axis.
        
    shape : string, optional
        The 2D shape of the transducer footprint. Default is 'rectangular',
        the other allowed options are: 'elliptical' in which case specifying
        a yWidth is needed to specify the other axis, 'circular' in which
        case yWidth is set equal to xWidth. Currently only 'rectangular' is
        supported.
        
    Returns
    -------
    nodeIndices : array, int
        The indices of the nodes which are in the transducer.
    '''
    nDims, nNodes = np.shape(nodes)
    
    if nDims not in [2,3]:
        raise ValueError('nodes must have 2 or 3 rows (coordinates).')
    
    if nDims != len(transducerCentre):
        raise ValueError('transducerCentre must have the same number of dimensions as the node coordinates')
        
    x_min = transducerCentre[0] - xWidth/2.
    x_max = transducerCentre[0] + xWidth/2.
    
    if nDims == 3:
        if yWidth == None:
            yWidth = xWidth
            
        if shape == 'rectangular':
            y_min = transducerCentre[1] - yWidth/2.
            y_max = transducerCentre[1] + yWidth/2.
            
        elif shape == 'circular':
            pass
        
        elif shape == 'elliptical':
            pass
        
        else:
            raise ValueError('shape not recognised. Must be rectanguplar, elliptical \
            or circular.')
            
            
    onSurf = np.where(nodes[-1] == transducerCentre[-1])[0]
    
    if angle == None:
        
        if nDims == 2:
            nodeIndices = np.where((x_min <= nodes[0, onSurf]) & (nodes[0, onSurf] <= x_max))[0]
            
        else:
            nodeIndices = np.where((x_min <= nodes[0, onSurf]) & (nodes[0, onSurf] <= x_max) & (y_min <= nodes[1, onSurf]) & (nodes[1, onSurf] <= y_max))[0]
    
    else:
        if nDims != 3:
            raise ValueError('Rotation only works in 3D')

        A = [x_max, y_max]
        B = [x_max, y_min]
        D = [x_min, y_max]
        
        Ar = rotate2D(A, angle, transducerCentre[:2])
        Br = rotate2D(B, angle, transducerCentre[:2])
        Dr = rotate2D(D, angle, transducerCentre[:2])
        
        AB = Br - Ar
        AD = Dr - Ar
        
        vector1 = xWidth**2
        vector2 = yWidth**2
        
        AM = np.copy(nodes[:2, onSurf]) - np.array([[Ar[0]],[Ar[1]]])
        alpha = AM[0]*AD[0] + AM[1]*AD[1]
        beta = AM[0]*AB[0] + AM[1]*AB[1]
        
        nodeIndices = np.where((0<beta) & (beta<vector2) & (0<alpha) & (alpha<vector1))[0]
        
    return (onSurf[nodeIndices]).astype('int32')
        

def rotate2D(position, angle, centreOfRotation):
    '''
    Function to rotate a point through space about another point through a
    given angle in 2 dimensions.
    
    Parameters
    ----------
    position : array, float
        The point which is to be rotated.
        
    angle : float
        The angle through which the point is to be rotated anti-clockwise, in
        degrees.
        
    centreOfRotation : array, float
        The coordinates about which the point is to be rotated.
        
    Returns
    -------
    coords : float, array
        The coordinates of the rotated point.
    '''
    if len(position) != 2 or len(centreOfRotation) != 2:
        raise ValueError('All coordinates must be in 2D.')
        
    angle *= np.pi/180.

    x = position[0] - centreOfRotation[0]
    y = position[1] - centreOfRotation[1]

    xP = np.cos(angle)*x - np.sin(angle)*y + centreOfRotation[0]
    yP = np.sin(angle)*x + np.cos(angle)*y + centreOfRotation[1]
    
    return np.array([xP, yP])
        
def gaussianAngledBeam(transducerNodes, verticalAngle, dt, nt, velocity, f0, xAngle=None, centre=None):
    '''
    Function to generate the loads needed for each node in a transducer
    surface to create an angled beam in the specimen. The pulse is a 5 cycle
    Gaussian tone burst.
    
    Parameters
    ----------
    transducerNodes : array, float
        The coordinates of the nodes in the transducer in the shape
        (nDims, nNodes).
        
    verticalAngle : float
        The angle of the beam relative to the negative z axis (vertical) in
        degrees.
        
    dt : float
        The time step used in the model in seconds.
        
    nt : int
        The number of time steps in the model.
        
    velocity : float
        The wave velocity of the desired wave type in the material.
        
    f0 : float
        The centre frequency of the pulse.
        
    xAngle : float, optional
        The angle in degrees of the beam in the x-y plane relative to the +ve
        x axis if it is not parallel to the x axis. This can be used if the 
        probe is rotated on the surface. Default is None, the beam is 
        parallel to the x axis.
        
    centre : sequence, float, optional
        The centre of the transducer.
        
    Returns
    -------
    beams : array, float
        The amplitude trace for each node in the transducer surface, of shape
        (nt, nNodes).
    '''
    
    if np.shape(transducerNodes)[0] not in [2,3]:
        raise ValueError('Nodes must have 2 or 3 coordinates.')
        
    nNodes = np.shape(transducerNodes)[1]
    traces = np.zeros((nt, nNodes))
    time = np.linspace(0, nt-1, nt)*dt
    verticalAngle *= np.pi/180.
    #if xAngle != None:
    #    xAngle *= np.pi/180.
    
    if xAngle == None:
        xMin = np.min(transducerNodes[0])
        distances = transducerNodes[0] - xMin
        taus = distances*np.sin(verticalAngle)*1./velocity + 1E-6 #1E-6 is a shift
        
    else:
        if centre == None:
            raise ValueError('Centre must be supplied if xAngle is set.')
        #import matplotlib.pyplot as plt
        centre = centre[:2]
        xPrime, yPrime = rotate2D(transducerNodes[:2], -1.*xAngle, centre)
        #plt.scatter(xPrime, yPrime)
        #plt.gca().set_aspect('equal')
        xMin = np.min(xPrime)
        distances = xPrime - xMin
        taus = distances*np.sin(verticalAngle)*1./velocity + 1E-6 #1E-6 is a shift
    
    for c1 in range(0, nNodes):
        traces[:,c1] = gauss_tone(np.copy(time), taus[c1], f0)
        
    return traces


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    r = np.linspace(0,10,11).astype(int)
    x,y,z = np.meshgrid(r,r,r)
    x = x.flatten()
    y = y.flatten()
    z = z.flatten()
    nodes = np.vstack((x,y,z))
    t = [5,5,10]
    inds = findNodesInTransducer(nodes, t, 2, yWidth=2, angle=45.)
    traces = gaussianAngledBeam(nodes[:,inds], 45., 2E-9, 5000, 3150.,4E6, xAngle=45., centre=t)
    plt.imshow(traces, interpolation='None', aspect='auto')
    for a in inds:
        print nodes[:,a]
