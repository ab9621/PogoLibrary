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
gaussTone
findNodesInTransducer
rotate2D
gaussianAngledBeam
tukeyWindow
cartesianCombinations
gaussianBeamProfile


10/10/2016 Version 1.0
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
import itertools

def __criticalCalc__(a,b):
    '''
    Function used in the Snell's Law calculations
    '''
    tmp = a*b
    if tmp > 1.0:
        return np.inf
    else:
        return np.arcsin(tmp)
        
def loadNodeFile(fileName, nDims=3, regionAttributes=False):
    '''
    Function to load in a .node file as used with tetgen and triangle, store
    the nodes in an array of shape (nDimensions, nNodes) and return it.
    
    Parameters
    ----------
    fileName : string
        The name of the .node file to read. The extension '.node' must be
        passed in fileName.
        
    regionAttributes : boolean, optional
        Whether region attributes were specified in the generation of the 
        mesh.
        
    Returns
    -------
    nodes : array, float
        The array of shape (nDimensions, nNodes) of the nodes loaded in.
    '''
    if fileName[-5:] != '.node':
        raise ValueError('File supplied is not .node file')
    
    if nDims not in [2,3]:
        raise ValueError('The nodes must have 2 or 3 dimensions, not {}'.format(nDims))
    nodes = (np.genfromtxt(fileName, dtype='float64', skip_header=1)).T
    
    if regionAttributes == False:
        if nDims == 3:
            nodes = nodes[1:]
        elif nDims == 2:
            nodes = nodes[1:-1]
    else:
        if nDims == 3:
            nodes = nodes[1:-2]
        elif nDims == 2:
            nodes = nodes[1:-3]
    return nodes
    
def loadElementFile(fileName, regionAttributes=False):
    '''
    Function to load in a .ele file as used with tetgen and triangle, store
    the nodes in an array of shape (nNodesPerElement, nElements) and return
    it.
    
    Parameters
    ----------
    fileName : string
        The name of the .ele file to read.
        
    regionAttributes : boolean, optional
        Whether region attributes were specified in the generation of the 
        mesh.
        
    Returns
    -------
    elements : array, float
        The array of shape (nNodesPereElement, nElements) of the nodes loaded
        in.
    '''
    if fileName[-4:] != '.ele':
        raise ValueError('File supplied is not .node file')
    
    elements = (np.genfromtxt(fileName, dtype='int32', skip_header=1)).T
    
    if regionAttributes == False:
        elements = elements[1:]
        return elements.astype('int32')

    else:
        regionAttributes = np.copy(elements[-1])
        elements = elements[1:-1]
        return elements.astype('int32'), regionAttributes.astype('int32')
    
def gaussTone(t, tau, f0):
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

def findNodesInTransducer(nodes, 
                          transducerCentre, 
                          xWidth, 
                          yWidth=None, 
                          angle=None, 
                          shape='rectangular'):
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
        
def gaussianAngledBeam(transducerNodes, 
                       verticalAngle, 
                       dt, 
                       nt, 
                       velocity, 
                       f0, 
                       xAngle=None,
                       centre=None):
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
    
    if xAngle == None:
        xMin = np.min(transducerNodes[0])
        distances = transducerNodes[0] - xMin
        taus = distances*np.sin(verticalAngle)*1./velocity + 1E-6 #1E-6 is a shift
        
    else:
        if centre == None:
            raise ValueError('Centre must be supplied if xAngle is set.')
        centre = centre[:2]
        xPrime, yPrime = rotate2D(transducerNodes[:2], -1.*xAngle, centre)
        xMin = np.min(xPrime)
        distances = xPrime - xMin
        taus = distances*np.sin(verticalAngle)*1./velocity + 1E-6 #1E-6 is a shift
    
    for c1 in range(0, nNodes):
        traces[:,c1] = gaussTone(np.copy(time), taus[c1], f0)
        
    return traces

def tukeyWindow(left, right, taperWidth, nWindow, rightTaperWidth=None):
    '''
    Function to return a Tukey window.
    
    Parameters
    ----------
    left : int
        The index of the position where the window first reaches 1
        
    right : int
        The index of the position where the window last reaches 1
        
    taperWidth : int
        The width of the taper regions
        
    nWindow : int
        The length of the whole array in which the window sits
        
    rightTaperWidth : optional
        If this is None then the taper on both sides are of equal length. If
        set to an int, this is the width of the right taper and taperWidth is
        the width of the left taper. Default is None in which case the window
        is symmetric
        
    Returns
    -------
    window : array, float
        The tukey window
    '''
    lWidth = taperWidth
    if rightTaperWidth == None:
        rWidth = taperWidth
        
    else:
        try:
            rWidth = int(rightTaperWidth)
        except:
            raise ValueError('rightTaperWidth must be an int.')
            
    window = np.zeros(nWindow)
    window[left:right] = 1.0
    window[left-lWidth:left] = 0.5*(1.-np.cos(np.pi*np.linspace(0, lWidth-1, lWidth)/(lWidth-1)))
    window[right:right+rWidth] = 0.5*(1+np.cos(np.pi*(np.linspace(0, rWidth-1, rWidth))/(rWidth-1)))
    
    return window
    
def cartesianCombinations(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """
    
    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesianCombinations(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

def gaussianBeamProfile(nodes, 
                        transducerCentre, 
                        sigmaX, 
                        sigmaY,
                        sigmaZ=None,
                        amplitude=1,
                        xAngle=None,
                        plotting=False):
    '''
    Function to calculate the amplitudes needed to generate a Gaussian beam
    profile across the foot print of a transducer. It can calculate both 1
    and 2 dimensional Gaussian functions.
    
    Parameters
    -----------
    nodes : array, float
        The coordinates of the nodes in the transducer passed in an array of
        shape (nDims, nNodes). It must have 3 >= nDims >= 2 for a 1 
        dimensional Gaussian profile and nDims = 3 for a 2 dimensional
        profile.
        
    transducerCentre : float
        The coordinates of the centre of the transducer. This acts as the
        mean of the distribution.
        
    sigmaX : float
        The standard deviation of the profile in the x-dimension, defined as
        the first coordinate in nodes. Can be None if sigmaY is specified in
        which case it is a 1 dimension distribution in the y axis (second 
        coordinate) or sigmaZ is specified in which case it is a 1 dimension
        distribution in the z axis (third coordinate).
        
    sigmaY : float
        The standard deviation of the profile in the y-dimension, defined as
        the second coordinate in nodes. See documentation for sigmaX for more
        details.
        
    sigmaZ : float, optional
        The standard deviation of the profile in the z-dimension, defined as
        the thrid coordinate in nodes. See documentation for sigmaX for more
        details. Default is None in which case it is not used.
    
    amplitude : float
        The maximum amplitude of the Gaussian function.
        
    xAngle : float
        This allows for a rotation of the transducer on the surface of the
        specimen. Angle must be in degrees and is anti-clockwise around the
        positive z axis. 0 degrees points in the direction of the positive
        x axis.
        
    plotting : boolean
        Whether to plot the generated distribution or not. This is a designed
        as a convenient way of quickly visualising the beam profile. It works
        for 1 dimension and 2 dimension Gaussian functions.
    
    Returns
    --------
    amplitudes : array
        The amplitudes needed for each node in nodes, in the same order, to
        create the Gaussian profile.
    '''
    if np.shape(nodes)[0] not in [2,3]:
        raise ValueError('nodes must have 2 or 3 coordinates')
    
    if np.shape(nodes)[0] != len(transducerCentre):
        raise ValueError('transducerCentre must have the same number of \
        dimensions as nodes.')
    if sigmaX == None and sigmaY == None and sigmaZ == None:
        raise ValueError('One or two of sigmaX, sigmaY and sigmaZ must be \
        specified')
        
    if sigmaX != None and sigmaY != None and sigmaZ != None:
        raise ValueError('A maximum of two of sigmaX, sigmaY and sigmaZ can \
        be specified')

    nNodes = len(nodes[0])
        
    vals = [sigmaX, sigmaY, sigmaZ]
    axes = [[pos,val] for pos, val in enumerate(vals) if val != None]

    amplitudes = np.zeros(nNodes)
    
    if len(axes) == 1:
        amplitudes = amplitude*np.exp(-1*np.power((nodes[axes[0][0]]-transducerCentre[axes[0][0]])/axes[0][1], 2))
        
    elif len(axes) == 2:
        if xAngle != None:
            if len(nodes) != 3:
                raise ValueError('Rotation only works in 3D')
                
            centre = transducerCentre[:2]
            xPrime, yPrime = rotate2D(nodes[:2], -1.*xAngle, centre)
            nodes[0] = xPrime
            nodes[1] = yPrime

        amplitudes = amplitude*np.exp(-1*(np.power((nodes[axes[0][0]]-transducerCentre[axes[0][0]])/axes[0][1], 2) \
                    + np.power((nodes[axes[1][0]]-transducerCentre[axes[1][0]])/axes[1][1], 2)))
    
    else:
        raise ValueError('Three axis calculations are not yet implemented')
    
    if plotting == True:
        plt.figure()
        if len(axes) == 1:
            plt.plot(nodes[axes[0][0]], amplitudes, 'o', ms=4)
            plt.xlabel('{} coordinate'.format(axes[0][0]))
            plt.ylabel('Amplitude (Arbitrary Units)')
            
        elif len(axes) == 2:
            plt.scatter(nodes[axes[0][0]], nodes[axes[1][0]], c=amplitudes)
            plt.xlabel('{} coordinate'.format(axes[0][0]))
            plt.ylabel('{} coordinate'.format(axes[1][0]))
            plt.colorbar(label='Amplitude (Arbitrary Units)')
    
    return amplitudes
    
def hanningBeamProfile(nodes, 
                        transducerCentre,
                        transXWidth,
                        transYWidth,
                        transZWidth = None,
                        amplitude=1,
                        plotting=False):
    '''
    Function to calculate the amplitudes needed to generate a Gaussian beam
    profile across the foot print of a transducer. It can calculate both 1
    and 2 dimensional Gaussian functions.
    
    Parameters
    -----------
    nodes : array, float
        The coordinates of the nodes in the transducer passed in an array of
        shape (nDims, nNodes). It must have 3 >= nDims >= 2 for a 1 
        dimensional Gaussian profile and nDims = 3 for a 2 dimensional
        profile.
        
    transducerCentre : float
        The coordinates of the centre of the transducer.
        
    transXWidth : float
        The width of the probe in the x-dimension, defined as the first
        coordinate in nodes. Can be None if transYWidth is specified in
        which case it is a 1 dimension distribution in the y axis (second 
        coordinate) or transZWidth is specified in which case it is a 1 
        dimension distribution in the z axis (third coordinate). If two
        widths are supplied then the profile is the outer product of two
        hanning windows.
        
    transYWidth : float
        The width of the transducer in the y-dimension, defined as
        the second coordinate in nodes. See documentation for transXWidth for
        more details.
        
    transZWidth : float, optional
        The width of the transducer in the z-dimension, defined as
        the thrid coordinate in nodes. See documentation for transXWidth for 
        more details. Default is None in which case it is not used.
    
    amplitude : float
        The maximum amplitude of the Gaussian function.
        
    plotting : boolean
        Whether to plot the generated distribution or not. This is a designed
        as a convenient way of quickly visualising the beam profile. It works
        for 1 dimension and 2 dimension Gaussian functions.
    
    Returns
    --------
    amplitudes : array
        The amplitudes needed for each node in nodes, in the same order, to
        create the Gaussian profile.
    '''
    if np.shape(nodes)[0] not in [2,3]:
        raise ValueError('nodes must have 2 or 3 coordinates')
    
    if np.shape(nodes)[0] != len(transducerCentre):
        raise ValueError('transducerCentre must have the same number of \
        dimensions as nodes.')

    nNodes = len(nodes[0])
    
    vals = [transXWidth, transYWidth, transZWidth]
    axes = [[pos,val] for pos, val in enumerate(vals) if val != None]

    amplitudes = np.zeros(nNodes)
    
    if len(axes) == 1:
        nPointsWindow = 101
        h1 = np.hanning(nPointsWindow)
        base = np.linspace(transducerCentre[axes[0][0]] - axes[0][1]/2.,
                           transducerCentre[axes[0][0]] + axes[0][1]/2., 
                           nPointsWindow)
        interp = si.interp1d(base, h1, kind='cubic')
        try:
            amplitudes = amplitude*interp(nodes[axes[0][0]])
        except:
            raise ValueError('Error doing interpolation. Possible that some \
            nodes are outsie convex hull of transducer footprint.')
        amplitudes[amplitudes<0.0] = 0.0
        
    elif len(axes) == 2:
        nPointsWindow = 101
        h1 = np.hanning(nPointsWindow)
        h2 = np.hanning(nPointsWindow)
        base1 = np.linspace(transducerCentre[axes[0][0]] - axes[0][1]/2.,
                           transducerCentre[axes[0][0]] + axes[0][1]/2., 
                           nPointsWindow)
        base2 = np.linspace(transducerCentre[axes[1][0]] - axes[1][1]/2.,
                           transducerCentre[axes[1][0]] + axes[1][1]/2., 
                           nPointsWindow)
                           
        interp1 = si.interp1d(base1, h1, kind='cubic')
        interp2 = si.interp1d(base2, h2, kind='cubic')
        
        a1 = interp1(nodes[axes[0][0]])
        a1[a1<0.0] = 0.0
        a2 = interp2(nodes[axes[1][0]])
        a2[a2<0.0] = 0.0
        amplitudes = amplitude*a1*a2
        
    else:
        raise ValueError('Three axis calculations are not yet implemented')
    if plotting == True:
        plt.figure()
        if len(axes) == 1:
            plt.plot(nodes[axes[0][0]], amplitudes, 'o', ms=4)
            plt.xlabel('{} coordinate'.format(axes[0][0]))
            plt.ylabel('Amplitude (Arbitrary Units)')
            
        elif len(axes) == 2:
            plt.scatter(nodes[axes[0][0]], nodes[axes[1][0]], c=amplitudes)
            plt.xlabel('{} coordinate'.format(axes[0][0]))
            plt.ylabel('{} coordinate'.format(axes[1][0]))
            plt.colorbar(label='Amplitude (Arbitrary Units)')
    
    return amplitudes
    
def waveVelocity(E, nu, rho):
    '''
    Convenience function to calculate the longitudinal and shear velocities
    given some material properties. All units are SI.
    
    Parameters
    ----------
    E : float
        The Young's modulus of the material.
        
    nu : float
        The Poisson's ratio of the material.
        
    rho : float
        The density of the material.
        
    Returns
    -------
    cp : float
        The longitudinal wave velocity.
        
    cs : float
        The shear wave velocity.
    '''
    cp = np.sqrt((E*(1.-nu))/(rho*(1.+nu)*(1-2.*nu)))
    cs = np.sqrt(E/(2.*(1+nu)*rho))
    return cp, cs
   
def minEdgeLength(nodes, elements):
    '''
    Function to find the minimum edge length of any element in a mesh.
    
    Parameters
    ----------
    nodes : array
        The position of the nodes, either in 2 or 3 dimensions.
        
    elements : array
        The nodes which make up each element.
        
    Returns
    -------
    minDistance : float
        The minimum length of any element side of any element.
    '''
    minDistance = 1000000000.0 #pointlessly large number
    
    nNodesPerElement = len(elements)
    print '{} nodes per element'.format(nNodesPerElement)
    
    combinations = itertools.combinations(range(1, nNodesPerElement+1), 2)
    
    for a in combinations:
        p1 = nodes[:, elements[a[0]-1, :]-1 ]
        p2 = nodes[:, elements[a[0]-2, :]-1 ]
        
        tmpMin = np.min(np.sum((p1-p2)**2, axis=0))
        
        if tmpMin < minDistance:
            minDistance = tmpMin
    
    minDistance = np.sqrt(minDistance)        
    return minDistance
    
def snellsLawModeConversion(vL1,vL2,vS1,vS2,theta,mode=0):
    '''
    Function to calculate the angles of reflection and refraction for a wave
    impinging on a boundary.
    
    Parameters
    ----------
    vL1 : float
        The longitudinal wave velocity in the material in which the wave is
        incident.
        
    vL2 : float
        The longitudinal wave velocity in the material in which the wave is
        travelling into.
        
    vS1 : float
        The shear wave velocity in the material in which the wave is 
        incident.
        
    vS2 : float
        The shear wave velocity in the material in which the wave is
        travelling into.
        
    theta : float
        The incident angle of the wave in degrees.
        
    mode : int
        Whether the incidnet wave is longitudinal (0) or shear (1). Default
        is 0, a longitudinal wave.
        
    Returns
    -------
    reflectedLong : float
        The angle of the reflected longitudinal wave. None indicates that
        the incident angle is beyond the critical angle therefore no wave
        of this type is generated. Angle is is degrees.
        
    reflectedShear : float
        The angle of the reflected shear wave. None indicates that the 
        incident angle is beyond the critical angle therefore no wave of this
        type is generated. Angle is is degrees.
    
    refractedLong : float
        The angle of the refracted longitudinal wave. None indicates that
        the incident angle is beyond the critical angle therefore no wave
        of this type is generated. Angle is is degrees.
    
    refractedShear : float
        The angle of the refracted shear wave. None indicates that the 
        incident angle is beyond the critical angle therefore no wave of this
        type is generated. Angle is is degrees.        
    
    '''
    if mode not in [0,1]:
        raise ValueError('mode must be 0 or 1, not {}'.format(mode))
    
    theta = np.deg2rad(theta)
    
    if mode == 0:
        v0 = vL1
    elif mode == 1:
        v0 = vS1
        
    prod = np.sin(theta)/v0
    
    refractedLong = np.rad2deg(__criticalCalc__(prod, vL2))
    refractedShear = np.rad2deg(__criticalCalc__(prod, vS2))
    
    if mode == 0:
        reflectedLong = np.rad2deg(theta)
        reflectedShear = np.rad2deg(__criticalCalc__(prod, vS1))
        
    if mode == 1:
        reflectedLong = np.rad2deg(__criticalCalc__(prod, vL1))
        reflectedShear = np.rad2deg(theta)
        
    return reflectedLong, reflectedShear, refractedLong, refractedShear
    
def plotFieldData(fieldData, increment, component='magnitude',
                  returnFig=False):
    '''
    Function to plot the field data for a given increment.
    
    Parameters
    ----------
    fieldData : instance of __fieldObject__ class. 
        The data to be plotted contained in an instance of the 
        __fieldObject__ class. This is the class generated by the function in
        this library to load in a field history file.
        
    increment : int
        The increment at which the data is to be plotted. Must be 1 indexed.
        
    component : string, optional
        The component of the displacement to be plotted. Default is
        'magnitude' in which case the magnitude of the displacement vector is
        plotted. Other valid values are 'ux', 'uy' or 'uz'.
        
    returnFig : boolean, optional
        Whether or not to return the generated figure. Default is False in
        which case it is not returned.
        
    Returns
    -------
    fig : matplotlib figure, optional
        The figure instance of the plot.
    '''
    if component not in ['ux', 'uy', 'uz', 'magnitude']:
        raise ValueError('Invalid displacement component.')
        
    if increment<1 and increment>fieldData.nFieldIncs:
        raise ValueError('Invalid increment, must be in range [1,{}].'.format(fieldData.nFieldInc))
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    if component == 'magnitude':
        data = np.sqrt(fieldData.ux[:,increment]**2 + fieldData.uy[:,increment]**2 + fieldData.uz[:,increment]**2)
        
    elif component == 'ux':
        data = fieldData.ux[:,increment]
        
    elif component == 'uy':
        data = fieldData.uy[:,increment]
        
    elif component == 'uz':
        data = fieldData.uz[:,increment]
    
    if fieldData.nDims == 2:
        #try an interpolator - this works
        print 'Doing interpolation for plotting'
        inputCoords = np.vstack((fieldData.nodePos[0], fieldData.nodePos[1])).T
        interp = si.LinearNDInterpolator(inputCoords, data, fill_value=0.0)
        nx = 400
        ny = 340
        xBase = np.linspace(np.min(fieldData.nodePos[0]), np.max(fieldData.nodePos[0]), nx)
        yBase = np.linspace(np.min(fieldData.nodePos[1]), np.max(fieldData.nodePos[1]), ny)
        
        xs, ys = np.meshgrid(xBase, yBase)
        
        finalData = interp(xs, ys)
        
        print 'Plotting'
        extent_ = [xBase[0], xBase[-1], yBase[0], yBase[-1]]
        ax.imshow(finalData, origin='lower', extent=extent_, aspect='auto')
        
    elif fieldData.nDims == 3:
        print 'Not implemented yet...'
        
    if returnFig == True:
        return fig
        
    else:
        return
    
def createRectOrientation(phis,origin = np.array([[0,0,0],]),addRotDim=3,addRotAngle=0):
    '''
    Function to create list of orientations based on angles relative to the x-axis. Currently 2D only
    
    Parameters
    ----------
    phis: iterable containing desired rotation angles
    origin: origin of the coordinate system
    addRotDim: dimension of any additional rotation (default is 3; z-axis)
    addRotAngle: scalar defining additional rotation about axis defined by addRotDim
    NOTE: the additional rotation is a hangover from more complex rotation types (e.g. cylindrical);
    everything should be acheivable with the first two arguments.
    
    Returns
    -------
    orOutList: list of rectangular orientations defined by [ax, ay, az, bx, by, bz, ox, oy, oz] 
               where (ax,ay,az) is the transformed x-axis, (bx,by,bz) is the transformed y-axis
               and (ox,oy,oz) is the transformed origin.
    
    '''
    orOutList = []
    for ii in range(len(phis)):
        phi = phis[ii]
        xAx = np.array([1,0,0])
        yAx = np.array([0,1,0])
        
        R = np.array([[np.cos(phi), -np.sin(phi), 0],
                    [np.sin(phi),  np.cos(phi), 0],
                    [0,                  0, 1]])
        
        xPrime = np.matmul(R,xAx)
        yPrime = np.matmul(R,yAx)
        orOut = np.hstack((xPrime, yPrime, origin[ii], addRotDim, addRotAngle))
        orOutList.append(orOut)
    return orOutList
    
def plot2Dmesh(nodes, elements, returnFig=False):
    '''
    Function to plot a graph of a 2D mesh. WARNING this will break for large
    number of elements.
    
    Parameters
    ----------
    nodes : array, float
        The coordinates of the nodes passed in an array of shape
        (nDims, nNodes). It must have 2 coordinates for each node.
        
    elements : array, int
        The defintion of the elements in the shape (nNodes per element,
        nElements). Each column should have the nodes that make up that
        element and should be 1 indexed.
        
    returnFig : boolean, optional
        Whether or not to return the generated figure. Default is False in
        which case it is not returned.
    '''
    nDims, nNodes = np.shape(nodes)
    if nDims != 2:
        raise ValueError('This will only plot for 2D.')
        
    nNodesPerElement, nElements = np.shape(elements)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.scatter(nodes[0], nodes[1], c='b', s=10.0)
    ax.set_aspect('equal')
    dx = np.max(nodes[0]) - np.min(nodes[0])
    dy = np.max(nodes[1]) - np.min(nodes[1])
    ax.set_xlim([np.min(nodes[0])-0.05*dx, np.max(nodes[0])+0.05*dx])
    ax.set_ylim([np.min(nodes[1])-0.05*dy, np.max(nodes[1])+0.05*dy])
    
    for c1 in range(nElements):
        for c2 in range(nNodesPerElement):
            ax.plot([nodes[0, elements[c2, c1]-1], nodes[0, elements[(c2+1)%nNodesPerElement, c1]-1]],
                    [nodes[1, elements[c2, c1]-1], nodes[1, elements[(c2+1)%nNodesPerElement, c1]-1]],
                     'k')
    
    if returnFig == True:
        return fig
    else:
        return