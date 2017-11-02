'''
Author : Alexander Ballisat

Script of various functions that are used in Pogo geometry creation and other
functions functions needed for these. These should provide a good basis for
creating most models.
'''

import numpy as np
import pogoFunctions as pf
import matplotlib.pyplot as plt

def _facetLine_(string):
    '''
    NEEDS COMMENTING
    '''
    string = string.split(' ')
    n = int(string[0])
    return [int(a) for a in string[1:1+n]]

def _genString_(nodeNumbers, nMax=1024):
    '''
    NEEDS COMMENTING
    '''
    n = len(nodeNumbers)
    tmp = '{}'.format(n)
    tmp2 = ' {}'*n
    tmp2 = tmp2.format(*nodeNumbers)
    tmp += tmp2 + '\n'
    tmp = splitString(tmp, nMax, '    ')
    return tmp

def _normaliseVector_(vector):
    '''
    Function to normalise a vector
    
    Parameters
    ----------
    vector : array, float
        The vector to be normalised. May have any length.
        
    Returns
    -------
    vector : array, float
        The normalised vector.
    '''
    return vector*1./np.sqrt(np.sum(vector*vector))
    
def _plotLine_(nodes, a, b):
    '''
    NEEDS COMMENTING
    '''
    a -=1
    b -=1
    x = [nodes[a, 1], nodes[b, 1]]
    y = [nodes[a, 2], nodes[b, 2]]
    plt.plot(x,y, 'o-')
    
def _plotPath_(nodes, path, inds=[0,1]):
    '''
    NEEDS COMMENTING
    '''
    path = [a-1 for a in path]
    print path
    x = nodes[path, inds[0]+1]
    y = nodes[path, inds[1]+1]
    plt.plot(x,y,'o-')
    return
    
def _solveQuadratic_(a,b,c):
    '''
    Function to solve a quadratic equation of form ax^2 + bx + c = 0. Returns
    the minimum value or None is there is no real solution. Designed to be
    used in ray tracing calculations.
    '''
    det = b*b - 4*a*c
    
    if det < 0:
        return None
        
    else:
        x1 = (-1.*b + np.sqrt(det))/(2.*a)
        x2 = (-1.*b - np.sqrt(det))/(2.*a)
        return min([x1,x2])
     
        
def block3DPolyFormat(x,y,z,origin=[0.,0.,0.]):
    '''
    Function to generate the nodes and faces of a 3D block in the format used
    in a .poly file.
    
    Parameters
    ----------
    x : float
        The x extent of the block.
        
    y : float
        The y extent of the block.
        
    z : float
        The z extent of the block.
        
    origin : iterable, float, optional
        The coordinate of the node nearest the origin of the global
        coordinate system. The default is the origin of the global coordinate
        system, [0., 0., 0.].
    
    Returns
    -------
    nodes : array, float
        The coordinates of the nodes of the block.
        
    faces : array, int
        The indices of the nodes that make the faces of the block.
    '''
    nodes = np.array([origin,
               [x, origin[1], origin[2]],
               [origin[0], y, origin[2]],
               [x, y, origin[2]],
               [origin[0], origin[1], z],
               [x, origin[1], z],
               [origin[0], y, z],
               [x, y, z]])
             
    faces = np.array([[0,1,3,2], 
              [4,5,7,6],
              [0,4,6,2], 
              [1,5,7,3], 
              [0,1,5,4],
              [2,3,7,6]]) + 1
    
    return nodes, faces

def circleGen(centre,r,dTheta=None,nPoints=None):
    '''
    Function to generate the points on a circle its centre and radius along
    with one other parameter, either dTheta or nPoints. One of the latter two
    must be set.
    
    Use for circles, for cylinders use cylinderGeneration as it is more
    general.
    
    Input
    -----
    centre : iterable
        An iterable containing atleast the x and y coordinates of the centre of
        the circle. It can contain more, they will be unused.
        
    r : float
        The radius of the circle.
        
    dTheta : float, optional
        The change in angle between two points on the boundary of the circle.
        The default is None in which case nPoints is used to determine the
        change in angle.
        
    nPoints : int, optional
        The number of points on the boundary of the circle. The dedault is
        None in which case dTheta determines the number of points.
        
    Output
    -----
    points : float array
        Array of shape (2, nPoints) of the coordinates of the circumference of
        the circle.
        
    facets : int array
        Array of shape (4, nPoints) that describes the facets of the sides that
        make a cylinder. The z coordinates of such points need to be set
        afterwards.
    '''
    
    if dTheta == None and nPoints == None:
        raise ValueError('\nOne of dTheta or nPoints must be set\n')
        
    if len(centre) < 2:
        raise ValueError('\nThe coordinates of the centre of the circle were not set')
        
    if nPoints == None:
        nPoints = int(np.pi*2/dTheta)
        
    if dTheta == None:
        dTheta= 2*np.pi/nPoints
        
    points = np.zeros((2, nPoints))
    c1 = np.linspace(0, nPoints-1, nPoints)*dTheta
    
    points[0] = centre[0] + r*np.cos(c1)
    points[1] = centre[1] + r*np.sin(c1)
    
    facets = np.zeros((4,nPoints), dtype=int)
    c1 = (np.linspace(1, nPoints, nPoints)).astype(int)
    
    facets[0] = np.copy(c1)
    facets[1] = (np.copy(c1) - 1 + 1) % nPoints + 1
    facets[2] = np.copy(facets[1]) +nPoints
    facets[3] = np.copy(c1) + nPoints
    
    return points, facets

def cylinderGeneration(r, z, theta=None, dTheta=None, nPoints=None,
                       origin=[0.,0.,0.], theta0=0.0):
    '''
    Funciton to generate the nodes and facets for a cylinder or part of a
    cylinder. The cylinder is generated such that its length is in the z
    axis. One of dTheta and nPoints must be specified.
    
    Parameters
    ----------
    r : float
        The radius of the cylinder.
        
    z : float
        The height of the cylinder.
        
    theta : float, optional
        The proportion of the cylinder to generate in degrees. Default is
        None in which case a complete cylinder, 360 degrees is generated.
        
    dTheta : float, optional
        The change in angle between two points on the boundary of the circle.
        The default is None in which case nPoints is used to determine the
        change in angle.
        
    nPoints : int, optional
        The number of points on the boundary of the circle. The dedault is
        None in which case dTheta determines the number of points.
        
    origin : iterable, float
        The global origin of the centre of the base of the cylinder.
        
    theta0 : float, optional
        The start of the circle around the z axis anticlockwise from [r,0].
        Useful for generating arcs at arbitrary points around the z axis.
        Angle is in degrees.
        
    Output
    -----
    points : float array
        Array of shape (2, nPoints) of the coordinates of the circumference of
        the circle.
        
    facets : int array
        Array of shape (4, nPoints) that describes the facets of the sides that
        make a cylinder. The z coordinates of such points need to be set
        afterwards.
    '''
    if dTheta == None and nPoints == None:
        raise ValueError('\nOne of dTheta or nPoints must be set\n')
        
    if len(origin) != 3:
        raise ValueError('\norigin must have 3 coordinates')
    
    if theta == None:
        theta = 360.
        
    theta = np.deg2rad(theta)
        
    if theta > 2.*np.pi:
        raise ValueError('theta must be in (0, 2.*pi]')
        
    if nPoints == None:
        nPoints = int(theta/dTheta)
    
    if dTheta == None:
        dTheta= theta/nPoints
    
    theta0_ = np.deg2rad(theta0)
    
    if theta != 2.*np.pi:
        nPoints += 1
    
    points = np.zeros((3, nPoints*2))    
    c1 = np.linspace(0, nPoints-1, nPoints)*dTheta + theta0_
    
    xs = origin[0] + r*np.cos(c1)
    ys = origin[1] + r*np.sin(c1)
    points[0] = np.hstack((xs, xs))
    points[1] = np.hstack((ys, ys))
    points[2] = np.hstack((np.ones(nPoints)*(origin[2]+z), np.ones(nPoints)*origin[2]))
    
    facets = np.zeros((4,nPoints), dtype=int)
    c1 = (np.linspace(1, nPoints, nPoints)).astype(int)
    
    facets[0] = np.copy(c1)
    facets[1] = (np.copy(c1)) % nPoints + 1
    facets[2] = np.copy(facets[1]) + nPoints
    facets[3] = np.copy(c1) + nPoints
    
    if theta != 2.*np.pi:
        facets = facets[:,:-1]
    
    return points, facets
    
def delayLineOnBlock2DPolyFile(fileName, x1, x2, y1, y2, angle,
                               probeCentre=None, origin=[0., 0.], size1=None,
                               size2=None):
    '''
    Function to generate a poly file for a delay line (similar to a wedge) on
    top of a block. NEEDS COMMENTING
    
    Parameters
    ----------
    x1 : float
        The width of the lower block.
        
    x2 : float
        The width of the delay line
        
    y1 : float
        The height of the lower block
        
    y2 : float
        The length of the delay line as if it were a rectangle (maximum length)
        
    angle : float
        The angle of the delay line relative to the negative y-axis, in
        degrees, in the anti-clockwise direction
    
    probeCentre : iterable, float
        The centre of the delay line on the surface relative to the origin,
        easiest to think of as if it were vertical
    '''
    if probeCentre == None:
        probeCentre = [x1/2., y1]
    
    theta_ = np.deg2rad(angle)
    widthOnFace = x2/np.cos(theta_)
    ## Define the lower block nodes - cyclic anticlockwise from origin
    lowerNodes = [origin,
                  [origin[0] + x1, origin[1]],
                  [origin[0] + x1, origin[1] + y1],
                  [origin[0] + probeCentre[0] + x2/2., origin[1] + y1],
                  [origin[0] + probeCentre[0] + x2/2. - widthOnFace, origin[1] + y1],
                  [origin[0], origin[1] + y1]]
    
    ## Define the lower block segments
    lowerSegments = [[1,2],
                     [2,3],
                     [3,4],
                     [4,5],
                     [5,6],
                     [6,1]]
    
    ## Define the delay line nodes - only the top two needed
    centreOfRot = np.array([origin[0] + probeCentre[0] + x2/2., origin[1]+y1])
    rightNode = np.array([origin[0] + probeCentre[0] + x2/2., origin[1]+y1+y2])
    leftNode = np.array([origin[0] + probeCentre[0] - x2/2., origin[1]+y1+y2])
    
    upperNodes = [pf.rotate2D(rightNode, angle, centreOfRot),
                  pf.rotate2D(leftNode, angle, centreOfRot)]
           
    ## Define the delay line segments
    upperSegments = [[4,7],
                     [7,8],
                     [8,5]]
                 
    corners = lowerNodes + upperNodes
    segments = lowerSegments + upperSegments
    
    if fileName[-5:] != '.poly':
        fileName += '.poly'
    
    with open(fileName, 'w') as f:
        # Write out the corners
        f.write('# <Number of vertices> <Number of dimensions> <Number of attributes> <Boundary markers 0 or 1>\n')
        f.write('{} 2 2 0\n'.format(len(corners)))
        
        for c1 in range(len(corners)):
            string = '{}'.format(c1+1)
            string += ' {} {}\n'.format(*corners[c1])
            f.write(string)
        
        # Write out the segments
        f.write('# <Number of segments> <Boundary markers 0 or 1>\n')
        f.write('{} 0\n'.format(len(segments)))
        
        for c1 in range(len(segments)):
            string = '{}'.format(c1+1)
            string += ' {} {}\n'.format(*segments[c1])
            f.write(string)
        
        # Write out the holes - none in this case
        f.write('0\n')
        
        #Write out the number of regions and the region attributes
        if size1 == None:
            size1 = -1
        if size2 == None:
            size2 = -1
            
        f.write('2\n')
        f.write('1 {} {} 0 {:.12f}\n'.format(x1/2, y1/2, size1))
        point2_ = [origin[0] + probeCentre[0], origin[1]+probeCentre[1]+y2/2.]
        point2 = pf.rotate2D(point2_, angle, centreOfRot)
        f.write('2 {} {} 1 {:.12f}\n'.format(point2[0], point2[1], size2))
        
        
    return

def rayCylinderIntersectionPoint(rayOrigin, rayVector, cylinderCentre,
                                 cylinderRadius, cylinderVector=None,
                                 cylinderHeight=None):
    '''
    Function to find the intersection of a ray with an arbitrarily angled
    cylinder.
    
    NEEDS COMMENTING
    '''
    if cylinderVector == None:
        cylinderVector = np.array([0., 0., 1.])
        
    else:
        cylinderVector = _normaliseVector_(cylinderVector)
    
    tmp1 = rayVector - np.dot(rayVector, cylinderVector)*cylinderVector
    dp = rayOrigin - cylinderCentre
    tmp2 = dp - np.dot(dp, cylinderVector)*cylinderVector
    
    A = np.dot(tmp1, tmp1)
    B = 2.*np.dot(tmp1, tmp2)
    C = np.dot(tmp2, tmp2) - cylinderRadius**2
    
    infCylinderIntersectionT = _solveQuadratic_(A, B, C)
    
    if infCylinderIntersectionT == None or infCylinderIntersectionT < 0:
        ## No intersection
        print 'No intersection found.'
        return None
    
    elif cylinderHeight == None:
        ## Infinite cylinder
        return rayOrigin + rayVector*infCylinderIntersectionT
        
    else:
        ## Finite cylinder
    
        # Define the caps
        pTop = cylinderCentre + cylinderVector*cylinderHeight/2.
        pBottom = cylinderCentre - cylinderVector*cylinderHeight/2.
        vTop = np.copy(cylinderVector)
        vBottom = np.copy(cylinderVector)*-1.
        
        # Test if intersection lies within the caps
        if np.dot(cylinderVector, rayOrigin+infCylinderIntersectionT*rayVector-pBottom) > 0:
            if np.dot(cylinderVector, rayOrigin+infCylinderIntersectionT*rayVector-pTop) < 0:
                infCylinderIntersectionT *= -1
        
        # Test if intersection is on caps
        tTop = rayPlaneIntersection(rayOrigin, rayVector, pTop, vTop)
        tBottom = rayPlaneIntersection(rayOrigin, rayVector, pBottom, vBottom)

        ts = np.array([tTop, tBottom, infCylinderIntersectionT])
        try:
            tMin = np.min(ts[ts>0])
        except:
            tMin = None
        
        return rayOrigin + tMin*rayVector

def rayPlaneIntersection(rayOrigin, rayVector, planePoint, planeVector,
                         returnTime=True):
    '''
    Function to calculate the ray plane intersection using the point normal
    method.
    
    NEEDS COMMENTING
    '''
    rayVector = _normaliseVector_(rayVector)
    planeVector = _normaliseVector_(planeVector)
    
    denom = np.dot(rayVector, planeVector)
    if denom == 0.:
        print 'ray is parallel to plane'
        return None
    else:
        t = np.dot(planePoint-rayOrigin, planeVector)/denom
        
    if returnTime == True:
        return t
    
    else:
        return rayOrigin + rayVector*t
        
    
def raySphereInstersetion(circleOrigin, circleRadius, rayOrigin, rayPath):
    '''
    Function to find the intersection points of a ray with a sphere or circle
    
    NEEDS COMMENTING
    '''
    if (len(circleOrigin) != len(rayOrigin) and len(rayOrigin) != len(rayPath)):
        raise ValueError('Circle origin, ray origin and ray path must all have the same number of dimensions')
    
    nDims = len(circleOrigin)
    B = 0
    C = 0
    for c1 in range(nDims):
        B += 2. * rayPath[c1] * (rayOrigin[c1] - circleOrigin[c1])
        C += np.power(rayOrigin[c1] - circleOrigin[c1], 2)
        
    C -= circleRadius*circleRadius
    
    t1 = (-1.*B - np.sqrt(B*B-4.*C))/2.
    t2 = (-1.*B + np.sqrt(B*B-4.*C))/2.
    
    tmin = min((t1, t2))
    
    intersectionPoint = rayOrigin + rayPath * tmin
    
    return intersectionPoint
    
def pointsNotInCrack2D(leftIntersection, rightIntersection, circlePoints):
    '''
    Function to find the points which are not region where the defect joins 
    the circular hole.
    
    IMPORTANT Assumes that the crack propagates in the +y direction and that
    the leftIntersection has a larger x coordinate than rightIntersection. 
    Also assumes that the width of the crack is less than the radius of the 
    circle
    '''
    
    nCirclePoints = len(circlePoints[0])
    print 'Number of points on circle = {}\n'.format(nCirclePoints)
    inCircle = np.ones(nCirclePoints) # will get set to 0 if it is
    inCircle = inCircle.astype(int)
    pointInds = np.linspace(0, nCirclePoints-1, nCirclePoints)
    
    for c1 in range(0, nCirclePoints):        
        if (leftIntersection[0] > circlePoints[0,c1] > rightIntersection[0]
        and leftIntersection[1] < circlePoints[1,c1]
        and rightIntersection[1] < circlePoints[1,c1]):
            inCircle[c1] = 0
    
    outCircle = np.where(inCircle != 0)[0]
    outCircle = outCircle.astype(int)
    nRemoved = nCirclePoints - len(outCircle)
    print 'Number of points removed = {}\n'.format(nRemoved)
    
    if nRemoved != 0:
        firstRemoved = np.where(inCircle==0)[0][0]
        lastRemoved = np.where(inCircle==0)[0][0]
        rhs = pointInds[lastRemoved+1:]
        lhs = pointInds[:firstRemoved]
        orderedBoundary = np.hstack((rhs, lhs))
    
    else:
        firstRemoved = None
        lastRemoved = None
        orderedBoundary = pointInds
        
    return orderedBoundary, nRemoved
    
def splitString(string, n, insertString=None):
    '''
    Function to split a string into chunks and insert something at the start
    the next line
    
    Input
    -----
    string : string
        The string to be split up
        
    n : int
        The number of characters in each line
        
    insertString : string
        The text to be inserted into each new line created
        
    Output
    -----
    newString : string
        The final formatted string
    '''
    
    newString = ''
    nChar = len(string)
    nInsert = len(insertString)
    insertString = '\n' + insertString
    chunkSize = n - nInsert
    count = 0
    while count < nChar:
        if count > nChar-chunkSize:
            toAdd = string[count:]
        else:
            toAdd = string[count:count+chunkSize]
            while toAdd[-1].isalnum() == True:
                toAdd = toAdd[:-1]
        if count == 0:
            newString = toAdd
        else:
            newString = insertString.join((newString, toAdd))
        count += len(toAdd)
    
    return newString

def rect2DPolyFormat(x,y,origin=[0.,0.,0.]):
    '''
    Function to generate the information of a rectangle. Needs documenting.
    '''
    corners = [origin,
               [origin[0] + x, origin[1]],
               [origin[0] + x, origin[1] + y],
               [origin[0], origin[1] + y]]
               
    segments = [[1,2],
                [2,3],
                [3,4],
                [4,1]]
                
    return corners, segments
    
def write2DRectanglePolyFile(x,y,fileName,origin=[0.,0.]):
    '''
    Function to generate a .poly file of a 2D rectangle.
    
    Parameters
    ----------
    x : float
        The x extent of the block.
        
    y : float
        The y extent of the block.
        
    fileName : string
        The name of the .poly file to be saved.
        
    origin : iterable, float, optional
        The coordinate of the node nearest the origin of the global
        coordinate system. The default is the origin of the global coordinate
        system, [0., 0., 0.]
    
    Returns
    -------
    None
    
    '''
    if fileName[-5:] != '.poly':
        fileName += '.poly'
    
    corners, segments = rect2DPolyFormat(x,y,origin)
    
    with open(fileName, 'w') as f:
        # Write out the corners
        f.write('# <Number of vertices> <Number of dimensions> <Number of attributes> <Boundary markers 0 or 1>\n')
        f.write('{} 2 0 0\n'.format(len(corners)))
        
        for c1 in range(len(corners)):
            string = '{}'.format(c1+1)
            string += ' {} {}\n'.format(*corners[c1])
            f.write(string)
        
        # Write out the segments
        f.write('# <Number of segments> <Boundary markers 0 or 1>\n')
        f.write('{} 0\n'.format(len(segments)))
        
        for c1 in range(len(segments)):
            string = '{}'.format(c1+1)
            string += ' {} {}\n'.format(*segments[c1])
            f.write(string)
        
        # Write out the holes - none in this case
        f.write('0\n')
        
        #Write out the number of regions
        f.write('0\n')
        
    return

def write3DBlockPolyFile(x,y,z,fileName,origin=[0.,0.,0.]):
    '''
    Function to generate a .poly file of a 3D block.
    
    Parameters
    ----------
    x : float
        The x extent of the block.
        
    y : float
        The y extent of the block.
        
    z : float
        The z extent of the block.
        
    fileName : string
        The name of the .poly file to be saved.
        
    origin : iterable, float, optional
        The coordinate of the node nearest the origin of the global
        coordinate system. The default is the origin of the global coordinate
        system, [0., 0., 0.]
    
    Returns
    -------
    None
    '''
    
    if fileName[-5:] != '.poly':
        fileName += '.poly'
    
    nodes, faces = block3DPolyFormat(x,y,z,origin)
    
    with open(fileName, 'w') as f:
        # Nodes
        f.write('# <Number of nodes> <Number of dimensions>\n')
        f.write('8 3\n')
        for c1 in range(1, 9):
            s = '{} '.format(c1) + '{} {} {}\n'.format(*nodes[c1])
            f.write(s)
        
        # Faces
        f.write('\n# <Number of facets> <Boundary markers 0 or 1>\n' )
        f.write('6 0\n')
        for c1 in range(0, len(faces)):
            f.write('1 0\n')
            f.write('4 {} {} {} {}\n'.format(*faces[c1]))
            
        # Holes
        f.write('# <Number of holes>\n')
        f.write('0\n')
        
        # Regions
        f.write('# <Number of regions>\n')
        f.write('0')
    return
    
def writeDelayLineOnPartCylinderWithCrackPolyFile(fileName, 
                                                  innerRadius,
                                                  outerRadius,
                                                  arcAngle,
                                                  blockHeight,
                                                  probeWidth,
                                                  probeHeight,
                                                  probeVAngle,
                                                  probeRotation,
                                                  probeCentre,
                                                  size1,
                                                  size2,
                                                  crack=False,
                                                  crackLength=None,
                                                  crackHeight=None,
                                                  crackRotation=None,
                                                  returnTransProperties=False,
                                                  blockRotation=None):
    '''
    Function to write the poly file for a  delay line on part of a chunk of
    a cylinder. The idea is to not model the whole cylinder but only a chunk
    of it with the inner radius and a crack coming off the inner radius.
    
    Parameters
    ----------
    
    Returns
    -------
    None
    '''
    if fileName[-5:] != '.poly':
        fileName += '.poly'
    
    nodeCount = 0 # keep track of the number of nodes so far
    tempStore = np.copy(probeCentre)
    
    if blockRotation == None:
        blockRotation = arcAngle/2.
    ##### Inner radius - the upper nodes are always first
    nodes1, faces1 = cylinderGeneration(innerRadius,
                                        blockHeight,
                                        dTheta=np.pi/72.,
                                        theta=arcAngle,
                                        theta0=-1.*blockRotation)
    nTop = len(nodes1[0])/2
    innerCorners = [1, nTop, nTop+1, 2*nTop]
    nodeCount += len(nodes1[0])
    innerBaseNodes = np.linspace(1, nTop, nTop)+nTop
    innerTopNodes = np.linspace(1, nTop, nTop)
    
    ##### Outer radius
    nodes2, faces2 = cylinderGeneration(outerRadius,
                                        blockHeight,
                                        dTheta=np.pi/72.,
                                        theta=arcAngle,
                                        theta0=-1.*blockRotation)
    nTop = len(nodes2[0])/2
    outerCorners = [1, nTop, nTop+1, 2*nTop]
    outerCorners = [a+nodeCount for a in outerCorners]
    faces2 += nodeCount
    
    outerBaseNodes = np.linspace(1, nTop, nTop)+nTop+nodeCount
    outerBaseNodes = outerBaseNodes[::-1]
    outerTopNodes = np.linspace(1, nTop, nTop)+nodeCount
    outerTopNodes = outerTopNodes[::-1]
    
    nodeCount += len(nodes2[0])
    
    ##### Bottom face
    bottom = np.hstack((innerBaseNodes, outerBaseNodes))
    
    ##### Side faces
    sides = np.array([[innerCorners[0], innerCorners[2],
                           outerCorners[2], outerCorners[0]],
                           [innerCorners[1], innerCorners[3],
                           outerCorners[3], outerCorners[1]]]).T
    
    ##### Top cylinder
    nodes3, faces3 = cylinderGeneration(probeWidth/2., 
                                         probeHeight,
                                         nPoints=360, 
                                         origin=probeCentre)
    
    nTop = len(nodes3[0])/2
    
    ## The bottom face of the cylinder rotation - elongate
    xMin = np.min(nodes3[0])
    nodes3[0,nTop:] = xMin + (nodes3[0,nTop:] - xMin)/np.cos(np.deg2rad(probeVAngle))
    
    ## Rotate the top nodes
    rotNodes = np.vstack((np.copy(nodes3[0,:nTop]), np.copy(nodes3[2, :nTop])))
    rotNodes = pf.rotate2D(rotNodes, -1.*probeVAngle, np.array([xMin, blockHeight]))
    nodes3[0, :nTop] = np.copy(rotNodes[0])
    nodes3[2, :nTop] = np.copy(rotNodes[1])
    
    ## Rotate the whole face of the transducer - Make sure there is enough material for the rotation
    nodes3[:2, :] = pf.rotate2D(nodes3[:2, :], probeRotation, np.array([xMin, tempStore[1]]))
    
    probeCentre[:2] = pf.rotate2D(probeCentre[:2], probeRotation, np.array([xMin, tempStore[1]]))
    
    faces3 += nodeCount
    
    probeBottom = np.linspace(1, nTop, nTop) + nTop + nodeCount
    probeTop = np.linspace(1, nTop, nTop) + nodeCount
    
    nodeCount += len(nodes3[0])
    
    ##### Top face
    topOuterPoly = np.hstack((innerTopNodes, outerTopNodes))
    topInnerPoly = np.copy(probeBottom)
    
    ##### Crack
    if crack == True:
        ## Check information supplied
        if crackLength == None and crackHeight == None:
            raise ValueError('One of crackLength and crackHeight must be specified.')
            
        if crackLength != None and crackHeight == None:
            crackHeight = crackLength
            
        elif crackHeight != None and crackLength == None:
            crackLength = crackHeight
            
        if crackLength > outerRadius - innerRadius:
            raise ValueError('Crack is wider than the part. Your part is really broken.')
        
        if crackHeight > blockHeight:
            raise ValueError('Crack is taller than the part. Your part is really broken.')
            
        ## Set crack origin - the position of the tip
        crackOrigin = np.array([innerRadius + crackLength, 0., 0.])
        if crackRotation != None:
            crackOrigin[:2] = pf.rotate2D(crackOrigin[:2], crackRotation,
                                          np.array([innerRadius, 0.]))
        
        ## Set the crack top join to the part
        crackTopCentre = np.array([innerRadius, 0., crackHeight])
        
        ## Work out what faces to delete from the inner radius and delete them
        # Find where the centre nodes are
        pT = np.array([innerRadius, 0., blockHeight])
        pB = np.array([innerRadius, 0., 0.])
        indTop = np.where((nodes1[0]==pT[0]) & (nodes1[1]==pT[1]) &(nodes1[2]==pT[2]))[0][0]+1
        indBottom = np.where((nodes1[0]==pB[0]) & (nodes1[1]==pB[1]) &(nodes1[2]==pB[2]))[0][0]+1
        
        facetsToDelete = np.where(faces1 == indTop)[1]
        faces1 = np.delete(faces1, facetsToDelete, 1)
        
        ## Fix the crack width for now
        crackWidth = 0.1E-3
        leftVector = np.array([innerRadius-crackOrigin[0], -1.*crackWidth/2., 0.])
        rightVector = np.array([innerRadius-crackOrigin[0], crackWidth/2., 0.])
        leftSurfNorm = np.cross(pB-pT, (nodes1[:,indTop-2]-nodes1[:, indTop-1]).T)
        rightSurfNorm = np.cross(pB-pT, (nodes1[:,indTop]-nodes1[:, indTop-1]).T)

        rightIntersect = rayPlaneIntersection(crackOrigin, rightVector,
                                              nodes1[:,indTop-1], rightSurfNorm,
                                              returnTime=False)
                                              
        leftIntersect = rayPlaneIntersection(crackOrigin, leftVector,
                                             nodes1[:,indTop-1], leftSurfNorm,
                                             returnTime=False)

        crackNodes = np.vstack((crackOrigin, leftIntersect, rightIntersect,
                                crackTopCentre)).T
        crackFaces = [np.array([1,2,4])+nodeCount,
                       np.array([1,3,4])+nodeCount,
                       np.array([indTop-1, indTop, nodeCount+4, nodeCount+2, indBottom-1]),
                       np.array([indTop, indTop+1, indBottom+1, nodeCount+3, nodeCount+4])]
                                     
        # Udpdate the bottom face - this is a much neater way of doing it
        insertInd = np.where(bottom == indBottom)[0]
        bottom[insertInd] = 2+nodeCount
        bottom = np.insert(bottom, insertInd+1, np.array([1+nodeCount, 3+nodeCount]))
        
    if crack == True:                           
        allNodes = (nodes1, nodes2, nodes3, crackNodes)
    else:
        allNodes = (nodes1, nodes2, nodes3)
    allFaces = (faces1, faces2, sides, faces3) #but not the bottom or top
    nodes = np.hstack(allNodes)
    faces = np.hstack(allFaces)
    
    nodes[np.abs(nodes) < 1E-15] = 0.0
    
    nNodes = len(nodes[0])
    #special faces are top and bottom of block and probe and 4 for the crack
    nFaces = len(faces[0]) + 1 + 1 + 1 + 1 #should be +8 total
    if crack == True:
        nFaces += 4
    nBlitzWrite = len(faces[0])
    

    print 'nNodes = {}'.format(nNodes)
    print 'nFaces = {}'.format(nFaces)
    
    ##### Write it all out
    with open(fileName, 'w') as f:
        # Nodes
        f.write('# <Number of nodes> <Number of dimensions>\n')
        f.write('{} 3 0 0\n'.format(nNodes))
        for c1 in range(nNodes):
            s = '{} '.format(c1+1) + '{} {} {}\n'.format(*nodes[:, c1])
            f.write(s)
        
        # Faces
        f.write('\n# <Number of facets> <Boundary markers 0 or 1>\n' )
        f.write('{} 0\n'.format(nFaces))
        fCount = 0
        for c1 in range(0, nBlitzWrite):
            fCount += 1
            f.write('1 0\n')
            f.write('4 {} {} {} {}\n'.format(*faces[:, c1]))
            
        # Bottom face
        f.write('1 0\n')
        f.write(_genString_(bottom))
        
        # Top face of block
        f.write('2 1\n')
        f.write(_genString_(topOuterPoly))
        f.write(_genString_(topInnerPoly))        
        f.write('1 {} {} {}\n'.format(*probeCentre))
        
        # Bottom of cylinder
        f.write('1 0\n')
        f.write(_genString_(probeBottom))
        
        # Top of cylinder
        f.write('1 0\n')
        f.write(_genString_(probeTop))
        
        if crack == True:
            # Crack faces
            for c1 in range(len(crackFaces)):
                f.write('1 0\n')
                f.write(_genString_(crackFaces[c1]))
        
        # Holes
        f.write('\n# <Number of holes>\n')
        f.write('0\n')
        
        # Regions
        f.write('\n# <Number of regions>\n')
        f.write('2\n')
        if crack == True:
            f.write('1 {} {} {} 1 {:.12f}\n'.format(crackOrigin[0], crackOrigin[1],
            blockHeight/2., size1))
        else:
            midBlock = 0.5*(np.average(nodes1, axis=1) + np.average(nodes2, axis=1))
            f.write('1 {} {} {} 1 {:.15f}\n'.format(midBlock[0], midBlock[1],
            blockHeight/2., size1))
        midProbe = np.average(nodes3, axis=1)
        f.write('2 {} {} {} 2 {:.15f}\n'.format(midProbe[0], midProbe[1],
                midProbe[2], size2))
    
    if returnTransProperties == True:
        # Find the centre of the top of the probe and swing it round as before
        point = np.copy(tempStore)
        point[2] += probeHeight
        z = np.array([point[0], point[2]])
        z = pf.rotate2D(z, -1.*probeVAngle, np.array([xMin, blockHeight]))
        point[0] = z[0]
        point[2] = z[1]
        point[:2] = pf.rotate2D(point[:2], probeRotation, np.array([xMin, tempStore[1]]))
        
        # Create a normal vector for the surface
        v1 = nodes3[:,0] - nodes3[:,1]
        v2 = nodes3[:,1] - nodes3[:,2]
        vector = np.cross(v1, v2)
        
        # Return it all
        return point, _normaliseVector_(vector)
        
    else:
        return
    
def writeTwoLayerRectanglePolyFile(x, y1, y2, fileName, origin=[0.0, 0.0],
                                   size1=None, size2=None):
    '''
    Function to generate a .poly file of a 2D rectangle with 2 layers in the
    y axis of equal width.
    
    Parameters
    ----------
    x : float
        The x extent of the block.
        
    y1 : float
        The y extent of the first layer of the block.
        
    y2 : float
        The y extent of the second layer of the block.
        
    fileName : string
        The name of the .poly file to be saved.
        
    origin : iterable, float, optional
        The coordinate of the node nearest the origin of the global
        coordinate system. The default is the origin of the global coordinate
        system, [0., 0., 0.]
        
    size1 : float, optional
        The target mesh size for first layer. Default is None in which case
        this is set to -1. When the area is negative it is ignored by the
        Triangle meshing program.
    
    size12: float, optional
        The target mesh size for second layer. Default is None in which case
        this is set to -1. When the area is negative it is ignored by the
        Triangle meshing program.
    
    Returns
    -------
    None
    
    '''
    if fileName[-5:] != '.poly':
        fileName += '.poly'
    
    origin1 = [origin[0], origin[1]]
    origin2 = [origin[0], origin[1]+y1]
    corners1, segments1 = rect2DPolyFormat(x,y1,origin1)
    corners2, segments2 = rect2DPolyFormat(x,y2,origin2)
    
    corners2 = corners2[2:]
    segments2 = [[3,5],
                 [5,6],
                 [6,4]]
    
    corners = corners1 + corners2
    segments = segments1 + segments2
    
    with open(fileName, 'w') as f:
        # Write out the corners
        f.write('# <Number of vertices> <Number of dimensions> <Number of attributes> <Boundary markers 0 or 1>\n')
        f.write('{} 2 1 0\n'.format(len(corners)))
        
        for c1 in range(len(corners)):
            string = '{}'.format(c1+1)
            string += ' {} {}\n'.format(*corners[c1])
            f.write(string)
        
        # Write out the segments
        f.write('# <Number of segments> <Boundary markers 0 or 1>\n')
        f.write('{} 0\n'.format(len(segments)))
        
        for c1 in range(len(segments)):
            string = '{}'.format(c1+1)
            string += ' {} {}\n'.format(*segments[c1])
            f.write(string)
        
        # Write out the holes - none in this case
        f.write('0\n')
        
        #Write out the number of regions and the region attributes
        if size1 == None:
            size1 = -1
        if size2 == None:
            size2 = -1
            
        f.write('2\n')
        f.write('1 {} {} 1 {:.12f}\n'.format(x/2, y1/2, size1))
        f.write('2 {} {} 2 {:.12f}\n'.format(x/2, y1+y2/2, size2))
        
        
    return
            

    


    
    