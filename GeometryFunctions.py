'''
Script of various functions that are used in Pogo geometry creation and other
functions
'''

import numpy as np
import pogoFunctions as pf

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
    
    
def circleGen(centre,r,dTheta=None,nPoints=None):
    '''
    Function to generate the points on a circle its centre and radius along
    with one other parameter, either dTheta or nPoints. One of the latter two
    must be set.
    
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
    
def sphereRayInstersetion(circleOrigin, circleRadius, rayOrigin, rayPath):
    '''
    Function to find the intersection points of a ray with a sphere or circle
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
        f.write('1 {} {} 0 {:.12f}\n'.format(x/2, y1/2, size1))
        f.write('2 {} {} 1 {:.12f}\n'.format(x/2, y1+y2/2, size2))
        
        
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
            s = '{} '.format(c1) + '{} {} {}\n'.format(*nodes[c1-1])
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
                  [origin[0] + probeCentre[0] + widthOnFace/2., origin[1] + y1],
                  [origin[0] + probeCentre[0] - widthOnFace/2., origin[1] + y1],
                  [origin[0], origin[1] + y1]]
    
    ## Define the lower block segments
    lowerSegments = [[1,2],
                     [2,3],
                     [3,4],
                     [4,5],
                     [5,6],
                     [6,1]]
    
    ## Define the delay line nodes - only the top two needed
    centreOfRot = np.array([origin[0] + probeCentre[0] + widthOnFace/2., origin[1]+y1]) 
    rightNode = np.array([origin[0] + probeCentre[0] + widthOnFace/2., origin[1]+y1+y2])
    leftNode = np.array([origin[0] + probeCentre[0] - widthOnFace/2., origin[1]+y1+y2])
    
    upperNodes = [pf.rotate2D(rightNode, angle, centreOfRot),
                  pf.rotate2D(leftNode, angle, centreOfRot)]
    
    lowerNodes[4] = [origin[0] + probeCentre[0] + widthOnFace/2. - widthOnFace/np.cos(theta_),
                     origin[1]+y1] #Update the intersection point
           
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
    
    