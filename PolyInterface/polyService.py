# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 15:39:56 2017
Rearranging of dxfToPoly functions into a class which contains information for writing a poly file, 
as well as plotting the PSLG for sanity checks. 

Planned as a parent class for a polyGui (not yet implemented).

Methods are:
createEmptyGraph: Default in case of incorrect file path, or in case user plans to create one from scratch
polyLinesToPSLG: Converts dxf lwpolylines to PSLG (vertices and edges)
writePoly2D: Writes the poly class to a .poly file
findOuterBoundaryIndex: From lwpolylines, finds the (probable) outer boundary. NOT FULLY FUNCTIONAL
findPlines: Finds entites and stores them as lwpolylines. In case of circles/arcs, calls method to convert to lwpolyline.
closeLine: Closes lwpolyline. Not usually needed.
convertToPolyline: converts circles and arcs to lwpolylines
findHoles: Finds holes in the PSLG based on points in the dxf file.
plotPoly: Plots the PSLG
readPoly: Reads in an existent poly file

@author: Callum
"""

#Make simple class for poly, and then make polyFile a service class.
import matplotlib.path as mplPath
import warnings
import numpy as np
import matplotlib.pyplot as plt
import pdb
from matplotlib.collections import LineCollection 
from itertools import compress
   
def polylinesToPSLG(pLines,isClosed,indexOfBoundary):
    jj=1
    numberOfVertices = sum(len(vertexList) for vertexList in pLines)
    vertices = np.zeros([numberOfVertices,2])
    boundaryFlags = np.zeros([numberOfVertices,2])
    edges = []
    for ii in range(len(pLines)):            
        for kk in range(len(pLines[ii].vertices)):
            if ii==indexOfBoundary:
                vertices[jj-1,:] = [pLines[ii].vertices[kk][0],pLines[ii].vertices[kk][1]]
                boundaryFlags[jj-1] = 1
            else:
                vertices[jj-1,:] = [pLines[ii].vertices[kk][0],pLines[ii].vertices[kk][1]]
                boundaryFlags[jj-1] = 0
            jj+=1
    jj = 1
    for ii in range(len(pLines)):
        jStart = jj
        for kk in range(len(pLines[ii])-1):
            edges.append([jj,jj+1])
            jj+=1
        if isClosed[ii]:
            edges.append([jj,jStart])
            jj+=1
    edges = np.array(edges)
    return vertices,boundaryFlags.astype(int),edges

def writePoly2d(polyInstance,filePath):
    '''
	 Function to write information loaded from the .dxf file to a .poly file.
	 Parameters
     ----------
     pLines : list
	    	List of matplotlib.path object containing polyline information obtained from .dxf file.
     isClosed : dict
		Dictionary of booleans corresponding to the pslg entities. 'true' forces the 
		pslg to close the entity, 'false' does not.
	boundaryKey : string
		String representing the pslg key corresponding to the external boundary
	filePath : string
		String indicating the destination to save the .poly file.
        
    Returns
    -------
	  0, regardless of success or failure. May implement error codes in the future.
    '''
#        filePath = filePath[:-4]
    if filePath[-5:]!='.poly':
        filePath += '.poly'
    filePath = filePath
    with open(filePath,'w') as f:
        f.write('<# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>\n')
        f.write("{} 2 ".format(polyInstance.numberOfVertices) + "0 {} \n".format(polyInstance.numberOfBoundaryVertices))
        for ii in range(polyInstance.numberOfVertices):
            polyInstance.vertexIDs[ii].tofile(f," ")
            f.write(" ")
            polyInstance.vertices[ii,:].tofile(f," ")
            #polyInstance.boundaryFlags[ii].tofile(f," ",format = "%d")
            f.write("\n")
        
        f.write("{} 0\n".format(polyInstance.numberOfEdges))
        for ii in range(polyInstance.numberOfEdges):
            polyInstance.edgeIDs[ii].tofile(f," ")
            f.write(" ")
            polyInstance.edges[ii,:].tofile(f," ")
            f.write("\n")
        
        f.write("{}\n".format(polyInstance.numberOfHoles))
        
        for ii in range(polyInstance.numberOfHoles):
            polyInstance.holeIDs[ii].tofile(f," ")
            f.write(" ")
            polyInstance.holes[ii,:].tofile(f," ")
            f.write("\n")
            
    return 0    

    
def findOuterBoundaryIndex(pLines):
    '''
	Function to automatically determine the external boundary in the PSLG. 
	Essentially checks for any of the entities which encloses all other entities.
	
	Parameters
    -------
	pLines : list
		List of matplotlib.path object containing polyline information obtained from .dxf file.
	
	Returns
    -------
    ii : integer
        The index corresponding to the enclosing boundary in the PSLG
    '''
    for ii in range(len(pLines)):
        jj = 0
        while np.all(pLines[ii].contains_points(pLines[jj].vertices)) or jj==ii:
            if jj == len(pLines)-1:
                return ii
            jj+=1
    errStr = "Must have one entity (LWPOLYLINE,CIRCLE) defining the outer boundary."
    raise ValueError(errStr)
    
def findPlines(entities,elementSize,precision = False):
    '''
	Function to determine the PSLG entites from the .dxf entities. 
	Reads LWPOLYLINE and LINE directly, and converts CIRCLES and 
	ARCS into LWPOLYLINES based on elementSize.
	
	Parameters
    ----------
    entities : EntitySection
        Object containing all entity information from the .dxf file.
		For detailed documentation see http://pythonhosted.org/dxfgrabber/#entity-section
	elementSize : float
		The length of the LWPOLYLINE segments generated from ARC and CIRCLE entities.
    Returns
    -------
    pLines : list
		List of matplotlib.path object containing polyline information obtained from .dxf file.
	isClosed : dict
		Dictionary of booleans corresponding to the pslg entities. 'true' forces the 
		pslg to close the entity, 'false' does not.		
    '''
    numberOfEntities = len(entities)
    pLines = []
    isClosed = []
    jj=0
    for ii in range(numberOfEntities):
        if entities[ii].dxftype == 'LWPOLYLINE':
            if entities[ii].is_closed:
                isClosed.append(True)
            else:
                isClosed.append(False)
            vertices = np.asarray(entities[ii].points)
            
            vertices = setPrecision(vertices,elementSize,precision)
            pLines.append(mplPath.Path(vertices))
            jj = jj+1
        elif entities[ii].dxftype == 'CIRCLE' or entities[ii].dxftype == 'ARC':
            pdb.set_trace()
            pLine = convertToPolyline(entities[ii],elementSize)
            pLine.vertices = setPrecision(pLine.vertices,elementSize,precision)
            pLines.append(pLine)
            if entities[ii].dxftype == 'CIRCLE':
                isClosed.append(True)
            else:
                isClosed.append(False)
            jj = jj+1
            
            
        elif entities[ii].dxftype == 'LINE':
            dy = entities[ii].end[1]-entities[ii].start[1]
            dx = entities[ii].end[0]-entities[ii].start[0]
            dr = np.sqrt(dy*dy + dx*dx)
            numberOfSubVertices = int(np.ceil(dr/elementSize))
            subVertices = np.zeros([numberOfSubVertices,2])
            subVertices[:,0] = np.linspace(entities[ii].start[0],entities[ii].end[0],numberOfSubVertices)
            subVertices[:,1] = np.linspace(entities[ii].start[1],entities[ii].end[1],numberOfSubVertices)
            pLine = mplPath.Path(subVertices)
            pLine.vertices = setPrecision(pLine.vertices,elementSize,precision)                                       
            pLines.append(pLine)
            isClosed.append(False)            
            jj = jj+1
        elif entities[ii].dxftype == 'POINT':
            continue
        else:
            warnStr="Currently only LWPOLYLINE,CIRCLE,ARC and LINE supported for boundary definition. Found {}".format(entities[ii].dxftype)
            warnings.warn(warnStr)
        
    successStr = "Found {} convertible entities.".format(jj)
    print(successStr)
    return (pLines,isClosed)

def getExponents(floatingPointNumbers):
    return np.round(np.log10(floatingPointNumbers)).astype(int)
    
def setPrecision(vertices,elementSize,precision):
    if not precision:
        return vertices
    exponent = getExponents(elementSize)
    return np.round(vertices*np.power(10,-exponent),precision)*np.power(10.,exponent)
    #return np.ceil((vertices/precision).astype(int)) * precision
        
#def joinPlines(pLines):
    
    
    
#    def splineToPline(entity):
#        degreeOfSpline = entity.degree
#        controlPoints = entity.control_points
#        knots = entity.knots
    
    
def closeLine(vertices):
    '''
	Function to close a polyline
	
	Parameters
    ----------
    vertices : float array
        Array of (x,y) pairs representing vertices on a polyline.
    Returns
    -------
    line : float array
        Array of (x,y) pairs representing vertices on a closed polyline.	
    '''
    closedVertices = np.zeros([len(vertices)+1,2])
    closedVertices[0:len(vertices)-1,:] = vertices
    closedVertices[-1,:] = vertices[0,:]
    return closedVertices
    
def convertToPolyline(entity,elementSize)    :
    '''
	Function to convert ARC and CIRCLE entities to polylines.
	
	Parameters
    ----------
    entity : Circle/Arc
        Entity object representing Circle or Arc. For detailed
		documentation see
		http://pythonhosted.org/dxfgrabber/#circle
		http://pythonhosted.org/dxfgrabber/#arc
	elementSize : float
		The length of the LWPOLYLINE segments generated from ARC and CIRCLE entities.
		
    Returns
    -------
    polyline : float array
        Array of (x,y) pairs representing vertices on a polyline representing the entity.	
    '''
    if entity.dxftype == 'CIRCLE':
        thetaBegin = 0
        thetaEnd = 2*np.pi
        isArc=False
    elif entity.dxftype == 'ARC':
        thetaBegin = entity.start_angle * np.pi / 180
        thetaEnd = entity.end_angle * np.pi / 180
        isArc=True
    else:
        errStr = 'Error: Invalid entity passed to convertToPolyline.'
        raise ValueError(errStr)
            
    angleRange = thetaEnd - thetaBegin            
    r = entity.radius
    dTheta = np.arcsin(elementSize/(2 * r)) * 2
    nTheta = int(np.ceil(angleRange/dTheta))
    
    theta = np.linspace(thetaBegin,thetaEnd,nTheta)
    if isArc:
        polyline = np.zeros([nTheta,2])
        polyline[:,0] = r * np.cos(theta) + entity.center[0]
        polyline[:,1] = r * np.sin(theta) + entity.center[1]
    else:
        polyline = np.zeros([nTheta-1,2])
        polyline[:,0] = r * np.cos(theta[:-1]) + entity.center[0]
        polyline[:,1] = r * np.sin(theta[:-1]) + entity.center[1]        
    polyline = mplPath.Path(polyline)
    return polyline

def joinPlines(pLines, isClosed):
    ignoreArray = np.zeros([len(pLines),1],dtype=bool)
    for ii,pLine in enumerate(pLines):
        if ignoreArray[ii]:
            continue
        continueLoop = True
        maxmimumPotentialConnections = 0
        while continueLoop:
            maxmimumPotentialConnections = len(pLines)-ii-1
            for jj in range(ii+1,len(pLines)):
                connectionId = _getConnectionId(pLine,pLines[jj])
                if connectionId == 0 or ignoreArray[jj]:
                    maxmimumPotentialConnections-=1
                    if maxmimumPotentialConnections == 0:
                        continueLoop=False
                elif connectionId == 1:
                    pLine.vertices = np.vstack((pLine.vertices,pLines[jj].vertices[1:-1,:]))
                    isClosed[ii] = True
                    ignoreArray[jj] = True
                    continueLoop = False
                elif connectionId == 2:
                    pLine.vertices = np.vstack((pLine.vertices,pLines[jj].vertices[-2::-1,:]))
                    isClosed[ii] = True
                    ignoreArray[jj] = True
                    continueLoop=False
                elif connectionId == 3:
                    pdb.set_trace()
                    pLine.vertices = np.vstack((pLines[jj].vertices[-1:1:-1,:],pLine.vertices))
                    ignoreArray[jj] = True
                    break
                elif connectionId == 4:
                    pLine.vertices = np.vstack((pLine.vertices,pLines[jj].vertices[-1:1:-1,:]))
                    ignoreArray[jj] = True
                    break
                elif connectionId == 5:
                    pLine.vertices = np.vstack((pLines[jj].vertices[:-1,:],pLine.vertices))
                    ignoreArray[jj] = True
                    break
                elif connectionId == 6:
                    pLine.vertices = np.vstack((pLine.vertices,pLines[jj].vertices[1:,:]))
                    ignoreArray[jj] = True
                    break
                    
    pLines = list(compress(pLines,~ignoreArray))
    isClosed = list(compress(isClosed,~ignoreArray))
    return pLines,isClosed

def _getConnectionId(line1, line2):
    '''
	Function to return the type of connection between polylines
	
	Parameters
    ----------
    line1EndPoints : 
        Matplotlib.path object containing vertices of the first line
    line2EndPoints : 
        Matplotlib.path object containing vertices of the second line
    sameLine :
        Boolean showing whether
    Returns
    -------
    connectId : int
        Integer representing the type of connection between the two lines.
                -1: Unassigned (should never return this)
                0 : No connection
                1 : Line1 and Line2 are the same polyline
                2 : Line1 and Line2 are the same polyline with vertices in reverse order
                3 : The first vertex of Line1 is coincident with the first vertex of Line2
                4 : The last vertex of Line1 is coincident with the last vertex of Line2
                5 : The first vertex of Line1 is coincident with the last vertex of Line2
                6 : The last vertex of Line1 is coincident with the first vertex of Line2
    '''
    line1EndPoints = np.array([[line1.vertices[0,:]],[line1.vertices[-1,:]]])
    line2EndPoints = np.array([[line2.vertices[0,:]],[line2.vertices[-1,:]]])
    sameLine = line1.vertices==line2.vertices or line1.vertices == line2.vertices[:,-1:None:-1]
    firstMatch = np.all(line1EndPoints[0] == line2EndPoints[0])
    lastMatch = np.all(line1EndPoints[1] == line2EndPoints[1])
    beginToEndMatch = np.all(line1EndPoints[0] == line2EndPoints[1])
    endToBeginMatch = np.all(line1EndPoints[1] == line2EndPoints[0])
    
    connectId = -1
    if not firstMatch and not lastMatch and not beginToEndMatch and not endToBeginMatch:
        connectId =  0
    elif firstMatch and lastMatch and sameLine:
        connectId =  1
    elif beginToEndMatch and endToBeginMatch and sameLine:
        connectId =  2
    elif firstMatch:
        connectId =  3
    elif lastMatch:
        connectId =  4
    elif beginToEndMatch:
        connectId =  5
    elif endToBeginMatch:
        connectId =  6
    return connectId
        
def findHoles(entities):
    numberOfEntities = len(entities)
    holeCoords = np.zeros([0,2])
    for ii in range(numberOfEntities):
        if entities[ii].dxftype == 'POINT':
            holeCoords = np.vstack((holeCoords,entities[ii].point[0:2]))
    
    holes = np.zeros([len(holeCoords),3])
    holes[:,1:] = holeCoords
    holes[:,0] = np.arange(len(holeCoords))    
    return holes


def plotPoly(polyInstance,filePath=None):    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(polyInstance.vertices[:,1],polyInstance.vertices[:,2],'r.','MarkerSize',1)
    testList = []
    for ii in range(polyInstance.numberOfEdges):
        x0=polyInstance.vertices[polyInstance.edges[ii,1]-1,1]
        y0=polyInstance.vertices[polyInstance.edges[ii,1]-1,2]
        x1=polyInstance.vertices[polyInstance.edges[ii,2]-1,1]
        y1=polyInstance.vertices[polyInstance.edges[ii,2]-1,2]
        line = ((x0,y0),(x1,y1))
        testList.append(line)
        
        
    #lineData=LineCollection((((0,1),(1,1)),((0,0),(1,1))))
    lineData = LineCollection(tuple(testList))
    pdb.set_trace()
    ax.add_collection(lineData)
    
    plt.axis('equal')
#        fig.canvas.mpl_connect('pick_event',onpick)
    plt.show()
    
    
def readPoly(filePath):
    fileContents = open(filePath,'r')    
    fileLines=fileContents.readlines()
    vertexInformation = fileLines[1].split()
    numberOfVertices = int(vertexInformation[0])
    edgeInformation = fileLines[2+numberOfVertices].split()
    numberOfEdges = int(edgeInformation[0])
    holeInformation = fileLines[3+numberOfVertices+numberOfEdges].split()
    numberOfHoles = int(holeInformation[0])
       
    vertices = _readVerticesFromFile(numberOfVertices,fileLines)    
    edges = _readEdgesFromFile(numberOfVertices,numberOfEdges,fileLines)    
    holes = _readHolesFromFile(numberOfVertices,numberOfEdges,numberOfHoles,fileLines)
    fileContents.close()
    return vertices,edges,holes

def _readVerticesFromFile(numberOfVertices,fileLines):
    vertices = np.zeros([numberOfVertices,4])
    for ii in range(numberOfVertices):
        vertString = fileLines[ii+2].split()
        vertices[ii,0] = float(vertString[0])
        vertices[ii,1] = float(vertString[1])
        vertices[ii,2] = float(vertString[2])
        if len(vertString)==3:
            vertices[ii,3] = 0
        else:
            vertices[ii,3] = float(vertString[3])
    return vertices

def _readEdgesFromFile(numberOfVertices,numberOfEdges,fileLines):
    edges = np.zeros([numberOfEdges,3])
    
    for ii in range(numberOfEdges):
        edgeString = fileLines[ii+3+numberOfVertices].split()
        edges[ii,0] = int(edgeString[0])
        edges[ii,1] = int(edgeString[1])
        edges[ii,2] = int(edgeString[2])
    edges = edges.astype(int)
    return edges

def _readHolesFromFile(numberOfVertices,numberOfEdges,numberOfHoles,fileLines):
    #pdb.set_trace()
    holes = np.zeros([numberOfHoles,3])
    for ii in range(numberOfHoles):
        holeString = fileLines[ii+4+numberOfVertices+numberOfEdges].split()
        holes[ii,0] = int(ii+1)
        holes[ii,1] = holeString[0]
        holes[ii,2] = holeString[1]
    return holes

















#def polyDefinition(,polyInstance,filePath = None,elementSize = 2e-5,writeFile=True):
#    '''    
#    Parameters
#    ----------
#    filePath : string
#    The path of the .dxf file to read.
#    	 elementSize : float
#    		The approximate target element size. NOTE: This only dictates the size
#    		of the straight line segments that approximate any curves. Actual mesh
#    		size is passed directly to triangle (for example, using the -a command 
#    		line switch)
#    
#    Currently compatible DXF entities are:
#    LWPOLYLINE: Treated as external/internal boundaries
#    CIRCLE: Treated as external/internal boundaries
#    ARC: Treated as external/internal boundaries
#    LINE: Treated as external/internal boundaries
#    POINT: Used as the centre for non-meshed areas (e.g. holes). 
#    Triangle will delete the medh from the point to the nearest boundary
#    
#    External boundaries are assumed to be one entitity; currently doesn't 
#    support combined boundaries for the external surface    
#    
#    This is not yet fully tested
#    '''        
#    if not filePath:
#        print('Creating empty graph...\n')
#        createEmptyGraph()
#    elif filePath[-4:]=='.dxf': 
#        try:
#            dxfFile = dxf.readfile(filePath)
#        except IOError:
#            print('No such DXF file. Maybe you meant .poly? Creating empty graph.')
#            createEmptyGraph()
#            return
#        entities = dxfFile.entities
#        pLines,isClosed = findPlines(entities,elementSize)
#        holes = findHoles(entities)
#        indexOfBoundary = findOuterBoundaryIndex(pLines)
#        numberOfVertices = sum(len(v) for v in pLines)
#        polylinesToPSLG(pLines,isClosed,indexOfBoundary)
#        
#        #writePoly2d(pLines,isClosed,holes,indexOfBoundary,filePath)
#    elif filePath[-5:]=='.poly':
#        try:
#            readPoly(filePath)
#        except IOError:
#            print('No such poly file. Maybe you meant .dxf? Creating empty graph.')
#            createEmptyGraph()
#    polyInstance.emptyPoly()
#    polyInstance.addVertices(vertices)
#    polyInstance.addEdges(edges)
#    polyInstance.addHoles(holes)
#    if writeFile:
#        writePoly2d(polyInstance,filePath)
#    
#
#def createEmptyGraph():
#    holes = np.array([])
#    vertices = np.empty([0,4])         
#    edges = np.empty([0,3])
#    numberOfVertices=0



#    def writePoly2d(,pslg,isClosed, holes, boundaryKey, filePath):
#       '''
#    	Function to write information loaded from the .dxf file to a .poly file.
#    	Parameters
#        ----------
#        pLines : list
#    		List of matplotlib.path object containing polyline information obtained from .dxf file.
#        isClosed : dict
#    		Dictionary of booleans corresponding to the pslg entities. 'true' forces the 
#    		pslg to close the entity, 'false' does not.
#    	boundaryKey : string
#    		String representing the pslg key corresponding to the external boundary
#    	filePath : string
#    		String indicating the destination to save the .poly file.
#            
#        Returns
#        -------
#    	0, regardless of success or failure. May implement error codes in the future.
#        '''
#        filePath = filePath[:-4]
#        filePath += '.poly'
#        filePath = filePath
#        with open(filePath,'w') as f:
#            f.write('<# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>\n')
#            numberOfVertices=sum(len(v) for v in pslg)
#            f.write("{} 2 ".format(numberOfVertices) + "0 {} \n".format(len(pslg[boundaryKey])))
#            
#            jj=1
#            for ii in range(len(pslg)):
#                
#                if ii==boundaryKey:
#                    writeStr = "{} {} 1\n"
#                else:
#                    writeStr = "{} {}\n"
#                for kk in range(len(pslg[ii].vertices)):
#                    f.write("{} ".format(jj)+writeStr.format(*pslg[ii].vertices[kk]))
#                    jj+=1
#            
#            f.write("{} 0\n".format(numberOfVertices - len(pslg) + sum(v for v in isClosed)))
#            
#            jj = 1
#            for ii in range(len(pslg)):
#                jStart = jj
#                for kk in range(len(pslg[ii])-1):
#                    f.write("{} {} {}\n".format(jj,jj,jj+1))
#                    jj+=1
#                if isClosed[ii]:
#                    f.write("{} {} {}\n".format(jj,jj,jStart))
#                    jj+=1
#            f.write("{}\n".format(len(holes)))
#            
#            for ii in range(len(holes)):
#                f.write("{} ".format(ii))
#                f.write("{} {}\n".format(*holes[ii]))
#        
#        return 0