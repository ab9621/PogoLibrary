"""
@author: Callum White

This is a library of functions related to converting AutoCAD .dxf files to
 .poly files for use with the Triangle mesh generation program. 
https://www.cs.cmu.edu/~quake/triangle.html
https://www.cs.cmu.edu/~quake/triangle.poly.html
     
 The functions included are:
 dxfToPoly
 writePoly2d
 findOuterBoundaryKey
 findPSLG
 closeLine
 convertToPolyline
 findHoles
 
 08/02/2017 
"""
import dxfgrabber as dxf
import matplotlib.path as mplPath
import warnings
import numpy as np




def dxfToPoly(filePath=None,elementSize = 2e-5) :
    '''
    
    Function to convert .dxf file to .poly file, ready for use with Triangle.
	
	Parameters
    ----------
    filePath : string
        The path of the .dxf file to read.
	elementSize : float
		The approximate target element size. NOTE: This only dictates the size
		of the straight line segments that approximate any curves. Actual mesh
		size is passed directly to triangle (for example, using the -a command 
		line switch)
        
    Returns
    -------
	0, regardless of success or failure. May implement error codes in the future.
    
    Currently compatible DXF entities are:
    LWPOLYLINE: Treated as external/internal boundaries
    CIRCLE: Treated as external/internal boundaries
    ARC: Treated as external/internal boundaries
    LINE: Treated as external/internal boundaries
    POINT: Used as the centre for non-meshed areas (e.g. holes). 
        Triangle will delete the medh from the point to the nearest boundary
        
    External boundaries are assumed to be one entitity; currently doesn't 
        support combined boundaries for the external surface    
    
    This is not yet fully tested
    '''
    

    dxfFile = dxf.readfile(filePath)
    entities = dxfFile.entities
    boundaries,isClosed = findPSLG(entities,elementSize)
    holes = findHoles(entities)
    key = findOuterBoundaryKey(boundaries)
    writePoly2d(boundaries,isClosed,holes,key,filePath)
    return 0

def writePoly2d(pslg,isClosed, holes, boundaryKey, filePath):
    '''
	Function to write information loaded from the .dxf file to a .poly file.
	Parameters
    ----------
    pslg : dict
		Dictionary of .dxf entities, represented as (x,y) pair arrays constituting a polyline.
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
    filePath = filePath[:-4]
    filePath += '.poly'
    
    with open(filePath,'w') as f:
        f.write('<# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>\n')
        numberOfVertices=sum(len(v) for v in pslg.values())
        f.write("{} 2 ".format(numberOfVertices) + "0 {} \n".format(len(pslg[boundaryKey])))
        boundaryNames = list(pslg.keys())
        jj=1
        for boundaryName in boundaryNames:
            if boundaryName==boundaryKey:
                writeStr = "{} {} 1\n"
            else:
                writeStr = "{} {}\n"
            for ii in range(len(pslg[boundaryName])):
                f.write("{} ".format(jj)+writeStr.format(*pslg[boundaryName].vertices[ii]))
                jj+=1
        
        f.write("{} 0\n".format(numberOfVertices - len(boundaryNames) + sum(v for v in isClosed.values())))
        
        jj = 1
        for boundaryName in boundaryNames:
            jStart = jj
            for ii in range(len(pslg[boundaryName])-1):
                f.write("{} {} {}\n".format(jj,jj,jj+1))
                jj+=1
            if isClosed[boundaryName]:
                f.write("{} {} {}\n".format(jj,jj,jStart))
                jj+=1
        f.write("{}\n".format(len(holes)))
        
        for ii in range(len(holes)):
            f.write("{} ".format(ii))
            f.write("{} {}\n".format(*holes[ii]))
    
    return 0
    
def findOuterBoundaryKey(pslg):
    '''
	Function to automatically determine the external boundary in the PSLG. 
	Essentially checks for any of the entities which encloses all other entities.
	
	Parameters
    -------
	pslg : dict
		Dictionary of .dxf entities, represented as (x,y) pair arrays constituting a polyline.
	
	Returns
    -------
    key : String
        The key corresponding to the enclosing boundary in the PSLG
    '''
    keys = list(pslg.keys())
    for ii in range(len(keys)):
        jj = 0
        key = keys[ii]
        print(key)
        while np.all(pslg[key].contains_points(pslg[keys[jj]].vertices)) or jj==ii:
            if jj == len(keys)-1:
                return key
            jj+=1
    errStr = "Must have one entity (LWPOLYLINE,CIRCLE) defining the outer boundary."
    raise ValueError(errStr)
    
def findPSLG(entities,elementSize):
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
    pslg : dict
		Dictionary of .dxf entities, represented as (x,y) pair arrays constituting a polyline.
	isClosed : dict
		Dictionary of booleans corresponding to the pslg entities. 'true' forces the 
		pslg to close the entity, 'false' does not.		
    '''
    numberOfEntities = len(entities)
    pslg = {}
    isClosed = {}
    jj=0
    for ii in range(numberOfEntities):
        boundaryName = "{}".format(entities[ii].dxftype) + ' {}'.format(jj)
        if entities[ii].dxftype == 'LWPOLYLINE':
            if entities[ii].is_closed:
                isClosed[boundaryName] = True
            else:
                isClosed[boundaryName] = False
            vertices = np.asarray(entities[ii].points)
                
            pslg[boundaryName] = mplPath.Path(vertices);
            jj = jj+1
        elif entities[ii].dxftype == 'CIRCLE' or entities[ii].dxftype == 'ARC':
            pslg[boundaryName]=convertToPolyline(entities[ii],elementSize)
            if entities[ii].dxftype == 'CIRCLE':
                isClosed[boundaryName] = True
            else:
                isClosed[boundaryName] = False
            jj = jj+1
        elif entities[ii].dxftype == 'LINE':
            pslg[boundaryName] = mplPath.Path(np.asarray([[entities[ii].start[0],entities[ii].start[1]],
                                                   [entities[ii].end[0]  ,entities[ii].end[1] ]]))
            isClosed[boundaryName] = False
            jj = jj+1
        elif entities[ii].dxftype == 'POINT':
            continue
        else:
            warnStr="Currently only LWPOLYLINE,CIRCLE,ARC and LINE supported for boundary definition. Found {}".format(entities[ii].dxftype)
            warnings.warn(warnStr)
    successStr = "Found {} covnertible entities.".format(jj)
    print(successStr)
    return (pslg,isClosed)

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
    elif entity.dxftype == 'ARC':
        thetaBegin = entity.start_angle * np.pi / 180
        thetaEnd = entity.end_angle * np.pi / 180
    angleRange = thetaEnd - thetaBegin            
    r = entity.radius
    dTheta = np.arcsin(elementSize/(2 * r)) * 2
    nTheta = angleRange/dTheta
    polyline = np.zeros([nTheta-1,2])
    theta = np.linspace(thetaBegin,thetaEnd,nTheta)
    polyline[:,0] = r * np.cos(theta[:-1]) + entity.center[0]
    polyline[:,1] = r * np.sin(theta[:-1]) + entity.center[1]
        
    polyline = mplPath.Path(polyline)
    return polyline
    
def findHoles(entities):
    numberOfEntities = len(entities)
    holes= np.zeros([0,2])
    for ii in range(numberOfEntities):
        if entities[ii].dxftype == 'POINT':
            holes = np.vstack((holes,entities[ii].point[0:2]))
    return holes