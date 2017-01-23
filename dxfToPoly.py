import dxfgrabber as dxf
import matplotlib.path as mplPath
import warnings
import numpy as np




def dxfToPoly(filePath=None,elementSize = 2e-5) :
    """
    dxfToPoly and associated functions
    by Callum White 15/12/2016
    
    This set of functions is aimed to enable easy conversion between AutoCAD DXF files and .poly files, 
    for use with the "triangle" meshing program.
    https://www.cs.cmu.edu/~quake/triangle.html
    https://www.cs.cmu.edu/~quake/triangle.poly.html
    
    Currently accepts file path and desired element size. 
    
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
    """
    

    dxfFile = dxf.readfile(filePath)
    entities = dxfFile.entities
    boundaries,isClosed = findBoundaries(entities,elementSize)
    holes = findHoles(entities)
    key = findOuterBoundaryKey(boundaries)
    writePoly2d(boundaries,isClosed,holes,key,filePath)
    return 0

def writePoly2d(boundaries,isClosed, holes, boundaryKey, filePath):
    filePath = filePath[:-4]
    filePath += '.poly'
    
    with open(filePath,'w') as f:
        f.write('# <Number of nodes> <Number of dimensions> <Number of Boundary markers>\n')
        numberOfVertices=sum(len(v) for v in boundaries.values())
        f.write("{} 2 ".format(numberOfVertices) + "{} \n".format(len(boundaries[boundaryKey])))
        boundaryNames = list(boundaries.keys())
        jj=1
        for boundaryName in boundaryNames:
            if boundaryName==boundaryKey:
                writeStr = "{} {} 1\n"
            else:
                writeStr = "{} {}\n"
            for ii in range(len(boundaries[boundaryName])):
                f.write("{} ".format(jj)+writeStr.format(*boundaries[boundaryName].vertices[ii]))
                jj+=1
        
        f.write("{} 0\n".format(numberOfVertices - len(boundaryNames) + sum(v for v in isClosed.values())))
        
        jj = 1
        for boundaryName in boundaryNames:
            jStart = jj
            for ii in range(len(boundaries[boundaryName])-1):
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
    
def findOuterBoundaryKey(boundaries):
    keys = list(boundaries.keys())
    for ii in range(len(keys)):
        jj = 0
        key = keys[ii]
        print(key)
        while np.all(boundaries[key].contains_points(boundaries[keys[jj]].vertices)) or jj==ii:
            if jj == len(keys)-1:
                return key
            jj+=1
    errStr = "Must have one entity (LWPOLYLINE,CIRCLE) defining the outer boundary."
    raise ValueError(errStr)
    
def findBoundaries(entities,elementSize):
    numberOfEntities = len(entities)
    boundaries = {}
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
                
            boundaries[boundaryName] = mplPath.Path(vertices);
            jj = jj+1
        elif entities[ii].dxftype == 'CIRCLE' or entities[ii].dxftype == 'ARC':
            boundaries[boundaryName]=convertToPolyline(entities[ii],elementSize)
            if entities[ii].dxftype == 'CIRCLE':
                isClosed[boundaryName] = True
            else:
                isClosed[boundaryName] = False
            jj = jj+1
        elif entities[ii].dxftype == 'LINE':
            boundaries[boundaryName] = mplPath.Path(np.asarray([[entities[ii].start[0],entities[ii].start[1]],
                                                   [entities[ii].end[0]  ,entities[ii].end[1] ]]))
            isClosed[boundaryName] = False
            jj = jj+1
        elif entities[ii].dxftype == 'POINT':
            continue
        else:
            warnStr="Currently only LWPOLYLINE,CIRCLE,ARC and LINE supported for boundary definition. Found {}".format(entities[ii].dxftype)
            warnings.warn(warnStr)
    successStr = "Found {} boundaries.".format(jj)
    print(successStr)
    return (boundaries,isClosed)

def closeLine(line):
    vertices = np.zeros([len(line)+1,2])
    vertices[0:len(line)-1,:] = line
    vertices[-1,:] = line[0,:]
    return vertices
    
def convertToPolyline(entity,elementSize)    :

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
    boundary = np.zeros([nTheta-1,2])
    theta = np.linspace(thetaBegin,thetaEnd,nTheta)
    boundary[:,0] = r * np.cos(theta[:-1]) + entity.center[0]
    boundary[:,1] = r * np.sin(theta[:-1]) + entity.center[1]
        
    boundary = mplPath.Path(boundary)
    return boundary
    
def findHoles(entities):
    numberOfEntities = len(entities)
    holes= np.zeros([0,2])
    for ii in range(numberOfEntities):
        if entities[ii].dxftype == 'POINT':
            holes = np.vstack((holes,entities[ii].point[0:2]))
    return holes