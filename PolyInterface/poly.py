# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 10:12:28 2017

@author: Callum
"""
import numpy as np
import dxfgrabber as dxf
import polyService as pS
import pdb
import matplotlib.path as mplPath
import arcConversions as aC
class Poly:
    vertices = np.empty([0,2])
    vertexIDs = np.empty([0,1])
    boundaryFlags = np.empty([0,1],dtype='int32')
    holes = np.empty([0,3])
    holeIDs = np.empty([0,1])
    edges = np.empty([0,2]).astype(int)
    edgeIDs = np.empty([0,1])
    numberOfVertices = 0
    numberOfEdges = 0    
    numberOfHoles = 0
    numberOfBoundaryVertices = 0
    elementSize = 0
    def __init__(self,filePath = None,elementSize = 2e-5,writeFile=True,addLines = False):
        '''    
        Parameters
        ----------
        filePath : string
        The path of the .dxf file to read.
        	 elementSize : float
        		The approximate target element size. NOTE: This only dictates the size
        		of the straight line segments that approximate any curves. Actual mesh
        		size is passed directly to triangle (for example, using the -a command 
        		line switch)
        
        Currently compatible DXF entities are:
        LWPOLYLINE: Treated as external/internal boundaries
        CIRCLE: Treated as external/internal boundaries
        ARC: Treated as external/internal boundaries
        LINE: Treated as external/internal boundaries
        POINT: Used as the centre for non-meshed areas (e.g. holes). 
        Triangle will delete the mesh from the point to the nearest boundary
        
        External boundaries are assumed to be one entitity; currently doesn't 
        support combined boundaries for the external surface    
        
        ''' 
        self.emptyPoly()
        self.elementSize = elementSize
        if not filePath:
            print('Creating empty graph...\n')
            return
        elif filePath[-4:]=='.dxf': 
            try:
                dxfFile = dxf.readfile(filePath)
            except IOError:
                print('No such DXF file. Maybe you meant .poly? Creating empty graph.')
                return
            entities = dxfFile.entities
            pLines,isClosed = pS.findPlines(entities,self.elementSize,precision=5)
            pLines,isClosed,_ = pS.joinPlines(pLines,pLines,isClosed,ignoreList = [])
            
            holes = pS.findHoles(entities)
            indexOfBoundary=pS.findOuterBoundaryIndex(pLines)
            [vertices,boundaryFlags,edges]=pS.polylinesToPSLG(pLines,isClosed,indexOfBoundary)
            [faces,sLines] = pS.findFaces(pLines,isClosed)
            regions = pS.splitFaces(faces,sLines)
        elif filePath[-5:]=='.poly':
            try:
                vertices,edges,holes=pS.readPoly(filePath)
            except IOError:
                print('No such poly file. Maybe you meant .dxf? Creating empty graph.')
                return
        
        self.addVertices(vertices,boundaryFlags)
        self.addEdges(edges)
        self.addHoles(holes)
        self.addRegions(regions)
        if writeFile:
            pS.writePoly2d(self,filePath)
            
        
    def _setNumberOfEdges(self):
        self.numberOfEdges = len(self.edges)
    
    def _setNumberOfVertices(self):
        self.numberOfVertices = len(self.vertices)
        
    def _setNumberOfHoles(self):
        self.numberOfHoles = len(self.holes)
    
    def _setNumberOfBoundaryVertices(self):
        self.numberOfBoundaryVertices = np.count_nonzero(self.vertices[:,-1])
        
    def _setVertexIDs(self):
        self.vertexIDs=np.arange(len(self.vertices)).astype(int) + 1
        
    def _setEdgeIDs(self):
        self.edgeIDs=np.arange(len(self.edges)).astype(int) + 1

    def _setHoleIDs(self):
        self.holeIDs=np.arange(len(self.holes)).astype(int) + 1                                
                                 
    def _setNumberOfRegions(self):
        self.numberOfRegions = len(self.regions)
        
    def addVertices(self,inputVertices,boundaryFlags):
        newVertices = np.zeros([len(self.vertices)+len(inputVertices),2])
        newVertices[:-len(inputVertices),:] = self.vertices
        newVertices[-len(inputVertices):,:] = inputVertices
        self.vertices = newVertices
        newBoundaryFlags = np.zeros([len(self.boundaryFlags)+len(inputVertices),1],dtype='int32')
        newBoundaryFlags[:-len(inputVertices)] = self.boundaryFlags
        newBoundaryFlags[-len(inputVertices):] = boundaryFlags
        self.boundaryFlags = newBoundaryFlags
        self._setNumberOfVertices()
        self._setNumberOfBoundaryVertices()
        self._setVertexIDs()
    
        
    def addEdges(self,inputEdges):
        newEdges = np.zeros([len(self.edges)+len(inputEdges),2])
        newEdges[:-len(inputEdges),:] = self.edges
        newEdges[-len(inputEdges):,:] = inputEdges
        self.edges = newEdges.astype(int)
        self._setNumberOfEdges()
        self._setEdgeIDs()
        
    def addHoles(self,inputHoles):
        newHoles = np.zeros([len(self.holes)+len(inputHoles),3])
        newHoles[:-len(inputHoles),:] = self.holes
        newHoles[-len(inputHoles):,:] = inputHoles
        self.holes = newHoles
        self._setNumberOfHoles()
        self._setHoleIDs()
    
    def addRegions(self,regions):
        self.regions = []
        for ii in range(len(regions)):
            region = regions[ii]
            
            self.regions.append(mplPath.Path(np.array(region.exterior.coords)))
        self._setNumberOfRegions()    
    
    def emptyPoly(self):
        self.vertices = np.empty([0,2])
        self.holes = np.empty([0,3])
        self.edges = np.empty([0,2])
        self.numberOfVertices = 0
        self.numberOfEdges = 0
        self.numberOfHoles = 0
        self.numberOfBoundaryVertices = 0
        
    def addArc(self,arcPoint1,arcPoint2,arcPoint3,type = '3pt',boundary = False):
        #http://paulbourke.net/geometry/circlesphere/
        # 
        #     _do3PtArc()
        #if type.lower() == 'sce':
        #     _doSCEArc()
        # elif type.lower() == 'sca'
        #     _doSCAArc()
        # elif type.lower() == 'scl'
        #     _doSCLArc()
        if type.lower() == '3pt':
            x,y,startTheta,endTheta = ac.ThreePointToCenterAndAngles(arcPoint1,arcPoint2,arcPoint3)
        elif type.lower() == 'sce':
            x,y,startTheta,endTheta = ac.SCEToCentreAndAngles(arcPoint1,arcPoint2,arcPoint3)
        elif type.lower() == 'sca':
            x,y,startTheta,endTheta = ac.SCAToCentreAndAngles(arcPoint1,arcPoint2,arcPoint3)
        elif type.lower() == 'scl':
            x,y,startTheta,endTheta = ac.SCLToCentreAndAngles(arcPoint1,arcPoint2,arcPoint3)
        if not boundary:
            indexOfBoundary = -1
        else:
            indexOfBoundary = 0
        r = np.sqrt((x-p1[0])**2+(y-p1[1])**2)
        polyline = pS.createArcPolyline(startTheta,endTheta,r,[x,y],self.elementSize)
        pointsNotInExisting = np.logical_not(mplPath.Path(self.vertices).contains_points(polyline.vertices))
        isNewBoundary = np.any(pointsNotInExisting)
        if isNewBoundary
        pdb.set_trace()
        [vertices,boundaryFlags,edges] = pS.polylinesToPSLG([polyline,],[False,],indexOfBoundary)
        self.addVertices(vertices,boundaryFlags)
        self.addEdges(edges)
        
def ThreePointToCenterAndAngles(p1,p2,p3):
    def _getDiffs():
        dy1 = float(p2[1]-p1[1])
        dx1 = float(p2[0]-p1[0])
        dy2 = float(p3[1]-p2[1])
        dx2 = float(p3[0]-p2[0])
        return dy1,dx1,dy2,dx2
    dy1,dx1,dy2,dx2 = _geqtDiffs()
    if dx1 == 0:
        p2,p3 = p3,p2
        dy1,dx1,dy2,dx2 = _getDiffs()
    elif dx2 == 0:
        p1,p2 = p2,p1
        dy1,dx1,dy2,dx2 = _getDiffs()
        
    grad1 = dy1/dx1
    grad2 = dy2/dx2
    denom = 2*(grad2-grad1)
    if denom == 0 or np.any([dx1,dx2]==0):
        errMsg = 'Points lie on parallel lines. No circle could be found'
        raise ValueError(errMsg)
    numerator = grad1*grad2*(p1[1]-p3[1])+grad2*(p1[0]+p2[0])-grad1*(p2[0]+p3[0])
    x = numerator / denom
    y = -1/grad1*(x-(p1[0]+p2[0])/2)+(p1[1]+p2[1])/2
    
    
    th1,th2,th3 = getAnglesRelativeToCenter(p1,[x,y]),getAnglesRelativeToCenter(p2,[x,y]),getAnglesRelativeToCenter(p3,[x,y])
    startTheta = min([th1,th2,th3])
    endTheta = max([th1,th2,th3])
    return x,y,startTheta,endTheta

def SCEToCentreAndAngles(p1,p2,p3):
    startTheta,endTheta = getAnglesRelativeToCentre(p1,p2),getAnglesRelativeToCentre(p3,p2)
    return p2[0],p2[1],startTheta,endTheta

def SCEToCentreAndAngles(p1,p2,p3):
    startTheta = getAnglesRelativeToCentre(p1,p2)
    return p2[0],p2[1],startTheta,p3
    
def SCLToCentreAndAngles(p1,p2,p3):
    startTheta = getAnglesRelativeToCentre(p1,p2)
    return p2[0],p2[1],startTheta,p3

def getAnglesRelativeToCenter(point,centre):
    dy = point[1]-centre[1]    
    dx = point[0]-centre[0]
    th = np.arctan2(dy,dx)
    return th