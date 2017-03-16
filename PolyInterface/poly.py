# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 10:12:28 2017

@author: Callum
"""
import numpy as np
import dxfgrabber as dxf
import polyService as pS
import pdb
class poly:
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
    def __init__(self,filePath = None,elementSize = 2e-5,writeFile=True):
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
            pLines,isClosed = pS.findPlines(entities,elementSize,precision=5)
            pLines,isClosed = pS.joinPlines(pLines,isClosed)
            holes = pS.findHoles(entities)
            indexOfBoundary=pS.findOuterBoundaryIndex(pLines)
            [vertices,boundaryFlags,edges]=pS.polylinesToPSLG(pLines,isClosed,indexOfBoundary)
        elif filePath[-5:]=='.poly':
            try:
                vertices,edges,holes=pS.readPoly(filePath)
            except IOError:
                print('No such poly file. Maybe you meant .dxf? Creating empty graph.')
                return
        self.addVertices(vertices,boundaryFlags)
        self.addEdges(edges)
        self.addHoles(holes)
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
    
    def emptyPoly(self):
        self.vertices = np.empty([0,2])
        self.holes = np.empty([0,3])
        self.edges = np.empty([0,2])
        self.numberOfVertices = 0
        self.numberOfEdges = 0
        self.numberOfHoles = 0
        self.numberOfBoundaryVertices = 0