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
from CustomLineString import CustomLineString
from collections import Iterable
DEBUG = False
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
            self.polyLines = pLines
            self.generatePSLG()
        elif filePath[-5:]=='.poly':
            try:
                vertices,edges,holes=pS.readPoly(filePath)
                self.addVertices(vertices,boundaryFlags)
                self.addEdges(edges)
                self.addHoles(holes)
            except IOError:
                print('No such poly file. Maybe you meant .dxf? Creating empty graph.')
                return
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
            
            self.regions.append(CustomLineString(np.array(region.exterior.coords)))
        self._setNumberOfRegions()    
    
    def emptyPoly(self):
        self.vertices = np.empty([0,2])
        self.holes = np.empty([0,3])
        self.edges = np.empty([0,2])
        self.numberOfVertices = 0
        self.numberOfEdges = 0
        self.numberOfHoles = 0
        self.numberOfBoundaryVertices = 0
        
    def generatePSLG(self):
        holes = self.holes
        polyLines = self.polyLines
        self.emptyPoly()
        polyLines = pS.joinPlines(polyLines,polyLines)
        self.polyLines = polyLines
        indexOfBoundary=pS.findOuterBoundaryIndex(polyLines)
        [vertices,boundaryFlags,edges]=pS.polylinesToPSLG(polyLines,indexOfBoundary)
        [faces,sLines] = pS.findFaces(polyLines)
        regions = pS.splitFaces(faces,sLines)
        self.addVertices(vertices,boundaryFlags)
        self.addEdges(edges)
        self.addHoles(holes)
        self.addRegions(regions)    
    def addArc(self,arcPoint1,arcPoint2,arcPoint3,arcType = '3pt',boundary = False):
        #http://paulbourke.net/geometry/circlesphere/
        
        if arcType.lower() == '3pt':
            x,y,startTheta,endTheta,radius = aC.ThreePointToCenterAndAngles(arcPoint1,arcPoint2,arcPoint3)
        elif arcType.lower() == 'sce':
            x,y,startTheta,endTheta,radius = aC.SCEToCentreAndAngles(arcPoint1,arcPoint2,arcPoint3)
        elif arcType.lower() == 'sca':
            x,y,startTheta,endTheta,radius = aC.SCAToCentreAndAngles(arcPoint1,arcPoint2,arcPoint3)
        elif arcType.lower() == 'scl':
            x,y,startTheta,endTheta,radius = aC.SCLToCentreAndAngles(arcPoint1,arcPoint2,arcPoint3)
            
        polyLine = pS.createArcPolyline(startTheta,endTheta,radius,[x,y],self.elementSize)
        polyLine = CustomLineString(pS.setPrecision(polyLine.vertices(),self.elementSize,precision=5))
        self.polyLines.append(polyLine)
        self.generatePSLG()