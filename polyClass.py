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
import dxfgrabber as dxf
import matplotlib.path as mplPath
import warnings
import numpy as np
import matplotlib.pyplot as plt
import pdb
class polyFile:
    def __init__(self,filePath = None,elementSize = 2e-5) :
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
            Triangle will delete the medh from the point to the nearest boundary
            
        External boundaries are assumed to be one entitity; currently doesn't 
            support combined boundaries for the external surface    
        
        This is not yet fully tested
        '''
        if not filePath:
            print('Creating empty graph...\n')
            self.createEmptyGraph()
        elif filePath[-4:]=='.dxf': 
            try:
                dxfFile = dxf.readfile(filePath)
            except IOError:
                print('No such DXF file. Maybe you meant .poly? Creating empty graph.')
                self.createEmptyGraph()
                return
            entities = dxfFile.entities
            pLines,isClosed = self.findPlines(entities,elementSize)
            self.holes = self.findHoles(entities)
            indexOfBoundary = self.findOuterBoundaryIndex(pLines)
            self.numberOfVertices = sum(len(v) for v in pLines)
            self.polylinesToPSLG(pLines,isClosed,indexOfBoundary)
            self.writePoly2d(pLines,isClosed,self.holes,indexOfBoundary,filePath)
        elif filePath[-5:]=='.poly':
            try:
                self.readPoly(filePath)
            except IOError:
                print('No such poly file. Maybe you meant .dxf? Creating empty graph.')
                self.createEmptyGraph()
                
        
    
    def createEmptyGraph(self):
        self.holes = np.array([])
        self.vertices = np.empty([0,4])         
        self.edges = np.empty([0,3])
        self.numberOfVertices=0
        
    def polylinesToPSLG(self,pLines,isClosed,indexOfBoundary):
        jj=1
        self.vertices = np.zeros([self.numberOfVertices,4])
        self.edges = []
        for ii in range(len(pLines)):            
            for kk in range(len(pLines[ii].vertices)):
                if ii==indexOfBoundary:
                    self.vertices[jj-1,:] = [jj,pLines[ii].vertices[kk][0],pLines[ii].vertices[kk][1],1]
                else:
                    self.vertices[jj-1,:] = [jj,pLines[ii].vertices[kk][0],pLines[ii].vertices[kk][1],0]
                jj+=1
        jj = 1
        for ii in range(len(pLines)):
            jStart = jj
            for kk in range(len(pLines[ii])-1):
                self.edges.append([jj,jj,jj+1])
                jj+=1
            if isClosed[ii]:
                self.edges.append([jj,jj,jStart])
                jj+=1
        self.edges = np.array(self.edges)
        
    def writePoly2d(self,pslg,isClosed, holes, boundaryKey, filePath):
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
        filePath = filePath[:-4]
        filePath += '.poly'
        self.filePath = filePath
        with open(filePath,'w') as f:
            f.write('<# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>\n')
            numberOfVertices=sum(len(v) for v in pslg)
            f.write("{} 2 ".format(numberOfVertices) + "0 {} \n".format(len(pslg[boundaryKey])))
            
            jj=1
            for ii in range(len(pslg)):
                
                if ii==boundaryKey:
                    writeStr = "{} {} 1\n"
                else:
                    writeStr = "{} {}\n"
                for kk in range(len(pslg[ii].vertices)):
                    f.write("{} ".format(jj)+writeStr.format(*pslg[ii].vertices[kk]))
                    jj+=1
            
            f.write("{} 0\n".format(numberOfVertices - len(pslg) + sum(v for v in isClosed)))
            
            jj = 1
            for ii in range(len(pslg)):
                jStart = jj
                for kk in range(len(pslg[ii])-1):
                    f.write("{} {} {}\n".format(jj,jj,jj+1))
                    jj+=1
                if isClosed[ii]:
                    f.write("{} {} {}\n".format(jj,jj,jStart))
                    jj+=1
            f.write("{}\n".format(len(holes)))
            
            for ii in range(len(holes)):
                f.write("{} ".format(ii))
                f.write("{} {}\n".format(*holes[ii]))
        
        return 0
        
    def findOuterBoundaryIndex(self,pLines):
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
        
    def findPlines(self,entities,elementSize):
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
        pslg = {}
        pLines = []
        isClosed = {}
        isClosedList = []
        jj=0
        for ii in range(numberOfEntities):
            boundaryName = "{}".format(entities[ii].dxftype) + ' {}'.format(jj)
            if entities[ii].dxftype == 'LWPOLYLINE':
                if entities[ii].is_closed:
                    isClosed[boundaryName] = True
                    isClosedList.append(True)
                else:
                    isClosed[boundaryName] = False
                    isClosedList.append(True)
                vertices = np.asarray(entities[ii].points)
                    
                pslg[boundaryName] = mplPath.Path(vertices)
                pLines.append(mplPath.Path(vertices))
                jj = jj+1
            elif entities[ii].dxftype == 'CIRCLE' or entities[ii].dxftype == 'ARC':
                pslg[boundaryName]=self.convertToPolyline(entities[ii],elementSize)
                pLines.append(self.convertToPolyline(entities[ii],elementSize))
                if entities[ii].dxftype == 'CIRCLE':
                    isClosed[boundaryName] = True
                    isClosedList.append(True)
                else:
                    isClosed[boundaryName] = False
                    isClosedList.append(True)
                jj = jj+1
            elif entities[ii].dxftype == 'LINE':
                pslg[boundaryName] = mplPath.Path(np.asarray([[entities[ii].start[0],entities[ii].start[1]],
                                                       [entities[ii].end[0]  ,entities[ii].end[1] ]]))
                pLines.append(mplPath.Path(np.asarray([[entities[ii].start[0],entities[ii].start[1]],
                                                       [entities[ii].end[0]  ,entities[ii].end[1] ]])))
                isClosed[boundaryName] = False
                isClosedList.append(False)
                jj = jj+1
            elif entities[ii].dxftype == 'POINT':
                continue
            else:
                warnStr="Currently only LWPOLYLINE,CIRCLE,ARC and LINE supported for boundary definition. Found {}".format(entities[ii].dxftype)
                warnings.warn(warnStr)
        successStr = "Found {} convertible entities.".format(jj)
        print(successStr)
        return (pLines,isClosedList)
    
    
#    def splineToPline(entity):
#        degreeOfSpline = entity.degree
#        controlPoints = entity.control_points
#        knots = entity.knots
        
        
    def closeLine(self,vertices):
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
        
    def convertToPolyline(self,entity,elementSize)    :
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
        else:
            errStr = 'Error: Invalid entity passed to convertToPolyline.'
            raise ValueError(errStr)
                
        angleRange = thetaEnd - thetaBegin            
        r = entity.radius
        dTheta = np.arcsin(elementSize/(2 * r)) * 2
        nTheta = int(angleRange/dTheta)
        polyline = np.zeros([nTheta-1,2])
        theta = np.linspace(thetaBegin,thetaEnd,nTheta)
        polyline[:,0] = r * np.cos(theta[:-1]) + entity.center[0]
        polyline[:,1] = r * np.sin(theta[:-1]) + entity.center[1]
            
        polyline = mplPath.Path(polyline)
        return polyline
        
    def findHoles(self,entities):
        numberOfEntities = len(entities)
        holes= np.zeros([0,2])
        for ii in range(numberOfEntities):
            if entities[ii].dxftype == 'POINT':
                holes = np.vstack((holes,entities[ii].point[0:2]))
        return holes
    
    def plotPoly(self,filePath=None):
        print('Plotting figure... (this might take a while for large graphs)')
        self.fig=plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.plotLines = []
        for ii in range(len(self.edges)):
            self.plotLines.append(self.ax.plot([self.vertices[self.edges[ii,1]-1,1],self.vertices[self.edges[ii,2]-1,1]], \
                                               [self.vertices[self.edges[ii,1]-1,2],self.vertices[self.edges[ii,2]-1,2]], \
                                               'b-'))
        self.ax.plot(self.vertices[:,1],self.vertices[:,2],'r.')
        plt.axis('equal')
#        fig.canvas.mpl_connect('pick_event',onpick)
        plt.show()
        
        
    def readPoly(self,filePath):
        fileContents = open(filePath,'r')    
        fileLines=fileContents.readlines()
        vertexInformation = fileLines[1].split()
        self.numberOfVertices = int(vertexInformation[0])
        lineInformation = fileLines[2+self.numberOfVertices].split()
        numberOfLines = int(lineInformation[0])
        self.vertices = np.zeros([self.numberOfVertices,4])
        self.edges = np.zeros([numberOfLines,3])
        
        for ii in range(self.numberOfVertices):
            vertString = fileLines[ii+2].split()
            self.vertices[ii,0] = float(vertString[0])
            self.vertices[ii,1] = float(vertString[1])
            self.vertices[ii,2] = float(vertString[2])
            if len(vertString)==3:
                self.vertices[ii,3] = 0
            else:
                self.vertices[ii,3] = float(vertString[3])
        
        for ii in range(numberOfLines):
            lineString = fileLines[ii+3+self.numberOfVertices].split()
            self.edges[ii,0] = int(lineString[0])
            self.edges[ii,1] = int(lineString[1])
            self.edges[ii,2] = int(lineString[2])
            fileContents.close()
        self.edges = self.edges.astype(int)
    
