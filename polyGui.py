# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 16:45:00 2017

@author: Callum
"""
import matplotlib.pyplot as plt
import numpy as np
import pdb
class polyGui:
    def __init__(self,polyFile):
        self.readPoly(polyFile)
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.plotLines = []
        self.plotPoly()
        
    
    
    
    def readPoly(self,filePath):
        fileContents = open(filePath)    
        fileLines=fileContents.readlines()
        vertexInformation = fileLines[1].split()
        numberOfVertices = int(vertexInformation[0])
        lineInformation = fileLines[2+numberOfVertices].split()
        numberOfLines = int(lineInformation[0])
        self.vertices = np.zeros([numberOfVertices,2])
        self.lines = np.zeros([numberOfLines,2])
        for ii in range(numberOfVertices):
            vertString = fileLines[ii+2].split()
            self.vertices[ii,0] = float(vertString[1])
            self.vertices[ii,1] = float(vertString[2])
        
        for ii in range(numberOfLines):
            lineString = fileLines[ii+3+numberOfVertices].split()
            self.lines[ii,0] = int(lineString[1])-1
            self.lines[ii,1] = int(lineString[2])-1
        fileContents.close()
        self.lines = self.lines.astype(int)
        
    def plotPoly(self):
        for ii in range(len(self.lines)):
            self.plotLines.append(self.ax.plot( [self.vertices[self.lines[ii,0],0],self.vertices[self.lines[ii,1],0]], \
                                                [self.vertices[self.lines[ii,0],1],self.vertices[self.lines[ii,1],1]], \
                                                'b-', \
                                                picker=5 \
                                                ))
        self.ax.plot(self.vertices[:,0],self.vertices[:,1],'r.')
        plt.show()
        self.connect()
        
    
    def connect(self):
        self.fig.canvas.mpl_connect('pick_event',self.onClick)
        
    def onClick(self,event):
        pdb.set_trace()
        print('testing')
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        self.points = tuple(zip(xdata[ind], ydata[ind]))
        print('onpick points:', self.points)
        