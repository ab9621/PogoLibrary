# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 16:45:00 2017

@author: Callum
"""
import matplotlib.pyplot as plt
import numpy as np
def readPoly(filePath):
    fileContents = open(filePath)    
    fileLines=fileContents.readlines()
    vertexInformation = fileLines[1].split()
    numberOfVertices = int(vertexInformation[0])
    lineInformation = fileLines[2+numberOfVertices].split()
    numberOfLines = int(lineInformation[0])
    vertices = np.zeros([numberOfVertices,2])
    lines = np.zeros([numberOfLines,2])
    for ii in range(numberOfVertices):
        vertString = fileLines[ii+2].split()
        vertices[ii,0] = float(vertString[1])
        vertices[ii,1] = float(vertString[2])
    
    for ii in range(numberOfLines):
        lineString = fileLines[ii+3+numberOfVertices].split()
        lines[ii,0] = int(lineString[1])-1
        lines[ii,1] = int(lineString[2])-1
    fileContents.close()
    lines = lines.astype(int)
    return vertices,lines
        
        
        
def onClick(event):
    print('testing')
    thisline = event.artist
    xdata = thisline.get_xdata()
    ydata = thisline.get_ydata()
    ind = event.ind
    points = tuple(zip(xdata[ind], ydata[ind]))
    print('onpick points:', points)
        
plotLines = []
filePath = r'D:\Callum\test.poly'    
fig = plt.figure()
ax = fig.add_subplot(111)
vertices,lines = readPoly(filePath)
for ii in range(len(lines)):
    plotLines.append(ax.plot( [vertices[lines[ii,0],0],vertices[lines[ii,1],0]], \
                                        [vertices[lines[ii,0],1],vertices[lines[ii,1],1]]
                                        ,'b-'
                                        ))
ax.plot(vertices[:,0],vertices[:,1],'r.')

picker = fig.canvas.mpl_connect('pick_event',onClick)
plt.show()

