# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 16:45:00 2017

@author: Callum
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
import pdb
class polyGui:
    def __init__(self,polyInstance):
        self.normalSelectedColor = np.array([[0, 0, 1, 1.0], [1, 0, 0, 1.0]])
        self.plotPoly(polyInstance)
        
    def plotPoly(self,polyInstance):
        lineList = []    
        for ii in range(polyInstance.numberOfEdges):
            x0=polyInstance.vertices[polyInstance.edges[ii,1]-1,1]
            y0=polyInstance.vertices[polyInstance.edges[ii,1]-1,2]
            x1=polyInstance.vertices[polyInstance.edges[ii,2]-1,1]
            y1=polyInstance.vertices[polyInstance.edges[ii,2]-1,2]
            line = ((x0,y0),(x1,y1))
            lineList.append(line)
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.lineData = LineCollection(tuple(lineList),pickradius=5)
        self.lineData.set_picker(True)
        self.ax.add_collection(self.lineData)
        self.ax.set_xbound(lower=-.02,upper=.02)
        self.ax.set_ybound(lower=-.02,upper=.02)
        self.selected = np.zeros(len(self.lineData.get_segments())).astype(int)
        #self.ax.plot(polyInstance.vertices[:,1],polyInstance.vertices[:,2],'r.')
        self.connect()
        plt.show()
    
    def connect(self):
        self.picker = self.fig.canvas.mpl_connect('pick_event',self.onClick)
        
    def onClick(self,event):
        print('testing')
        ind = event.ind[0]
        self.selected[ind] = 1 - self.selected[ind]
        self.lineData.set_color(self.normalSelectedColor[self.selected])
        self.fig.canvas.draw_idle()
#        thisline = event.artist
#        pdb.set_trace()
#        xdata = thisline.get_xdata()
#        ydata = thisline.get_ydata()
#        ind = event.ind
#        self.points = tuple(zip(xdata[ind], ydata[ind]))
#        print('onpick points:',self.points)
#        