# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 16:45:00 2017

@author: Callum
"""
import matplotlib
#matplotlib.use('pdf')
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.collections import LineCollection
from matplotlib.widgets import Lasso
from matplotlib import colors as mcolors
import pdb
class PolyGui:
    def __init__(self,polyInstance,nodes=False,**kwargs):
        self.normalSelectedColor = np.array([[0, 0, 1, 1.0], [1, 0, 0, 1.0]])
        self.plotPoly(polyInstance,nodes,**kwargs)
        
    def plotPoly(self,polyInstance,nodes,**kwargs):
        plt.ion()
        lineList = []
        for ii in range(polyInstance.numberOfEdges):
            x0=polyInstance.vertices[polyInstance.edges[ii,0]-1,0]
            y0=polyInstance.vertices[polyInstance.edges[ii,0]-1,1]
            x1=polyInstance.vertices[polyInstance.edges[ii,1]-1,0]
            y1=polyInstance.vertices[polyInstance.edges[ii,1]-1,1]
            line = ((x0,y0),(x1,y1))
            lineList.append(line)
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111,**kwargs)
        if np.any(nodes):
            plt.plot(nodes[:,0],nodes[:,1],'.',zorder=1,MarkerSize=3)
        colors = [mcolors.to_rgba(c) for c in plt.rcParams['axes.prop_cycle'].by_key()['color']]
        self.lineData = LineCollection(tuple(lineList),pickradius=5,colors='r',zorder=10)
        self.lineData.set_picker(True)
        self.ax.add_collection(self.lineData)
        plt.autoscale(tight=True)
        plt.axis('equal')
        self.selected = np.zeros(len(self.lineData.get_segments())).astype(int)
        #self.ax.plot(polyInstance.vertices[:,1],polyInstance.vertices[:,2],'r.')
        #self.connect()
        #plt.savefig('FigureOutput')
        plt.show()
#    def connect(self):
#        self.picker = self.fig.canvas.mpl_connect('pick_event',self.onClick)
#    
#    def callback(self,verts)        
#    def onClick(self,event):
#        if event.inaxes is None:
#            return
#        self.lasso = Lasso(event.inaxes,event.xdata,event.ydata,self.callback)
#        self.fig.widgetlock(self.lasso)
#        print('testing')
#        ind = event.ind[0]
#        self.selected[ind] = 1 - self.selected[ind]
#        self.lineData.set_color(self.normalSelectedColor[self.selected])
#        self.fig.canvas.draw_idle()
#        thisline = event.artist
#        pdb.set_trace()
#        xdata = thisline.get_xdata()
#        ydata = thisline.get_ydata()
#        ind = event.ind
#        self.points = tuple(zip(xdata[ind], ydata[ind]))
#        print('onpick points:',self.points)
#        
#
#class Line(object):
#    colorin = mcolors.to_rgba("red")
#    colorout = mcolors.to_rgba("blue")
#
#    def __init__(self, x, y, include=False):
#        self.x = x
#        self.y = y
#        if include:
#            self.color = self.colorin
#        else:
#            self.color = self.colorout
