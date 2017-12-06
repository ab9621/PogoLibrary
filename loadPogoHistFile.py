# -*- coding: utf-8 -*-
"""
Created on Fri Oct 07 12:07:56 2016

@author: Alexander Ballisat

Script to load in a pogo history file and store it in a useful format.

07/10/2016 Version 1.0

Please email any bugs to alexander.ballisat@bristol.ac.uk
"""

import numpy as np
import struct
import matplotlib.pyplot as plt

class pogoHistData:
    '''
    Class to store the pogo history data. It is called by the loadPogoHistory
    function.
    
    Attributes
    ----------
    nodeNums : array, int
        The numbers of the nodes at which history measurements were taken.
        If multiple DoFs for a single node were taken this array will contain
        duplicates. Has length nMeasurements.
        
    nodeDoFs : array, int
        The degree of freedom for which the history measurement was taken for
        each node.
    
    nodePos : array, float
        The positions of the nodes at which history measurements are taken.
        Has the shape (nDims, nNodes). If multiple DoFs for a single node
        were taken this array will contain duplicates.
    
    histTraces : array, float
        The history traces recorded. Has shape (nMeasurements, nt).
    
    nt : int
        The number of time measurements taken for each trace.
    
    dt : float
        The time between sequential history measurements in seconds.
    
    startMeasure : int
        The time at which the first history measurement was taken.
    
    precision : string
        The precision the model was run at.
    
    '''
    
    def __init__(self,
                 nodeNums,
                 nodePos,
                 histTraces,
                 nt,
                 dt,
                 startMeasure,
                 precision,
                 nodeDof):
                     
         self.nodeNums = nodeNums
         self.nodePos = nodePos
         self.histTraces = histTraces
         self.nt = nt
         self.dt = dt
         self.startMeasure = startMeasure
         self.precision = precision
         self.timeBase = np.linspace(0, self.nt-1, self.nt) * self.dt
         self.nNodes = len(np.unique(nodeNums))
         self.nodeDofs = nodeDof
         
         return

def loadPogoHistory(fileName):
    '''
    Function to load a pogo history file and return the data. Supports
    versions 1.0 and 1.01.
    
    Parameters
    ----------
    fileName : string
        The name of the file to be loaded. Should include the file extension
        '.pogo-hist'
        
    Returns
    -------
    pogoHistoryData : instance of pogoHistData class
        The pogo history data in an instance of the pogoHistData class.
    
    '''
    if fileName[-10:] != '.pogo-hist':
        raise ValueError('File must be a .pogo-hist file.')

    with open(fileName, 'rb') as f:
        header = struct.unpack('20s', f.read(20))
        header = header[0].replace('\x00','')
        print header
        
        if header == '%pogo-hist1.0':
            fileVer = 1.0
        
        elif header == '%pogo-hist1.01':
            fileVer = 1.01
            
        else:
            raise ValueError('pogo history file version not recognised.')
        
        precision = struct.unpack('l', f.read(4))[0]
        
        if precision not in [4,8]:
            raise ValueError('Precision {} not supported. Should be 4 or 8.'.format(precision))
            
        if precision == 4:
            precString = 'f'
            
        elif precision == 8:
            precString = 'd'
            
        else:
            raise ValueError('Unsupported precision.')
            
        nDims = struct.unpack('l', f.read(4))[0]
        nMeasured = struct.unpack('l', f.read(4))[0]
        ntMeasured = struct.unpack('l', f.read(4))[0]
        dtMeasured = struct.unpack(precString, f.read(precision))[0]
        
        if fileVer >= 1.01:
            startMeasure = struct.unpack('l', f.read(4))[0] + 1
        else:
            startMeasure = 1
            
        nodeNums = np.zeros(nMeasured)
        nodeDofs = np.zeros(nMeasured)
        nodePos = np.zeros((nDims, nMeasured), dtype=precString)
        histTraces = np.zeros((nMeasured, ntMeasured), dtype=precString)
        
        for c1 in range(0, nMeasured):
            nodeNums[c1] = struct.unpack('l', f.read(4))[0] + 1
            nodeDofs[c1] = struct.unpack('l', f.read(4))[0] + 1
            nodePos[:, c1] = struct.unpack('{}{}'.format(nDims, precString), f.read(precision*nDims))
            histTraces[c1, :] = np.array(struct.unpack('{}{}'.format(ntMeasured, precString), f.read(precision*ntMeasured)), dtype=precString)
    
    histInstance = pogoHistData(nodeNums, nodePos, histTraces, ntMeasured, dtMeasured, startMeasure, precision, nodeDofs)
    
    return histInstance
    
#if __name__ == '__main__':
#    fileName = 'testPogoHist.pogo-hist'
#    ans = loadPogoHistory(fileName)