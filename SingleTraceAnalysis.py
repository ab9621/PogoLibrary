# -*- coding: utf-8 -*-
"""
Created on Mon May 23 10:38:04 2016

@author: ab9621
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as SS

def gauss_tone(t, tau, f0):
    return np.exp(-1.*np.power((t-tau)/(1/f0), 2)) * np.sin(2*np.pi*f0*(t-tau))

def Analyse(filename,
            trace=False,
            delimiter=',', 
            simple=False,
            element=0, 
            hilbert=False, 
            save=None, 
            f0=None, 
            mapPlot=False,
            returnFFT=False):
    '''
    Function to plot and analyse a single time trace by calculating the FFT.
    
    Input
    -----
    filename : string or array
        Name of the file to be loaded if a string or a preloaded data set to
        be analysed if an array is passed.
        
    trace : boolean
        Whether a trace has been passed as an array or not, default is false.
    
    delimiter : string
        String of the delimiter of the file, defaults to comma for csv files.
        
    simple : boolean
        Whether the information contained in the file is simple time vs
        amplitude in columns (True) or time and position information too
        (False). The default is set to False
        
    element : int
        For array results this is the element that you want to plot. Default is
        zero for single element probes.
    
    hilbert : Boolean
        Whether or not to plot the hilbert transform of the signal. Default is
        false
        
    save : string
        The filename if the plots are to be saved. Default is None which does
        not save the plots.
        
    f0 : float
        not currently used
        
    mapPlot : Boolean
        For array data this will plot all the traces in a 2D map if set to
        true. Default is False
        
    returnFFT : Boolean
        If True, return the FFT and the frequency spectrum used to create the
        plots. Default is false
    Output
    -----
    fft : array, complex float
        The fft of the signal
        
    fftfreq : array, float
        The frequencies present in the signal
        
    
    '''
    # Load in the file if its not a trace
    if trace == False:
        if filename.endswith('.csv'):
            raw = np.loadtxt(filename, delimiter=delimiter)
            
        elif filename.endswith('.npy'):
            raw = np.load(filename)
			
        else:
            raise ValueError('Input file type not recognised.')
            
    else:
        raw = filename
        
    # Check if the data is stored column by column
    dims = np.shape(raw)
    
    if dims[0] < dims[1]:
        raw = raw.T
        
    if dims[1] < 2:
        raise ValueError('Time information must be supplied')
    
    if simple == False:
        time = raw[:,-1]
        nCoords = len(np.where(time == 0)[0]) - 1
        rawTraces = raw[nCoords:, :-1]
        toAnalyse = rawTraces[:,element]
        time = time[nCoords:]
    
    elif simple == True:
        time = raw[:,0]
        toAnalyse = raw[:,1]
        
    dt = time[1] - time[0]
    nPoints = len(time)
    nPoints2 = nPoints/2
    
    fft = np.fft.fft(toAnalyse)
    fft = fft[:nPoints2]
    fftfreq = np.fft.fftfreq(nPoints, dt)
    fftfreq = fftfreq[:nPoints2]
    
    plt.clf()
    plt.figure(1)
    
    if mapPlot == False:
        
        plt.subplot(211)
        if hilbert == True:
            plt.plot(time*1E6, np.abs(SS.hilbert(toAnalyse)))
        else:
            plt.plot(time*1E6, toAnalyse, label='Raw Trace')
            
        plt.xlabel('Time ($\mu$s)')
        plt.ylabel('Amplitude (AU)')
        plt.legend()
        
        plt.subplot(212)
        plt.plot(fftfreq*1E-6, np.abs(fft), label='Frequency Spectrum')
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('Amplitude (AU)')
        plt.legend()
    
    if mapPlot == True and simple == False:
        nElements =  len(rawTraces[0])
        if nElements < 2:
            raise ValueError('Not array data.\n')
            return
        plt.subplot(211)
        if hilbert == False:
            plt.imshow(rawTraces.T, interpolation='None', aspect='auto', origin='lower', extent=[time[0]*1E6, time[-1]*1E6, 0, nElements])
        else:
            plt.imshow((np.abs(SS.hilbert(rawTraces,axis=1))).T, interpolation='None', aspect='auto', origin='lower', extent=[time[0]*1E6, time[-1]*1E6, 0, nElements])
        plt.xlabel('Time ($\mu$s)')
        plt.ylabel('Element Number')
        plt.colorbar()
        
        plt.subplot(212)
        fft2 = np.fft.fft(rawTraces, axis=1)
        fft2 = fft2[:nPoints2, :]
        plt.imshow(np.abs(fft2).T, interpolation='None', aspect='auto', origin='lower', extent=[fftfreq[0]*1E-6, fftfreq[-1]*1E-6, 0, len(rawTraces[0])])
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('Element Number')
            
    if save != None:
        plt.savefig(save, bbox_inches='tight')
        
    if returnFFT == True:
        return fft, fftfreq
        
    return