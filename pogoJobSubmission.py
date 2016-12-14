# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 09:27:05 2016

@author: Alexander Ballisat

"""

import subprocess
import time
import os

def SubmitRun(input_name, path, abaqus=False, benchmark=0, cleanup=False, dimensions=2):
    '''
    Function to run Pogo jobs from Abaqus input files without having to go
    through the command line and does timings on the jobs as well. This also 
    does timing and removes intermediate files generated. This function is
    useful for doing parameter studies.
    
    Parameters
    ----------
    input_name : string
        The name of the Pogo job, it should be should be e.g. 'inputfile' NOT
        'inputfile.inp'. This is the name set as input_name in the
        WriteInput function.
        
    path : string
        The folder location where the files are located.
        
    abaqus : boolean, optional
        Whether an Abaqus input file is supplied which must be converted.
        Default is False.
        
    benchmark : boolean, optional
        Whether to return the timings, default is False which is not 
        
    cleanup : boolean, optional
        Whether to delete the intermediate files to save space or not.
        Defualt is False which does not delete anything, set to True to 
        delete the  .inp, .pogo-inp and .pogo-block files.
    
    dimensions : float, optional
        The number of dimensions of the problem, default is 2, the other
        option is 3.
        
        
    Returns
    -------
    convert_time : float, optional
        The time taken to convert the Abaqus input file into the Pogo format.
        
    block_time : float, optional
        The time taken to block the Pogo input file.
    
    solve_time : float, optional
        The time taken to solve the problem.
    
    
    '''
    
    if dimensions not in [2,3]:
        raise ValueError('Number of dimensions must be 2 or 3')
    
    t1 = time.clock()
    #input_name should be e.g. 'inputfile' NOT 'inputfile.inp'
    if abaqus == True:
        subprocess.call('pogoFromAb {}.inp'.format(input_name), cwd=path)
            
    if dimensions == 2:
        t2 = time.clock()
        if abaqus == True:
            print 'Convert time = {}'.format(t2-t1)
        
        subprocess.call('pogoBlock {}.pogo-inp'.format(input_name), cwd=path)
        t3 = time.clock()
        print 'Block time = {}'.format(t3-t2)

        subprocess.call('pogoSolve {} -o'.format(input_name), cwd=path)
        t4 = time.clock()
        print 'Solve time = {}'.format(t4-t3)
    
    if dimensions == 3:
        t2 = time.clock()
        if abaqus == True:
            print '\nConvert time = {}\n'.format(t2-t1)
            
        subprocess.call('pogoBlockGreedy3d {}.pogo-inp'.format(input_name), cwd=path)
        t3 = time.clock()
        print '\nBlock time = {}\n'.format(t3-t2)

        subprocess.call('pogoSolve3d {} -o'.format(input_name), cwd=path)
        t4 = time.clock()
        print 'Solve time = {}'.format(t4-t3)
        
    if cleanup == True:
        if abaqus == True:
            os.remove(r'{}\{}.inp'.format(path, input_name))
        os.remove(r'{}\{}.pogo-inp'.format(path, input_name))
        os.remove(r'{}\{}.pogo-block'.format(path, input_name))

    
    if benchmark != 0.:
        convert_time = t2-t1
        block_time = t3-t2
        solve_time = t4-t3
        
        return convert_time, block_time, solve_time