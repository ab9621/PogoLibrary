# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 09:34:36 2016

@author: Alexander Ballisat

Version: 1.0
Date: 09/08/2016


This is a set of functions to generate an Abaqus input file in 2D that can
be converted into the Pogo format. The main function WriteInput is designed
so that it can just be run with the defaults with the exception of setting
the folder path.

Information on poly files and related mesh files can be found at 
https://www.cs.cmu.edu/~quake/triangle.html

Please send any bugs to alexander.ballisat@bristol.ac.uk.
"""

import numpy as np
import subprocess
import os

def node_locs(nodes):
    '''
    Function to convert the location of nodes passed as strings into floats.
    
    Parameters
    ----------
    nodes : sequence
        Sequence of the node strings which have the form [number, x, y, z,
        boundary marker]. This x and y coordinates are taken and put into an
        array.
        
    Returns
    -------
    locations : array, float
        An array of the node locations.
    '''
    return (np.array([a.split() for a in nodes[1:-2]])).astype(float)
    
def box_gen(box_x, box_y):
    '''
    Function to generate the coordinates of the corners of a box.
    
    Parameters
    ----------
    box_x : sequence, float
        The x coordinates of the corners of the box.
        
    box_y : sequence, float
        The y coordinates of the corners of the box.
        
    Returns
    -------
    locs : array, float
        The corner coordinates of the box.
    '''
    
    locs = np.zeros((4,2))
    locs[0] = [box_x[0], box_y[0]]
    locs[1] = [box_x[1], box_y[0]]
    locs[2] = [box_x[1], box_y[1]]
    locs[3] = [box_x[0], box_y[1]]
    return locs
    
def crack_gen(box_x, box_y, crack_origin, crack_length, crack_angle=0.0, crack_width=1E-5):
    '''
    Function to generate a box with a triangular crack emanating from the
    bottom face.
    
    Parameters
    ----------
    box_x : sequence, float
        The x coordinates of the corners of the box.
        
    box_y : sequence, float
        The y coordinates of the corners of the box.
        
    crack_origin : float
        The x coordinate of the centre of crack on the bottom face.
        
    crack_length : float
        The length of the crack. This must be less than the height of the box
        otherwise the shape will self-intersect.
        
    crack_angle : float
        The angle of the crack relative to the positive y_axis in radianns.
        Must be in the range (-pi/2, pi/2).
        
    crack_width : float
        The width of the crack on the base of the box. This must be small
        enough such that the crack edges do not pass through the sides of
        the box.
        
    Returns
    -------
    locs : array, float
        The coordinates of the points that make up the box with a crack in
        the base.
    '''
    # To do: add in more error checking.
    
    if box_y[1] - box_y[0] <= crack_length:
        raise ValueError('Crack is taller than the box.')

    if crack_angle < -0.5*np.pi or crack_angle > 0.5*np.pi:
        raise ValueError('Crack angle must be in the range (-pi/2, pi/2).')
        
    locs = np.zeros((7,2))
    locs[0] = [box_x[0], box_y[0]]
    locs[1] = [crack_origin[0]-crack_width/2., box_y[0]]    
    locs[2] = [crack_origin[0]+crack_length*np.sin(crack_angle), box_y[0]+crack_length*np.cos(crack_angle)]
    locs[3] = [crack_origin[0]+crack_width/2., box_y[0]]
    locs[4] = [box_x[1], box_y[0]]
    locs[5] = [box_x[1], box_y[1]]
    locs[6] = [box_x[0], box_y[1]]
    return locs

    
def find_nodes_at_loc(nodes, target, width):
    #Find and return the nodes at the target in a given width
    #Its to find the nodes in a transducer face for 2D
    #Currently assuming the transducer is on the top surface
    nodes = node_locs(nodes)
    locs = [0]*len(nodes) #arbitrary long list
    count = 0
    x_min = target[0] - width/2.
    x_max = target[0] + width/2.
    for a in nodes:
        if a[2] == target[1]:
            if x_min <= a[1] <= x_max:
                locs[count] = int(a[0])
                count += 1
    return locs[:count]

def gauss_tone(t, tau, f0):
    return np.exp(-1.*np.power((t-tau)/(1/f0), 2)) * np.sin(2*np.pi*f0*(t-tau))

def angled_beam_gaussian_analytic(nodes, trans_nodes, theta, signal, velocity):
    nodes = node_locs(nodes)
    z = nodes[np.logical_or.reduce([nodes[:,0] == x for x in trans_nodes])]
    z = z[:,1]
    lhs = np.min(z)
    rhs = np.max(z)

    n_points = len(signal[:,0])
    n_nodes = len(trans_nodes)

    dt = signal[1,0]-signal[0,0]
    n_pad = int((rhs-lhs)*np.sin(theta)/velocity/dt) + 160
    
    n_total = n_points + n_pad
    
    time = np.linspace(0, n_total-1, n_total)*dt
    amps = np.zeros((n_total, n_nodes))
    tau = (z - lhs)*np.sin(theta)/velocity+5E-7 #shift the centre
    f0 = 5E6
    for c1 in range(0, n_nodes):
        amps[:,c1] = gauss_tone(np.copy(time), tau[c1], f0)
    
    return time, amps, tau
    
def angled_beam(nodes, trans_nodes, theta, signal, velocity):
    nodes = node_locs(nodes)
    z = nodes[np.logical_or.reduce([nodes[:,0] == x for x in trans_nodes])]
    z = z[:,1]
    lhs = np.min(z)
    rhs = np.max(z)

    n_points = len(signal[:,0])
    n_nodes = len(trans_nodes)

    dt = signal[1,0]-signal[0,0]
    n_pad = int((rhs-lhs)*np.sin(theta)/velocity/dt) + 160
    
    trace = np.zeros(n_pad + n_points)
    if theta > 0.:
        trace[10:n_points+10] = signal[:,1]
    if theta < 0.:
        trace[n_pad-10:-10] = signal[:,1]
    fft = (np.fft.fft(trace))
    n_fft = int(len(fft)/2)
    fft = fft[:n_fft]
    
    fft_freq = np.fft.fftfreq(len(trace), dt)
    fft_freq = fft_freq[:n_fft]
    
    time = np.linspace(0, n_fft-1, n_fft)*dt
    amps = np.zeros((n_pad+n_points, n_nodes))
    
    dt_store = np.zeros(n_nodes)
    
    for c1 in range(0, n_nodes):
        f = np.copy(fft)
        t = (z[c1] - lhs)*np.sin(theta)/velocity
        f *= np.exp(-1j*2*np.pi*fft_freq*t)
        dt_store[c1] = t
        amps[:,c1] = np.real(np.fft.ifft(f, n_pad+n_points))
    
    return time, amps, dt_store
    
def WriteInput(material_name='Aluminium',
                v = 0.,
                density = 2700.,
                modulus = 7e+10,
                poisson = 0.34,
                angle = 0.,
                f = 5E6,
                probe_centre = [0.01, 0.03], #centre in global coordinates
                probe_width = 5e-3,
                amp = None,
                n_per_wavelength = 10.,
                max_area = 0.,
                element_type = 'CPE3',
                time_step = 3e-9,
                step_time_length = 5e-5,
                mesh_file_name = None,
                poly_file_name = None,
                x0 = 0.0,
                x1 = 0.08,
                y0 = 0.0,
                y1 = 0.03,
                crack = 0, #set to 1 for there to be a crack in the box
                crack_length = 0.01,
                crack_angle = 0.0,
                crack_origin = [0.04, 0.0],
                crack_width = 1E-6,
                path = 'C:/Users/useradmin/Desktop/Pogo',
                input_name = '',
                job_name = 'PogoModel',
                benchmark = 0,
                return_input_name = 0,
                field_output = 0):

    '''
    Function to write a 2D model input file in the Abaqus format with keyword
    arguments to allow it to be used in parameter studies. This currently
    only supports an object made of a single material. All units are SI by
    default, make sure they are all consistent.
    
    Parameters
    ----------
    material_name : string, optional
        The name of the material that is being modelled.
        
    density : float, optional
        The material density.
        
    modulus : float, optional
        The elastic modulus of the material.
        
    poisson : float, optional
        Poisson's ratio for the material.
        
    angle : float, optional
        The angle of the beam being generated in sample in radians. Default 
        is 0 in which case a signal may be passed in using the amp kwarg. 
        If it is at a non zero angle a 5 cycle Gaussian tone burst is
        generated as the necessary delays that have to be calculated
        introduce errors using FFTs.
        
    f : float, optional
        The centre frequency of the input spectrum.
        
    probe_centre : sequence, float, optional
        The coordinates of the centre of the probe. There is no test if this
        is at an appropriate place, this is worth checking if the model is
        not working.
        
    probe_width : float, optional
        The physical width of the probe. There is no test if this is larger
        than the sample.
    
    amp : string, optional
        A CSV file can be supplied as an input amplitude with the first
        column being the time and the second being the amplitude. There is
        no check to see if it is consistent. The default is None and a 5
        cycle tine burst with centre frequency f is generated.
        
    n_per_wavelength : float, optional
        The number of elements per wavelength. This is used to calculate the
        upper limit for the side of an element if a mesh is generated as part
        of this function. Default is 10.
        
    max_area : float, optional
        A constraint on the maximum area of an element to be used in the
        generation of a mesh. The default is 0 in which case the
        n_per_wavelength is used to generate a maximum area constraint.
        
    element_type : string, optional
        The Abaqus element type that is used. Currently only a single element
        type is supported in this function. The element type must be
        compatible with Pogo, see the Pogo documentation for supported
        element types.
        
    time_step : float, optional
        The time step for the FE calculation. Default is 3E-9.
        
    step_time_length : float, optional
        The total simulation time. Default is 5E-5.
        
    mesh_file_name : string, optional
        The name of the .ele and .node files that are generated using
        Triangle or any other mesher that can generate such files without the
        .1 on the end. This overrides all other mesh generation methods. 
        Default is None which generates a mesh based on other supplied 
        kwargs.
        
    poly_file_name : string, optional
        The name of a 2D poly file which defines the geometry which may be
        meshed. Default is None which generates a geometry and mesh based on
        other supplied kwargs.
        
    x0 : float, optional
        The lower x coordinate of the box that can be generated using this
        function.
        
    x1 : float, optional
        The upper x coordinate of the box that can be generated using this
        function.
        
    y0 : float, optional
        The lower y coordinate of the box that can be generated using this
        function.
        
    y1 : float, optional
        The upper y coordinate of the box that can be generated using this
        function.
    
    crack : float, optional
        Whether a crack is put in the generated box. Default is 0 which does
        not generate a crack, set to 1 to generate a crack.
        
    crack_length : float, optional
        If a crack is generated, this is the crack length. There is no test
        to see if this is fits within the box. Default is 0.01.
    
    crack_angle : float, optional
        The angle in radians of the crack relative to the positive y axis.
        Default is 0. which is a vertical crack.
        
    crack_origin : sequence, float, optional
        The origin of the crack. Default is [0.04, 0.0] which works with the
        default box.
        
    crack_width : float, optional
        The width of the crack on the bottom face, defualt is 1E-6.
        
    path : string, optional
        The folder in which the input file is to be written and the mesh
        files either exist or will be written in. Default is 
        'C:/Users/useradmin/Desktop/Pogo' but does really need to be changed.
        
    input_name : string, optional
        The name of the input file to be written, default is None in which
        case an input file name is generated using 
        '{}_by_{}_{}_mesh_{:.2e}'.format((x1-x0), (y1-y0), material_name,
        max_area).

    job_name : string, optional
        The Abaqus job name, default is 'PogoModel' but it has no general use
        in the rest of this script.
        
    benchmark : float, optional
        Whether benchmarking is being done, in which case the number of
        elements is returned. Default is 0 when it is not returned, set to
        1 to turn on.
        
    return_input_name : float, optional
        Whether the generated input name is returned, useful if it is being
        used with other scripts. Default is 0 which is not returned, set to 1
        to return.
        
    field_output : float, optional
        Whether field output is requested, default is 0 which is not
        requested. This is not fully implemented as history output is
        normally most useful.
        
    Returns
    -------
    n_elements : int, optional
        The number of elements in the mesh if benchmark is set to 1.
        
    input_name : string, optional
        The input name if return_input_name is set to 1.
        
    '''
    
    #Put some checks in
    if len(probe_centre) < 2:
        raise ValueError('Number of dimensions of probe centre is less than 2.')
    
    #If no velocity is set calculate it
    if v == 0.:
        v = np.sqrt((modulus*(1-poisson))/(density*(1+poisson)*(1-2*poisson)))
    print 'longitudinal velocity = {}'.format(v)
    
    #Rough element length constraint if not set in kwargs
    wavelength = v/f
    
    #Set a maximum element area size if none is set
    if max_area == 0.:
        dx = wavelength/n_per_wavelength
        print 'dx = {}'.format(dx)
        
        #Work out the maximum allowed area for each triangle
        max_area = (dx*dx)/2.
        
    print 'max area = {}\n'.format(max_area)
    
    #Load in the amplitude spectrum
    if amp != None:
        amplitude = np.loadtxt(amp, delimiter=',')
        
    else:
        f_sample = 100E6 #sample frequency
        dt = 1./f_sample
        tau = 2.5*1/f*1.1 #put a shift in to make sure the amplitude is tiny at t=0
        t_needed = 2*tau
        n_total = int(t_needed/dt)
        time = np.linspace(0, n_total-1, n_total)*dt
        amplitude = gauss_tone(np.copy(time), tau, f)
        amplitude = (np.vstack((time, amplitude))).T
    
    ##### Generate the geometry and create the mesh
    if mesh_file_name == None and poly_file_name == None:
        #Corners of the box
        x = [x0, x1]
        y = [y0, y1]
        
        #Whether or not there is a crack
        if crack == 0:
            points = box_gen(x, y)
        else:
            points = crack_gen(x, y, crack_origin, crack_length, crack_angle, crack_width) 
        
        #Number of points defined in the geometry
        n_points = len(points)
        
        #Sort out the input and filenames
        if input_name == '':
            #Consider a better naming convention
            input_name = '{}_by_{}_{}_mesh_{:.2e}'.format((x1-x0), (y1-y0), material_name, max_area)
        
        #create the poly file for input into Triangle
        geometry = 'Geometry_{}'.format(input_name)
        
        #Write out the node numbers and locations
        count = 1
        f = open('{}.poly'.format(geometry), 'w')
        f.write('{} 2 0 {}\n'.format(n_points, n_points)) #n_points points, 2 dimensions, no attributes, no boundary markers
        for a in points:
            f.write('{}   {} {} 1\n'.format(count, a[0], a[1]))
            count += 1
        
        #Write out the number of segments and number of boundary markers
        f.write('{} 0\n'.format(n_points)) #number of segments, number of boundary markers
        for c1 in range(1, n_points+1):
            f.write('{}    {} {}\n'.format(c1, c1, ((c1)%n_points)+1))    
        f.write('0')
        
        f.close()
    
    if poly_file_name != None:
        geometry = poly_file_name
        
    if mesh_file_name == None:
        #Call Triangle to generate the tri mesh
        subprocess.call("triangle -p -a{:.15f} -j -q30 -C {}".format(max_area, geometry), cwd=path)
    
    if mesh_file_name != None:
        geometry = mesh_file_name
        
    ##### Write the input file
    
    #Open the input file and write the preamble
    i = open('{}.inp'.format(input_name), 'w')
    i.write('*Heading\n')
    i.write('** Job name: {} Model name: Model-1\n'.format(job_name))
    i.write('*Preprint, echo=NO, model=NO, history=NO, contact=NO\n')
    i.write('**\n')
    i.write('** PART INSTANCE: Part-1-1\n')
    
    #Write the node locations
    i.write('*Node\n')
    f = open('{}.1.node'.format(geometry), 'r')
    nodes = f.read().split('\n')
    n_lines = len(nodes)
    n_nodes = n_lines-3
    
    for c1 in range(1,n_lines-2): #1 and -2 get rid of the extra bits
        line = filter(None, (nodes[c1]).split(' '))
        i.write(', '.join(line[:3]) + '\n')
    f.close()
    
    #to save space delete the unwanted files
    os.remove('{}/{}.1.node'.format(path, geometry))
    
    #Write the element types
    i.write('*Element, type={}\n'.format(element_type))
    f = open('{}.1.ele'.format(geometry), 'r')
    elements = f.read().split('\n')
    n_lines = len(elements)
    n_elements = n_lines-3
    
    for c1 in range(1,n_lines-2): #1 and -2 get rid of the extra bits
        line = filter(None, (elements[c1]).split(' '))
        i.write(', '.join(line) + '\n')    
    f.close()
    
    #to save space delete the unwanted files
    os.remove('{}/{}.poly'.format(path, geometry))
    os.remove('{}/{}.1.ele'.format(path, geometry))
    os.remove('{}/{}.1.poly'.format(path, geometry))
    
    #Create the all nodes and all elements sets
    i.write('*Nset, nset=All_nodes, generate\n')
    i.write('1, {}, 1\n'.format(n_nodes))
    
    i.write('*Elset, elset=All_elements, generate\n')
    i.write('1, {}, 1\n'.format(n_elements))
    
    #Create the node set which corresponds to the transducer
    i.write('*Nset, nset=Transducer\n')
    transducer_nodes = find_nodes_at_loc(nodes, probe_centre, probe_width/np.cos(angle))
    for a in transducer_nodes:
        i.write('{},\n '.format(a))
    
    #Write the section definition
    i.write('*Solid Section, elset=All_elements, material={}\n'.format(material_name))
    i.write(',\n')
    
    #Write the amplitude definition

    if angle == 0.:
        i.write('*System\n')
        i.write('*Amplitude, name=Amp-1\n')
        for a in amplitude:
            i.write('{}\n'.format(', '.join(map(str, a))))
    
    if angle != 0.:
        angle = np.pi/4. #angle in rads
        time, amps, delay_times = angled_beam_gaussian_analytic(nodes, transducer_nodes, angle, amplitude, v)
        for c1 in range(0, len(transducer_nodes)):
            i.write('*Amplitude, name=Amp-{}\n'.format(transducer_nodes[c1]))
            for c2 in range(0, len(time)):
                i.write('{}, {}\n'.format(time[c2], amps[c2, c1]))
    
    #Write the material definition
    i.write('*Material, name={}\n'.format(material_name))
    i.write('*Density\n')
    i.write('{},\n'.format(density))
    i.write('*Elastic\n')
    i.write('{}, {}\n'.format(modulus, poisson))
    
    #Write the step
    i.write('*Step, name=Step-1, nlgeom=YES\n')
    i.write('*Dynamic, Explicit, direct user control\n')
    i.write('{}, {}\n'.format(time_step, step_time_length))
    
    #Write bulk viscosity
    i.write('*Bulk Viscosity\n')
    i.write('0.0, 0.0\n')
    
    #Loads
    if angle == 0.:
        i.write('*Cload, amplitude=Amp-1\n')
        i.write('Transducer, 2, -10\n')
    
    if angle != 0.:
        for a in transducer_nodes:
            i.write('*Cload, amplitude=Amp-{}\n'.format(a))
            i.write('{}, 2, -10\n'.format(a))
        
    #Outputs
    if field_output != 0: 
        i.write('*Output, field, variable=PRESELECT\n')
        
    i.write('*Output, history, frequency=2\n')
    i.write('*Node Output, nset=Transducer\n')
    i.write('U2,\n')
    
    #End of step
    i.write('*End Step\n')
    
    i.close()
    
    print 'Input file {} written'.format(input_name)
    
    if benchmark != 0:
        return n_elements
        
    if return_input_name != 0:
        return input_name
        
    return