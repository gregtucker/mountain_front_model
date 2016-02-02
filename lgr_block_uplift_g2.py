#!/usr/env/python
"""
Hillslope model with block uplift.

TODO NOTES:
- get a realistic set of parms to play with
- clean up the code to make it "externally runnable"
- make a master run script
- start a matrix of runs exploring different u, d, and L
- while that's going, do some profiling to find speed bottlenecks
"""

_DEBUG = False

import time
#import random
from numpy import zeros, bincount, arange, savetxt, sqrt, log10, mean, arctan, pi
from pylab import subplots, plot, show, xlabel, ylabel, title, axis, figure, clf, savefig
from landlab import HexModelGrid
from landlab.io.netcdf import write_netcdf
from landlab.ca.celllab_cts import Transition, CAPlotter
from landlab.ca.oriented_hex_cts import OrientedHexCTS
from landlab.ca.boundaries.hex_lattice_tectonicizer import LatticeUplifter
from scipy.optimize import curve_fit
import matplotlib

def add_weathering_and_disturbance_to_transition_list(xn_list, d=0.0, w=0.0):
    """
    Add transition rules representing weathering and/or grain disturbance
    to the list, and return the list.
    
    Parameters
    ----------
    xn_list : list of Transition objects
        List of objects that encode information about the link-state 
        transitions. Normally should first be initialized with lattice-grain
        transition rules, then passed to this function to add rules for
        weathering and disturbance.
    d : float (optional)
        Rate of transition (1/time) from fluid / resting grain pair to
        mobile-grain / fluid pair, representing grain disturbance.
    w : float (optional)
        Rate of transition (1/time) from fluid / rock pair to
        fluid / resting-grain pair, representing weathering.
    
    
    Returns
    -------
    xn_list : list of Transition objects
        Modified transition list.
    """
    
    # Disturbance rule
    xn_list.append( Transition((7,0,0), (0,1,0), d, 'disturbance') )
    xn_list.append( Transition((7,0,1), (0,2,1), d, 'disturbance') )
    xn_list.append( Transition((7,0,2), (0,3,2), d, 'disturbance') )
    xn_list.append( Transition((0,7,0), (4,0,0), d, 'disturbance') )
    xn_list.append( Transition((0,7,1), (5,0,1), d, 'disturbance') )
    xn_list.append( Transition((0,7,2), (6,0,2), d, 'disturbance') )

    if _DEBUG:
        print
        print 'setup_transition_list(): list has',len(xn_list),'transitions:'
        for t in xn_list:
            print '  From state',t.from_state,'to state',t.to_state,'at rate',t.rate,'called',t.name
        
    return xn_list

    
def get_profile_and_soil_thickness(grid, data):
    
    nr = grid.number_of_node_rows
    nc = grid.number_of_node_columns
    elev = zeros(nc)
    soil = zeros(nc)
    for c in range(nc):
        e = (c%2)/2.0
        s = 0
        r = 0
        while r<nr and data[c*nr+r]!=0:
            e+=1
            if data[c*nr+r]==7:
                s+=1
            r+=1
        elev[c] = e
        soil[c] = s
    return elev, soil

    
def flin(x, m, b):
    """
    Linear function for curve-fitting.
    """
    return m*x+b
    
    
def fquad(x, a, b, c):
    """
    Quadratic function for curve-fitting.
    """
    return a*x*x+b*x+c
    
    
def write_output(grid, outfilename, iteration):
    """Write output to file (currently netCDF)."""
    filename = outfilename + str(iteration).zfill(4) + '.nc'
    write_netcdf(filename, grid)

    
def run(uplift_interval, d): #d_ratio_exp):

    # INITIALIZE

    #uplift_interval = 1e7

    #filenm = 'test_output'
    #imagenm = 'Hill141213/hill'+str(int(d_ratio_exp))+'d'
    
    
    # Remember the clock time, and calculate when we next want to report
    # progress.
    current_real_time = time.time()
    next_report = current_real_time + report_interval
    next_uplift = uplift_interval

    # Create a grid
    hmg = HexModelGrid(nr, nc, 1.0, orientation='vertical', shape='rect',
                       reorient_links=True)

    # Close the right-hand grid boundaries
    #hmg.set_closed_nodes(arange((nc-1)*nr, hmg.number_of_nodes))

    # Set up the states and pair transitions.
    # Transition data here represent particles moving on a lattice: one state
    # per direction (for 6 directions), plus an empty state, a stationary
    # state, and a wall state.
    ns_dict = { 0 : 'empty', 
                1 : 'moving up',
                2 : 'moving right and up',
                3 : 'moving right and down',
                4 : 'moving down',
                5 : 'moving left and down',
                6 : 'moving left and up',
                7 : 'rest',
                8 : 'wall'}
    xn_list = setup_transition_list(g, f, d)
    #xn_list = []
    #xn_list.append( Transition((0,1,0), (0,7,0), g, 'gravity 1') )


    # Create data and initialize values.
    node_state_grid = hmg.add_zeros('node', 'node_state_grid', dtype=int)

    # Lower rows get resting particles
    if nc % 2 == 0:  # if even num cols, bottom right is...
        bottom_right = nc - 1
    else:
        bottom_right = nc // 2
    right_side_x = 0.866025403784 * (nc - 1)
    for i in range(hmg.number_of_nodes):
        if hmg.node_y[i] < 3.0:
            if hmg.node_x[i] > 0.0 and hmg.node_x[i] < right_side_x:
                node_state_grid[i] = 7
        #elif hmg.node_x[i]>((nc-1)*0.866):
        #    node_state_grid[i] = 8
    node_state_grid[0] = 8  # bottom left
    node_state_grid[bottom_right] = 8
    #for i in range(hmg.number_of_nodes):
    #    print i, hmg.node_x[i], hmg.node_y[i], node_state_grid[i]

    # Create an uplift object
    uplifter = LatticeUplifter(hmg, node_state_grid)

    # Create the CA model
    ca = OrientedHexCTS(hmg, ns_dict, xn_list, node_state_grid)

    # Create a CAPlotter object for handling screen display
    # potential colors: red3='#CD0000'
    #mob = 'r'
    #rock = '#5F594D'
    sed = '#A4874B'
    #sky = '#CBD5E1'
    #sky = '#85A5CC'
    sky = '#D0E4F2'
    rock = '#000000' #sky
    mob = '#D98859'
    #mob = '#DB764F'
    #mob = '#FFFF00'
    #sed = '#CAAE98'
    #clist = [(0.5, 0.9, 0.9),mob, mob, mob, mob, mob, mob,'#CD6839',(0.3,0.3,0.3)]
    clist = [sky,mob, mob, mob, mob, mob, mob,sed,rock]
    my_cmap = matplotlib.colors.ListedColormap(clist)
    ca_plotter = CAPlotter(ca, cmap=my_cmap)
    k=0

    # Plot the initial grid
    ca_plotter.update_plot()
    axis('off')
    #savefig(imagenm+str(k)+'.png')
    k+=1
    
    # Write output for initial grid
    #write_output(hmg, filenm, 0)
    #output_iteration = 1

    # Create an array to store the numbers of states at each plot interval
    #nstates = zeros((9, int(run_duration/plot_interval)))
    #k = 0

    # Work out the next times to plot and output
    next_output = output_interval
    next_plot = plot_interval

    # RUN
    current_time = 0.0
    while current_time < run_duration:
        
        # Figure out what time to run to this iteration
        next_pause = min(next_output, next_plot)
        next_pause = min(next_pause, next_uplift)
        next_pause = min(next_pause, run_duration)

        # Once in a while, print out simulation and real time to let the user
        # know that the sim is running ok
        current_real_time = time.time()
        if current_real_time >= next_report:
            print 'Current sim time',current_time,'(',100*current_time/run_duration,'%)'
            next_report = current_real_time + report_interval

        # Run the model forward in time until the next output step
        print('Running to...' + str(next_pause))
        ca.run(next_pause, ca.node_state) #, 
               #plot_each_transition=plot_every_transition, plotter=ca_plotter)
        current_time = next_pause
        
        # Handle output to file
        if current_time >= next_output:
            #write_output(hmg, filenm, output_iteration)
            #output_iteration += 1
            next_output += output_interval
            
        # Handle plotting on display
        if current_time >= next_plot:
            #node_state_grid[hmg.number_of_node_rows-1] = 8
            ca_plotter.update_plot()
            axis('off')
            next_plot += plot_interval

        # Handle uplift
        if current_time >= next_uplift:
            uplifter.uplift_interior_nodes(rock_state=7)
            ca.update_link_states_and_transitions(current_time)
            next_uplift += uplift_interval

    print('Finished with main loop')

    # FINALIZE

    # Plot
    #ca_plotter.update_plot()
    #ca_plotter.finalize()
    #axis('off')

    #(elev, soil) = get_profile_and_soil_thickness(hmg, node_state_grid)
    #print elev
    #print soil
    
    #print 'Mean thickness:',mean(soil[75:])
    
    #x = 0.5*sqrt(3)*arange(nc)
    
    # Calculate slope in upper half
    #(cf,cm1) = curve_fit(flin, x[75:], elev[75:])
    #slpgrad = cf[0]
    #slpdeg = 180.*arctan(slpgrad)/pi
    #print 'Upper slope gradient', slpgrad,'angle',slpdeg
    
    # Calculate (roughly) Laplacian
    #(pf,cm2) = curve_fit(fquad, x[75:], elev[75:] )
    #second_deriv = 2*pf[0]
    #print '2nd deriv',second_deriv
    
    #savetxt(filenm, (elev,soil))
    #savetxt(filenm+'-ns', node_state_grid)

    #figure(2)
    #clf()
    #plot(elev)
    #show()
    #figure(3)
    #plot(soil)
    #show()
    # Display the numbers of each state
    #fig, ax = subplots()
    #for i in range(1, 8):
    #    plot(arange(plot_interval, run_duration+plot_interval, plot_interval), nstates[i,:], label=ns_dict[i])
    #ax.legend()
    #xlabel('Time')
    #ylabel('Number of particles in state')
    #title('Particle distribution by state')
    #axis([0, run_duration, 0, 2*nstates[7,0]])
    #show()

    # Reports ups vs downs
    #n_up = nstates[1,:]+nstates[2,:]+nstates[6,:]
    #n_down = nstates[3,:]+nstates[4,:]+nstates[5,:]
    #print n_up
    #print n_down
    #print n_up/(n_up+n_down)
    #print n_down/(n_up+n_down)


def get_params():
    """Set and return the various model parameters."""

    params = {}
    params['num_rows'] = 113
    params['num_cols'] = 127
    params['grav_transition_rate'] = 1.0
    params['fric_transition_rate'] = 0.7
    params['plot_interval'] = 1.0e9
    params['output_interval'] = 1.0e99
    params['run_duration'] = 1.0e10
    params['report_interval'] = 5.0  # report interval, in real-time seconds
    params['plot_every_transition'] = False
    params['uplift_interval_yrs'] = 0.01
    params['uplift_interval'] = \
        params['uplift_interval_yrs'] / (3600 * 24 * 365.25)
    params['cell_size'] = 0.1
    params['disturbance_rate'] = 100.0 / (3600 * 24 * 365.25)

    return params
    
    
def report_run_params(params):
    """Display on screen some basic information about the run."""
    spy = 3600 * 24 * 365.25  # number of seconds in a year
    print('Parameters for this run:')
    print('  Baselevel rate: ' + \
        str(params['cell_size'] / params['uplift_interval_yrs']) + ' m/yr')
    print('  Disturbance rate: ' + \
        str(params['disturbance_rate'] * spy) + ' events/yr = ' + \
        str(params['disturbance_rate']) + ' events/s = ' + \
        str(params['cell_size'] * params['disturbance_rate'] * spy) + \
        ' m/yr = ' + \
        str(params['cell_size'] * params['disturbance_rate']) + 'm/s')
    print('  Dimensionless disturbance rate: ' + \
        str(params['disturbance_rate'] * params['uplift_interval']))


def initialize():
    """Initialize the model."""
    
    params = get_params()
    report_run_params(params)

def run_model():
    """Initialize, run, and finalize the hillslope model."""
    initialize()
    

def main():
    
    #for w_exp in range(-4, 5, 4):
    #    for d_exp in range(-4, 5, 4):
     #       run(1e7, w_exp, d_exp)
    #for d in range(0, 1):
    #run(1.0e9, 1.0e-8)
    run_model()


if __name__=='__main__':
    main()
