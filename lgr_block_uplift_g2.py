#!/usr/env/python
"""
Hillslope model with block uplift.
"""

_DEBUG = False

import time
#import random
from numpy import zeros, bincount, arange, savetxt, sqrt, log10, mean, arctan, pi
from pylab import subplots, plot, show, xlabel, ylabel, title, axis, figure, clf, savefig
from landlab import HexModelGrid
from landlab.ca.celllab_cts import Transition, CAPlotter
from landlab.ca.oriented_hex_cts import OrientedHexCTS
from landlab.ca.boundaries.hex_lattice_tectonicizer import LatticeUplifter
from scipy.optimize import curve_fit
import matplotlib

def setup_transition_list(g=0.0, f=0.0, d=0.0, w=0.0):
    """
    Creates and returns a list of Transition() objects to represent state
    transitions for simple granular mechanics model.
    
    Parameters
    ----------
    (none)
    
    Returns
    -------
    xn_list : list of Transition objects
        List of objects that encode information about the link-state transitions.

    Notes
    -----
    The transitions for this version of lattice gas have 11 pair-transition
    rules. The shorthand for the states is as follows:

        AR = air/empty
        IN = incoming particle (moving toward its neighbor)
        OU = outgoing particle (moving away from its neighbor)
        IO = incoming at an oblique angle
        OO = outgoing at an oblique angle
        RE = rest particle
        WA = wall particle
        op = oblique pair moving in opposite perpendicular direction
        sm = oblique pair moving in same perpendicular direction

    The 11 pairs with transitions are:

        1. AR-IN => IN-AR (move to empty particle)
        2. IN-IN => OO-OO-op (1/3 each dir), OU-OU (1/3) (head-on collision)
        3. IN-IO => OO-OU (oblique collision)
        4. IN-OO => IO-OU (oblique collision from behind)
        5. IN-OU => IO-OO (1/4 each of 2 directions) (collision from behind)
        6. IN-RE => RE-OU (1/3) RE-OO (1/3 each dir) (collision with rest)
        7. IN-WA => OU-WA (1/3) OO-WA (1/3 each dir) (wall collision)
        8. IO-IO-op => OO-OO-op (1/2 each dir) (glacing collision)
        9. IO-IO-sm => OO-OO-sm (30-degree collision)
        10. IO-RE => RE-OU (oblique collision with rest particle)
        11. IO-WA => OO-WA (oblique collision with wall)
    """
    xn_list = []
    
    p_elast = 1.0-f  # probability of elastic (non-dissipative) collision

    # Rule 1: Transitions for particle movement into an empty cell
    xn_list.append( Transition((1,0,0), (0,1,0), 1., 'motion') )
    xn_list.append( Transition((2,0,1), (0,2,1), 1., 'motion') )
    xn_list.append( Transition((3,0,2), (0,3,2), 1., 'motion') )
    xn_list.append( Transition((0,4,0), (4,0,0), 1., 'motion') )
    xn_list.append( Transition((0,5,1), (5,0,1), 1., 'motion') )
    xn_list.append( Transition((0,6,2), (6,0,2), 1., 'motion') )

    # Rule 2: Transitions for head-on collision: elastic
    xn_list.append( Transition((1,4,0), (4,1,0), p_elast/3, 'head-on collision') )
    xn_list.append( Transition((1,4,0), (3,6,0), p_elast/3, 'head-on collision') )
    xn_list.append( Transition((1,4,0), (5,2,0), p_elast/3, 'head-on collision') )
    xn_list.append( Transition((2,5,1), (5,2,1), p_elast/3, 'head-on collision') )
    xn_list.append( Transition((2,5,1), (4,1,1), p_elast/3, 'head-on collision') )
    xn_list.append( Transition((2,5,1), (6,3,1), p_elast/3, 'head-on collision') )
    xn_list.append( Transition((3,6,2), (6,3,2), p_elast/3, 'head-on collision') )
    xn_list.append( Transition((3,6,2), (1,4,2), p_elast/3, 'head-on collision') )
    xn_list.append( Transition((3,6,2), (5,2,2), p_elast/3, 'head-on collision') )

    # Rule 2: Transitions for head-on collision: frictional dissipation
    xn_list.append( Transition((1,4,0), (7,7,0), f, 'head-on collision') )
    xn_list.append( Transition((2,5,1), (7,7,1), f, 'head-on collision') )
    xn_list.append( Transition((3,6,2), (7,7,2), f, 'head-on collision') )

    # Rule 3: Transitions for oblique collision: elastic
    xn_list.append( Transition((1,3,0), (3,1,0), p_elast, 'oblique collision') )
    xn_list.append( Transition((1,5,0), (5,1,0), p_elast, 'oblique collision') )
    xn_list.append( Transition((2,4,0), (4,2,0), p_elast, 'oblique collision') )
    xn_list.append( Transition((6,4,0), (4,6,0), p_elast, 'oblique collision') )
    xn_list.append( Transition((2,4,1), (4,2,1), p_elast, 'oblique collision') )
    xn_list.append( Transition((2,6,1), (6,2,1), p_elast, 'oblique collision') )
    xn_list.append( Transition((1,5,1), (5,1,1), p_elast, 'oblique collision') )
    xn_list.append( Transition((3,5,1), (5,3,1), p_elast, 'oblique collision') )
    xn_list.append( Transition((3,1,2), (1,3,2), p_elast, 'oblique collision') )
    xn_list.append( Transition((3,5,2), (5,3,2), p_elast, 'oblique collision') )
    xn_list.append( Transition((2,6,2), (6,2,2), p_elast, 'oblique collision') )
    xn_list.append( Transition((4,6,2), (6,4,2), p_elast, 'oblique collision') )

    # Rule 3 frictional
    xn_list.append( Transition((1,3,0), (7,7,0), f, 'oblique collision') )
    xn_list.append( Transition((1,5,0), (7,7,0), f, 'oblique collision') )
    xn_list.append( Transition((2,4,0), (7,7,0), f, 'oblique collision') )
    xn_list.append( Transition((6,4,0), (7,7,0), f, 'oblique collision') )
    xn_list.append( Transition((2,4,1), (7,7,1), f, 'oblique collision') )
    xn_list.append( Transition((2,6,1), (7,7,1), f, 'oblique collision') )
    xn_list.append( Transition((1,5,1), (7,7,1), f, 'oblique collision') )
    xn_list.append( Transition((3,5,1), (7,7,1), f, 'oblique collision') )
    xn_list.append( Transition((3,1,2), (7,7,2), f, 'oblique collision') )
    xn_list.append( Transition((3,5,2), (7,7,2), f, 'oblique collision') )
    xn_list.append( Transition((2,6,2), (7,7,2), f, 'oblique collision') )
    xn_list.append( Transition((4,6,2), (7,7,2), f, 'oblique collision') )

    # Rule 4: Transitions for oblique-from-behind collisions
    xn_list.append( Transition((1,2,0), (2,1,0), p_elast, 'oblique') )
    xn_list.append( Transition((1,6,0), (6,1,0), p_elast, 'oblique') )
    xn_list.append( Transition((3,4,0), (4,3,0), p_elast, 'oblique') )
    xn_list.append( Transition((5,4,0), (4,5,0), p_elast, 'oblique') )
    xn_list.append( Transition((2,1,1), (1,2,1), p_elast, 'oblique') )
    xn_list.append( Transition((2,3,1), (3,2,1), p_elast, 'oblique') )
    xn_list.append( Transition((4,5,1), (5,4,1), p_elast, 'oblique') )
    xn_list.append( Transition((6,5,1), (5,6,1), p_elast, 'oblique') )
    xn_list.append( Transition((3,2,2), (2,3,2), p_elast, 'oblique') )
    xn_list.append( Transition((3,4,2), (4,3,2), p_elast, 'oblique') )
    xn_list.append( Transition((1,6,2), (6,1,2), p_elast, 'oblique') )
    xn_list.append( Transition((5,6,2), (6,5,2), p_elast, 'oblique') )

    # Rule 4 frictional
    xn_list.append( Transition((1,2,0), (7,1,0), f, 'oblique') )
    xn_list.append( Transition((1,6,0), (7,1,0), f, 'oblique') )
    xn_list.append( Transition((3,4,0), (4,7,0), f, 'oblique') )
    xn_list.append( Transition((5,4,0), (4,7,0), f, 'oblique') )
    xn_list.append( Transition((2,1,1), (7,2,1), f, 'oblique') )
    xn_list.append( Transition((2,3,1), (7,2,1), f, 'oblique') )
    xn_list.append( Transition((4,5,1), (5,7,1), f, 'oblique') )
    xn_list.append( Transition((6,5,1), (5,7,1), f, 'oblique') )
    xn_list.append( Transition((3,2,2), (7,3,2), f, 'oblique') )
    xn_list.append( Transition((3,4,2), (7,3,2), f, 'oblique') )
    xn_list.append( Transition((1,6,2), (6,7,2), f, 'oblique') )
    xn_list.append( Transition((5,6,2), (6,7,2), f, 'oblique') )
   
    # Rule 5: Transitions for direct-from-behind collisions
    xn_list.append( Transition((1,1,0), (2,6,0), p_elast/4, 'behind') )
    xn_list.append( Transition((1,1,0), (6,2,0), p_elast/4, 'behind') )
    xn_list.append( Transition((4,4,0), (3,5,0), p_elast/4, 'behind') )
    xn_list.append( Transition((4,4,0), (5,3,0), p_elast/4, 'behind') )
    xn_list.append( Transition((2,2,1), (1,3,1), p_elast/4, 'behind') )
    xn_list.append( Transition((2,2,1), (3,1,1), p_elast/4, 'behind') )
    xn_list.append( Transition((5,5,1), (4,6,1), p_elast/4, 'behind') )
    xn_list.append( Transition((5,5,1), (6,4,1), p_elast/4, 'behind') )
    xn_list.append( Transition((3,3,2), (2,4,2), p_elast/4, 'behind') )
    xn_list.append( Transition((3,3,2), (4,2,2), p_elast/4, 'behind') )
    xn_list.append( Transition((6,6,2), (1,5,2), p_elast/4, 'behind') )
    xn_list.append( Transition((6,6,2), (5,1,2), p_elast/4, 'behind') )

    # Rule 5 frictional
    xn_list.append( Transition((1,1,0), (7,1,0), f/4, 'behind') )
    xn_list.append( Transition((4,4,0), (4,7,0), f/4, 'behind') )
    xn_list.append( Transition((2,2,1), (7,2,1), f/4, 'behind') )
    xn_list.append( Transition((5,5,1), (5,7,1), f/4, 'behind') )
    xn_list.append( Transition((3,3,2), (7,3,2), f/4, 'behind') )
    xn_list.append( Transition((6,6,2), (6,7,2), f/4, 'behind') )

    # Rule 6: Transitions for direct collision with stationary (resting) particle
    xn_list.append( Transition((1,7,0), (7,1,0), p_elast/3., 'rest') )
    xn_list.append( Transition((1,7,0), (7,2,0), p_elast/3., 'rest') )
    xn_list.append( Transition((1,7,0), (7,6,0), p_elast/3., 'rest') )
    xn_list.append( Transition((7,4,0), (4,7,0), p_elast/3., 'rest') )
    xn_list.append( Transition((7,4,0), (3,7,0), p_elast/3., 'rest') )
    xn_list.append( Transition((7,4,0), (5,7,0), p_elast/3., 'rest') )
    xn_list.append( Transition((2,7,1), (7,2,1), p_elast/3., 'rest') )
    xn_list.append( Transition((2,7,1), (7,1,1), p_elast/3., 'rest') )
    xn_list.append( Transition((2,7,1), (7,3,1), p_elast/3., 'rest') )
    xn_list.append( Transition((7,5,1), (5,7,1), p_elast/3., 'rest') )
    xn_list.append( Transition((7,5,1), (4,7,1), p_elast/3., 'rest') )
    xn_list.append( Transition((7,5,1), (6,7,1), p_elast/3., 'rest') )
    xn_list.append( Transition((3,7,2), (7,3,2), p_elast/3., 'rest') )
    xn_list.append( Transition((3,7,2), (7,2,2), p_elast/3., 'rest') )
    xn_list.append( Transition((3,7,2), (7,4,2), p_elast/3., 'rest') )
    xn_list.append( Transition((7,6,2), (6,7,2), p_elast/3., 'rest') )
    xn_list.append( Transition((7,6,2), (1,7,2), p_elast/3., 'rest') )
    xn_list.append( Transition((7,6,2), (5,7,2), p_elast/3., 'rest') )

    # Rule 6 frictionl
    xn_list.append( Transition((1,7,0), (7,7,0), f, 'rest') )
    xn_list.append( Transition((7,4,0), (7,7,0), f, 'rest') )
    xn_list.append( Transition((2,7,1), (7,7,1), f, 'rest') )
    xn_list.append( Transition((7,5,1), (7,7,1), f, 'rest') )
    xn_list.append( Transition((3,7,2), (7,7,2), f, 'rest') )
    xn_list.append( Transition((7,6,2), (7,7,2), f, 'rest') )

    # Rule 7: Transitions for wall impact
    xn_list.append( Transition((1,8,0), (4,8,0), p_elast/3, 'wall rebound') )
    xn_list.append( Transition((1,8,0), (3,8,0), p_elast/3, 'wall rebound') )
    xn_list.append( Transition((1,8,0), (5,8,0), p_elast/3, 'wall rebound') )
    xn_list.append( Transition((2,8,1), (5,8,1), p_elast/3, 'wall rebound') )
    xn_list.append( Transition((2,8,1), (4,8,1), p_elast/3, 'wall rebound') )
    xn_list.append( Transition((2,8,1), (6,8,1), p_elast/3, 'wall rebound') )
    xn_list.append( Transition((3,8,2), (6,8,2), p_elast/3, 'wall rebound') )
    xn_list.append( Transition((3,8,2), (5,8,2), p_elast/3, 'wall rebound') )
    xn_list.append( Transition((3,8,2), (1,8,2), p_elast/3, 'wall rebound') )
    xn_list.append( Transition((8,4,0), (8,1,0), p_elast/3, 'wall rebound') )
    xn_list.append( Transition((8,4,0), (8,6,0), p_elast/3, 'wall rebound') )
    xn_list.append( Transition((8,4,0), (8,2,0), p_elast/3, 'wall rebound') )
    xn_list.append( Transition((8,5,1), (8,1,1), p_elast/3, 'wall rebound') )
    xn_list.append( Transition((8,5,1), (8,2,1), p_elast/3, 'wall rebound') )
    xn_list.append( Transition((8,5,1), (8,3,1), p_elast/3, 'wall rebound') )
    xn_list.append( Transition((8,6,2), (8,2,2), p_elast/3, 'wall rebound') )
    xn_list.append( Transition((8,6,2), (8,3,2), p_elast/3, 'wall rebound') )
    xn_list.append( Transition((8,6,2), (8,4,2), p_elast/3, 'wall rebound') )

    # Rule 7 frictional
    xn_list.append( Transition((1,8,0), (7,8,0), f, 'wall rebound') )
    xn_list.append( Transition((2,8,1), (7,8,1), f, 'wall rebound') )
    xn_list.append( Transition((3,8,2), (7,8,2), f, 'wall rebound') )
    xn_list.append( Transition((8,4,0), (8,7,0), f, 'wall rebound') )
    xn_list.append( Transition((8,5,1), (8,7,1), f, 'wall rebound') )
    xn_list.append( Transition((8,6,2), (8,7,2), f, 'wall rebound') )

    # Rule 8: Transitions for glancing oblique collision
    xn_list.append( Transition((2,5,0), (3,6,0), p_elast, 'glancing') )
    xn_list.append( Transition((6,3,0), (5,2,0), p_elast, 'glancing') )
    xn_list.append( Transition((3,6,1), (4,1,1), p_elast, 'glancing') )
    xn_list.append( Transition((1,4,1), (6,3,1), p_elast, 'glancing') )
    xn_list.append( Transition((4,1,2), (5,2,2), p_elast, 'glancing') )
    xn_list.append( Transition((2,5,2), (1,4,2), p_elast, 'glancing') )
    
    # Rule 8 frictional
    xn_list.append( Transition((2,5,0), (7,7,0), f, 'glancing') )
    xn_list.append( Transition((6,3,0), (7,7,0), f, 'glancing') )
    xn_list.append( Transition((3,6,1), (7,7,1), f, 'glancing') )
    xn_list.append( Transition((1,4,1), (7,7,1), f, 'glancing') )
    xn_list.append( Transition((4,1,2), (7,7,2), f, 'glancing') )
    xn_list.append( Transition((2,5,2), (7,7,2), f, 'glancing') )

    # Rule 9: Transitions for "near-on" collisions
    xn_list.append( Transition((6,5,0), (5,6,0), p_elast, 'near-on') )
    xn_list.append( Transition((2,3,0), (3,2,0), p_elast, 'near-on') )
    xn_list.append( Transition((1,6,1), (6,1,1), p_elast, 'near-on') )
    xn_list.append( Transition((3,4,1), (4,3,1), p_elast, 'near-on') )
    xn_list.append( Transition((2,1,2), (1,2,2), p_elast, 'near-on') )
    xn_list.append( Transition((4,5,2), (5,4,2), p_elast, 'near-on') )
    
    # Rule 9 frictional
    xn_list.append( Transition((6,5,0), (7,6,0), f/2, 'near-on') )
    xn_list.append( Transition((6,5,0), (5,7,0), f/2, 'near-on') )
    xn_list.append( Transition((2,3,0), (7,2,0), f/2, 'near-on') )
    xn_list.append( Transition((2,3,0), (3,7,0), f/2, 'near-on') )
    xn_list.append( Transition((1,6,1), (7,1,1), f/2, 'near-on') )
    xn_list.append( Transition((1,6,1), (6,7,1), f/2, 'near-on') )
    xn_list.append( Transition((3,4,1), (7,3,1), f/2, 'near-on') )
    xn_list.append( Transition((3,4,1), (4,7,1), f/2, 'near-on') )
    xn_list.append( Transition((2,1,2), (7,2,2), f/2, 'near-on') )
    xn_list.append( Transition((2,1,2), (1,7,2), f/2, 'near-on') )
    xn_list.append( Transition((4,5,2), (7,4,2), f/2, 'near-on') )
    xn_list.append( Transition((4,5,2), (5,7,2), f/2, 'near-on') )
    
    # Rule 10: Transitions for oblique collision with rest particle
    xn_list.append( Transition((2,7,0), (7,1,0), p_elast, 'oblique with rest') )
    xn_list.append( Transition((6,7,0), (7,1,0), p_elast, 'oblique with rest') )
    xn_list.append( Transition((7,3,0), (4,7,0), p_elast, 'oblique with rest') )
    xn_list.append( Transition((7,5,0), (4,7,0), p_elast, 'oblique with rest') )
    xn_list.append( Transition((3,7,1), (7,2,1), p_elast, 'oblique with rest') )
    xn_list.append( Transition((1,7,1), (7,2,1), p_elast, 'oblique with rest') )
    xn_list.append( Transition((7,6,1), (5,7,1), p_elast, 'oblique with rest') )
    xn_list.append( Transition((7,4,1), (5,7,1), p_elast, 'oblique with rest') )
    xn_list.append( Transition((4,7,2), (7,3,2), p_elast, 'oblique with rest') )
    xn_list.append( Transition((2,7,2), (7,3,2), p_elast, 'oblique with rest') )
    xn_list.append( Transition((7,5,2), (6,7,2), p_elast, 'oblique with rest') )
    xn_list.append( Transition((7,1,2), (6,7,2), p_elast, 'oblique with rest') )

    # Rule 10 frictional
    xn_list.append( Transition((2,7,0), (7,7,0), f, 'oblique with rest') )
    xn_list.append( Transition((6,7,0), (7,7,0), f, 'oblique with rest') )
    xn_list.append( Transition((7,3,0), (7,7,0), f, 'oblique with rest') )
    xn_list.append( Transition((7,5,0), (7,7,0), f, 'oblique with rest') )
    xn_list.append( Transition((3,7,1), (7,7,1), f, 'oblique with rest') )
    xn_list.append( Transition((1,7,1), (7,7,1), f, 'oblique with rest') )
    xn_list.append( Transition((7,6,1), (7,7,1), f, 'oblique with rest') )
    xn_list.append( Transition((7,4,1), (7,7,1), f, 'oblique with rest') )
    xn_list.append( Transition((4,7,2), (7,7,2), f, 'oblique with rest') )
    xn_list.append( Transition((2,7,2), (7,7,2), f, 'oblique with rest') )
    xn_list.append( Transition((7,5,2), (7,7,2), f, 'oblique with rest') )
    xn_list.append( Transition((7,1,2), (7,7,2), f, 'oblique with rest') )

    # Rule 11: Transitions for oblique collision with wall particle
    xn_list.append( Transition((2,8,0), (3,8,0), p_elast, 'oblique with wall') )
    xn_list.append( Transition((6,8,0), (5,8,0), p_elast, 'oblique with wall') )
    xn_list.append( Transition((8,3,0), (8,2,0), p_elast, 'oblique with wall') )
    xn_list.append( Transition((8,5,0), (8,6,0), p_elast, 'oblique with wall') )
    xn_list.append( Transition((1,8,1), (6,8,1), p_elast, 'oblique with wall') )
    xn_list.append( Transition((3,8,1), (4,8,1), p_elast, 'oblique with wall') )
    xn_list.append( Transition((8,4,1), (8,3,1), p_elast, 'oblique with wall') )
    xn_list.append( Transition((8,6,1), (8,1,1), p_elast, 'oblique with wall') )
    xn_list.append( Transition((4,8,2), (5,8,2), p_elast, 'oblique with wall') )
    xn_list.append( Transition((2,8,2), (1,8,2), p_elast, 'oblique with wall') )
    xn_list.append( Transition((8,1,2), (8,2,2), p_elast, 'oblique with wall') )
    xn_list.append( Transition((8,5,2), (8,4,2), p_elast, 'oblique with wall') )

    # Rule 11 frictional
    xn_list.append( Transition((2,8,0), (7,8,0), f, 'oblique with wall') )
    xn_list.append( Transition((6,8,0), (7,8,0), f, 'oblique with wall') )
    xn_list.append( Transition((8,3,0), (8,7,0), f, 'oblique with wall') )
    xn_list.append( Transition((8,5,0), (8,7,0), f, 'oblique with wall') )
    xn_list.append( Transition((1,8,1), (7,8,1), f, 'oblique with wall') )
    xn_list.append( Transition((3,8,1), (7,8,1), f, 'oblique with wall') )
    xn_list.append( Transition((8,4,1), (8,7,1), f, 'oblique with wall') )
    xn_list.append( Transition((8,6,1), (8,7,1), f, 'oblique with wall') )
    xn_list.append( Transition((4,8,2), (7,8,2), f, 'oblique with wall') )
    xn_list.append( Transition((2,8,2), (7,8,2), f, 'oblique with wall') )
    xn_list.append( Transition((8,1,2), (8,7,2), f, 'oblique with wall') )
    xn_list.append( Transition((8,5,2), (8,7,2), f, 'oblique with wall') )

    # Gravity rule 1: rising particles become rest particles
    xn_list.append( Transition((0,1,0), (0,7,0), g, 'gravity 1') )
    xn_list.append( Transition((1,1,0), (1,7,0), g, 'gravity 1') )
    xn_list.append( Transition((2,1,0), (2,7,0), g, 'gravity 1') )
    xn_list.append( Transition((3,1,0), (3,7,0), g, 'gravity 1') )
    xn_list.append( Transition((4,1,0), (4,7,0), g, 'gravity 1') )
    xn_list.append( Transition((5,1,0), (5,7,0), g, 'gravity 1') )
    xn_list.append( Transition((6,1,0), (6,7,0), g, 'gravity 1') )
    xn_list.append( Transition((7,1,0), (7,7,0), g, 'gravity 1') )
    xn_list.append( Transition((8,1,0), (8,7,0), g, 'gravity 1') )

    # Gravity rule 2: resting particles become falling particles (if not above
    # rest or wall?)
    xn_list.append( Transition((0,7,0), (0,4,0), g, 'gravity 2') )
    xn_list.append( Transition((1,7,0), (1,4,0), g, 'gravity 2') )
    xn_list.append( Transition((2,7,0), (2,4,0), g, 'gravity 2') )
    xn_list.append( Transition((3,7,0), (3,4,0), g, 'gravity 2') )
    xn_list.append( Transition((4,7,0), (4,4,0), g, 'gravity 2') )
    xn_list.append( Transition((5,7,0), (5,4,0), g, 'gravity 2') )
    xn_list.append( Transition((6,7,0), (6,4,0), g, 'gravity 2') )

    # Gravity rule 3: up/sideways particles become down/sideways particles
    xn_list.append( Transition((0,2,0), (0,3,0), g, 'gravity 3') )
    xn_list.append( Transition((1,2,0), (1,3,0), g, 'gravity 3') )
    xn_list.append( Transition((2,2,0), (2,3,0), g, 'gravity 3') )
    xn_list.append( Transition((3,2,0), (3,3,0), g, 'gravity 3') )
    xn_list.append( Transition((4,2,0), (4,3,0), g, 'gravity 3') )
    xn_list.append( Transition((5,2,0), (5,3,0), g, 'gravity 3') )
    xn_list.append( Transition((6,2,0), (6,3,0), g, 'gravity 3') )
    xn_list.append( Transition((7,2,0), (7,3,0), g, 'gravity 3') )
    xn_list.append( Transition((8,2,0), (8,3,0), g, 'gravity 3') )
    xn_list.append( Transition((0,6,0), (0,5,0), g, 'gravity 3') )
    xn_list.append( Transition((1,6,0), (1,5,0), g, 'gravity 3') )
    xn_list.append( Transition((2,6,0), (2,5,0), g, 'gravity 3') )
    xn_list.append( Transition((3,6,0), (3,5,0), g, 'gravity 3') )
    xn_list.append( Transition((4,6,0), (4,5,0), g, 'gravity 3') )
    xn_list.append( Transition((5,6,0), (5,5,0), g, 'gravity 3') )
    xn_list.append( Transition((6,6,0), (6,5,0), g, 'gravity 3') )
    xn_list.append( Transition((7,6,0), (7,5,0), g, 'gravity 3') )
    xn_list.append( Transition((8,6,0), (8,5,0), g, 'gravity 3') )
    
    # Gravity rule 4: down/side to straight down
    xn_list.append( Transition((0,3,0), (0,4,0), g, 'gravity 4') )
    xn_list.append( Transition((1,3,0), (1,4,0), g, 'gravity 4') )
    xn_list.append( Transition((2,3,0), (2,4,0), g, 'gravity 4') )
    xn_list.append( Transition((3,3,0), (3,4,0), g, 'gravity 4') )
    xn_list.append( Transition((4,3,0), (4,4,0), g, 'gravity 4') )
    xn_list.append( Transition((5,3,0), (5,4,0), g, 'gravity 4') )
    xn_list.append( Transition((6,3,0), (6,4,0), g, 'gravity 4') )
    xn_list.append( Transition((7,3,0), (7,4,0), g, 'gravity 4') )
    xn_list.append( Transition((8,3,0), (8,4,0), g, 'gravity 4') )
    xn_list.append( Transition((0,5,0), (0,4,0), g, 'gravity 4') )
    xn_list.append( Transition((1,5,0), (1,4,0), g, 'gravity 4') )
    xn_list.append( Transition((2,5,0), (2,4,0), g, 'gravity 4') )
    xn_list.append( Transition((3,5,0), (3,4,0), g, 'gravity 4') )
    xn_list.append( Transition((4,5,0), (4,4,0), g, 'gravity 4') )
    xn_list.append( Transition((5,5,0), (5,4,0), g, 'gravity 4') )
    xn_list.append( Transition((6,5,0), (6,4,0), g, 'gravity 4') )
    xn_list.append( Transition((7,5,0), (7,4,0), g, 'gravity 4') )
    xn_list.append( Transition((8,5,0), (8,4,0), g, 'gravity 4') )
    
    # Uncertain gravity rule!
    xn_list.append( Transition((7,0,2), (3,0,2), g/2.0, 'gravity??') )
    xn_list.append( Transition((0,7,1), (0,5,1), g/2.0, 'gravity??') )
    
    
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
    
    
def run(uplift_interval, d_ratio_exp):

    # INITIALIZE

    # User-defined parameters
    nr = 81
    nc = 151
    g = 1.0
    f = 0.7
    #d_ratio_exp = -8.0
    #w_ratio_exp = -8.0
    plot_interval = 0.1
    run_duration = 10.0
    report_interval = 5.0  # report interval, in real-time seconds
    plot_every_transition = False
    #uplift_interval = 1e7
    #filenm = 'Hill141213/h141213'
    #imagenm = 'Hill141213/hill'+str(int(d_ratio_exp))+'d'
    
    # Calculate d and w
    upliftrate = sqrt(3)/uplift_interval
    #w = (2.0**w_ratio_exp)*sliprate
    d = (2.0**d_ratio_exp)*upliftrate
    print 'RUNNING u=', upliftrate, 'd=', d
    
    # Assemble the file name
    #filenm = filenm+'u'+str(int(log10(uplift_interval)))+'d'+str(int(d_ratio_exp))
    #print filenm
    
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

    # Create data and initialize values.
    node_state_grid = hmg.add_zeros('node', 'node_state_grid', dtype=int)

    # Lower rows get resting particles
    for i in range(hmg.number_of_nodes):
        if hmg.node_y[i]< 2.0 and hmg.node_x[i]>0.0:
            node_state_grid[i] = 7
        #elif hmg.node_x[i]>((nc-1)*0.866):
        #    node_state_grid[i] = 8
    #node_state_grid[0] = 8
    node_state_grid[hmg.number_of_node_rows-1] = 8

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
    rock = sky
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

    # Create an array to store the numbers of states at each plot interval
    #nstates = zeros((9, int(run_duration/plot_interval)))
    #k = 0

    # RUN
    current_time = 0.0
    while current_time < run_duration:
        
        # Once in a while, print out simulation and real time to let the user
        # know that the sim is running ok
        current_real_time = time.time()
        if current_real_time >= next_report:
            print 'Current sim time',current_time,'(',100*current_time/run_duration,'%)'
            next_report = current_real_time + report_interval

        # Run the model forward in time until the next output step
        ca.run(current_time+plot_interval, ca.node_state, 
               plot_each_transition=plot_every_transition, plotter=ca_plotter)
        current_time += plot_interval

        # Plot the current grid
        ca_plotter.update_plot()
        axis('off')
        node_state_grid[hmg.number_of_node_rows-1] = 8
        #savefig(imagenm+str(k)+'.png')
        k+=1

        if current_time >= next_uplift:
            uplifter.uplift_interior_nodes(rock_state=7)
            ca.update_link_states_and_transitions(current_time)
            next_uplift += uplift_interval


    # FINALIZE

    # Plot
    ca_plotter.finalize()
    axis('off')
    
    (elev, soil) = get_profile_and_soil_thickness(hmg, node_state_grid)
    #print elev
    #print soil
    
    print 'Mean thickness:',mean(soil[75:])
    
    x = 0.5*sqrt(3)*arange(nc)
    
    # Calculate slope in upper half
    (cf,cm1) = curve_fit(flin, x[75:], elev[75:])
    slpgrad = cf[0]
    slpdeg = 180.*arctan(slpgrad)/pi
    print 'Upper slope gradient', slpgrad,'angle',slpdeg
    
    # Calculate (roughly) Laplacian
    (pf,cm2) = curve_fit(fquad, x[75:], elev[75:] )
    second_deriv = 2*pf[0]
    print '2nd deriv',second_deriv
    
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


def main():
    
    #for w_exp in range(-4, 5, 4):
    #    for d_exp in range(-4, 5, 4):
     #       run(1e7, w_exp, d_exp)
    for d in range(0, 1):
        run(1e7, d)


if __name__=='__main__':
    main()
