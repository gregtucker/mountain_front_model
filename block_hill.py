#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
block_hill.py: version of GrainHill that adds blocks.

Created on Sat Jun 24 10:54:53 2017

@author: gtucker
"""

from grain_hill_as_class import GrainHill
from lattice_grain import (lattice_grain_node_states,
                           lattice_grain_transition_list)
from landlab.ca.boundaries.hex_lattice_tectonicizer import LatticeUplifter
from landlab.ca.celllab_cts import Transition


BLOCK_ID = 9


class BlockHill(GrainHill):
    """
    Model hillslope evolution with 'block' particles that can be undermined
    and weathered but not disturbed/activated.
    """
    def __init__(self, grid_size, report_interval=1.0e8, run_duration=1.0, 
                 output_interval=1.0e99, settling_rate=2.2e8,
                 disturbance_rate=1.0, weathering_rate=1.0, 
                 uplift_interval=1.0, plot_interval=1.0e99, friction_coef=0.3,
                 rock_state_for_uplift=7, opt_rock_collapse=False,
                 block_layer_dip_angle=0.0, block_layer_thickness=1.0,
                 show_plots=True, **kwds):
        """Call the initialize() method."""
        self.initialize(grid_size, report_interval, run_duration,
                        output_interval, settling_rate, disturbance_rate,
                        weathering_rate, uplift_interval, plot_interval,
                        friction_coef, rock_state_for_uplift,
                        opt_rock_collapse, block_layer_dip_angle,
                        block_layer_thickness, show_plots,
                        **kwds)

    def initialize(self, grid_size, report_interval, run_duration,
                   output_interval, settling_rate, disturbance_rate,
                   weathering_rate, uplift_interval, plot_interval,
                   friction_coef, rock_state_for_uplift, opt_rock_collapse,
                   block_layer_dip_angle, block_layer_thickness,
                   show_plots, **kwds):
        """Initialize the BlockHill model."""
        
        # Set block-related variables
        self.block_layer_dip_angle = block_layer_dip_angle
        self.block_layer_thickness = block_layer_thickness

        # Call parent class init
        super(BlockHill, self).initialize(grid_size=grid_size, 
                                          report_interval=report_interval, 
                                          run_duration=run_duration,
                                          output_interval=output_interval,
                                          settling_rate=settling_rate,
                                          disturbance_rate=disturbance_rate,
                                          weathering_rate=weathering_rate,
                                          uplift_interval=uplift_interval,
                                          plot_interval=plot_interval,
                                          friction_coef=friction_coef,
                                          rock_state_for_uplift=rock_state_for_uplift,
                                          opt_rock_collapse=opt_rock_collapse,
                                          show_plots=show_plots, **kwds)
        
        self.uplifter = LatticeUplifter(self.grid, 
                                self.grid.at_node['node_state'],
                                block_ID=BLOCK_ID,
                                block_layer_dip_angle=block_layer_dip_angle,
                                block_layer_thickness=block_layer_thickness)

    def node_state_dictionary(self):
        """
        Create and return dict of node states.
        
        Overrides base-class method. Here, we call on a function in
        the lattice_grain module, and then add an additional state for blocks.
        """
        nsd = lattice_grain_node_states()
        nsd[BLOCK_ID] = 'block'
        return nsd

    def add_block_transitions(self, xn_list):
        """Adds transitions for block undermining."""
        xn_list.append( Transition((0,BLOCK_ID,0), (BLOCK_ID,0,0),
                                   self.settling_rate, 'block settling') )
        xn_list.append( Transition((0,BLOCK_ID,1), (BLOCK_ID,0,1),
                                   self.settling_rate, 'block settling') )
        xn_list.append( Transition((0,BLOCK_ID,2), (BLOCK_ID,0,2),
                                   self.settling_rate, 'block settling') )

    def transition_list(self):
        """
        Make and return list of Transition object.
        """
        xn_list = lattice_grain_transition_list(g=self.settling_rate,
                                                f=self.friction_coef,
                                                motion=self.settling_rate)
        xn_list = super(BlockHill, self).add_weathering_and_disturbance_transitions(xn_list,
                    self.disturbance_rate, self.weathering_rate,
                    collapse_rate=self.collapse_rate)
        
        xn_list = self.add_block_transitions(xn_list)
        return xn_list
        

if __name__ == '__main__':
    bh = BlockHill(grid_size=(3, 3))

        
