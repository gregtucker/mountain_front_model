#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
block_hill.py: version of GrainHill that adds blocks.

Created on Sat Jun 24 10:54:53 2017

@author: gtucker
"""

from grain_hill_as_class import GrainHill
from lattice_grain import lattice_grain_node_states


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
                 show_plots=True, **kwds):
        """Call the initialize() method."""
        self.initialize(grid_size, report_interval, run_duration,
                        output_interval, settling_rate, disturbance_rate,
                        weathering_rate, uplift_interval, plot_interval,
                        friction_coef, rock_state_for_uplift,
                        opt_rock_collapse, show_plots,
                        **kwds)

    def initialize(self, grid_size, report_interval, run_duration,
                   output_interval, settling_rate, disturbance_rate,
                   weathering_rate, uplift_interval, plot_interval,
                   friction_coef, rock_state_for_uplift, opt_rock_collapse,
                   show_plots, **kwds):
        """Initialize the BlockHill model."""
        # Call base class init
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
        
        # TODO: add transitions for blocks
        # Add setup of blocks in grid somewhere... maybe this goes in a version
        # of the lattice uplifter?

    def node_state_dictionary(self):
        """
        Create and return dict of node states.
        
        Overrides base-class method. Here, we call on a function in
        the lattice_grain module, and then add an additional state for blocks.
        """
        nsd = lattice_grain_node_states()
        nsd[9] = 'block'
        return nsd


if __name__ == '__main__':
    bh = BlockHill(grid_size=(3, 3))

        
