# -*- coding: utf-8 -*-
"""
run_grain_facet_from_params.py: demonstrates how to instantiate and run
a GrainFacetSimulator from another Python script, passing parameters via
a dictionary rather than using a separate input file.

A note on time scales, in seconds:

Duration
(sec)  (equiv)
--------------
1 s    ~     1 s
10 s   ~     1 min
100 s  ~     1 min
1000 s ~     1 hr
10,000 s   ~  1 hr
10^5 s     ~  1 day (28 hrs)
10^6 s     ~  1 week (12 days)
10^7 s     ~  3 months
10^8 s     ~  3 years

Created on Sun Jun 26 09:13:46 2016

@author: gtucker
"""

import time
from grain_facet_model import GrainFacetSimulator

params = {
    'number_of_node_rows' : 20,
    'number_of_node_columns' : 31,
    'report_interval' : 5.0,
    'run_duration' : 150.0,
    'output_interval' : 1000.0,
    'plot_interval' : 10.0,
    'uplift_interval' : 10.0,
    'disturbance_rate' : 0.01,
    'weathering_rate' : 0.002,
    'friction_coef' : 1.0,
    'fault_x' : 8.0,
    'cell_width' : 1.0
    }

start = time.time()
gridsize = (params['number_of_node_rows'], params['number_of_node_columns'])
gfs = GrainFacetSimulator(gridsize, **params)
gfs.run()
print('Run complete. Run time (sec):')
print(time.time() - start)

