# -*- coding: utf-8 -*-
"""


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

import grain_facet_model as gfm

params = {
    'number_of_node_rows' : 40,
    'number_of_node_columns' : 61,
    'report_interval' : 5.0,
    'run_duration' : 5.0e5,
    'plot_interval' : 5.0e4,
    'uplift_interval' : 5.0e3,
    'disturbance_rate' : 1.0e-4,
    'weathering_rate' : 1.0e-5,
    'friction_coef' : 1.0,
    'fault_x' : 8.0
    }
start = time.time()
gfm.main(params)
print 'Run time:', time.time() - start

