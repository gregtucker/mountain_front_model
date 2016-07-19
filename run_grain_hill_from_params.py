# -*- coding: utf-8 -*-
"""
Created on Sun Jun 26 09:13:46 2016

@author: gtucker
"""

from grain_hill_as_class import GrainHill

params = {
    'number_of_node_rows': 50,
    'number_of_node_columns': 115,
    'report_interval': 5.0,
    'run_duration': 1000.0,
    'output_interval': 1000.0,
    'disturbance_rate': 0.01,
    'uplift_interval': 100.0,
    'plot_interval': 100.0,
    'friction_coef': 1.0
    }

nr = params['number_of_node_rows']
nc = params['number_of_node_columns']
gh = GrainHill((nr, nc), **params)
