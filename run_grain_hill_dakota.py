#!/usr/bin/env python

import sys
from subprocess import call
from grain_hill_as_class import GrainHill
from landlab.core import load_params
import numpy as np


def two_node_diff(a):
    """Calculate and return diffs over two nodes instead of one."""
    N = len(a)
    return a[2:] - a[:(N-2)]

dx = 0.1  # assumed node spacing, m

input_file_template = '../my_inputs.template'
input_file = 'my_inputs.txt'

# Run dprepro to get right value into the file my_inputs.txt
call(['dprepro', sys.argv[1], input_file_template, input_file])

# Read parameters from the input file
params = load_params(input_file)

domain_length = 10.0 ** params['number_of_node_columns']
num_cols = int(np.round(domain_length / (dx * 0.866) + 1))
num_rows = int(np.round(0.5 * domain_length / dx))
params['number_of_node_columns'] = num_cols
params['number_of_node_rows'] = num_rows
params['disturbance_rate'] = 10.0 ** params['disturbance_rate']
params['uplift_interval'] = 10.0 ** params['uplift_interval']
params['run_duration'] = 0.5 * domain_length * params['uplift_interval'] / dx
params['plot_interval'] = 0.5 * params['run_duration']
params['output_interval'] = 0.25 * params['run_duration']


print params


# Run the model
gh = GrainHill((num_rows, num_cols), **params)
gh.run()


# Evaluate the results
(elev_profile, soil) = gh.get_profile_and_soil_thickness(gh.grid, 
                                                         gh.ca.node_state)
max_elev = np.amax(elev_profile)
N = len(elev_profile)
mean_grad_left = np.mean(two_node_diff(elev_profile[:((N+1)/2)])/1.73205)
mean_grad_right = np.mean(-two_node_diff(elev_profile[((N+1)/2):])/1.73205)
mean_grad = (mean_grad_left + mean_grad_right) / 2

myfile = open(sys.argv[2], 'w')
myfile.write(str(max_elev) + ' ' + str(mean_grad) + '\n')
myfile.close()
