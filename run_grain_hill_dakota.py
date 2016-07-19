#!/usr/env/python

import sys
from subprocess import call
import grain_hill_as_class
from landlab.io import load_params

input_file_template = 'my_inputs.template'
input_file = 'my_inputs.txt'

print 'hi from test_driver'
print sys.argv

# Run dprepro to get right value into the file my_inputs.txt
call(['dprepro', sys.argv[1], input_file_template, input_file])

# Read parameters from the input file
params = load_params(input_file)

grain_hill_as_class.main(params)



# Run the model

val = 0
myfile = open(sys.argv[2], 'w')
myfile.write( str(10*val) )


