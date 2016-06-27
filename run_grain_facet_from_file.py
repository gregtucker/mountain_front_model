# -*- coding: utf-8 -*-
"""
Created on Sun Jun 26 09:13:46 2016

@author: gtucker
"""

import grain_facet_model as gfm

params = gfm.get_params_from_input_file('grain_facet_sample_inputs.txt')
gfm.main(params)
