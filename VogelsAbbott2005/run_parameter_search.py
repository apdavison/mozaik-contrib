# -*- coding: utf-8 -*-
import sys
sys.path.append('/home/jan/projects/mozaik/')
from mozaik.meta_workflow.parameter_search import run_parameter_search

run_parameter_search({'l4_cortex_exc.L4ExcL4ExcConnection.weights' : [0.0041,0.004],'l4_cortex_inh.L4InhL4ExcConnection.weights' : [0.52,0.51]})
