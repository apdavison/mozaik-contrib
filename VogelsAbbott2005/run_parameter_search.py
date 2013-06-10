# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackend
import numpy
CombinationParameterSearch(SlurmSequentialBackend(num_processes=1),{'l4_cortex_exc.L4ExcL4ExcConnection.weights' : numpy.linspace(0.0,0.01,20),'l4_cortex_inh.L4InhL4ExcConnection.weights' : numpy.linspace(0.0,0.1,20)}).run_parameter_search()
