# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackend
import numpy
CombinationParameterSearch(SlurmSequentialBackend(slurm_options=['-c','4']),{'l4_cortex_exc.L4ExcL4ExcConnection.weights' : numpy.linspace(0.0015,0.0025,20),'l4_cortex_inh.L4InhL4ExcConnection.weights' : numpy.linspace(0.015,0.025,20)}).run_parameter_search()
