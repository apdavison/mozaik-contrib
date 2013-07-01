# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackend
import numpy
CombinationParameterSearch(SlurmSequentialBackend(num_threads=4,num_mpi=1),{'l4_cortex_exc.L4ExcL4ExcConnection.weights' : numpy.linspace(0.0,0.006,10),'l4_cortex_inh.L4InhL4ExcConnection.weights' : numpy.linspace(0.0,0.05,10)}).run_parameter_search()
