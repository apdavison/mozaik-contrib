# -*- coding: utf-8 -*-
import sys
import numpy
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch, SlurmSequentialBackend

CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=1),{'l4_cortex_exc.L4ExcL4ExcConnection.weights' : numpy.linspace(0.0,0.006,20),'l4_cortex_inh.L4InhL4ExcConnection.weights' : numpy.linspace(0.0,0.05,20)}).run_parameter_search()
