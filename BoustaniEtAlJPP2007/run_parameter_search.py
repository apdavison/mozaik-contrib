# -*- coding: utf-8 -*-
import sys
import numpy
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch, SlurmSequentialBackend

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=4,num_mpi=1),{'l4_cortex_exc.L4ExcL4ExcConnection.weights' : numpy.linspace(0.0,0.005,20),'l4_cortex_inh.L4InhL4ExcConnection.weights' : numpy.linspace(0.0,0.05,20)}).run_parameter_search()
CombinationParameterSearch(SlurmSequentialBackend(num_threads=8,num_mpi=1),{'l4_cortex_exc.L4ExcL4ExcConnection.connection_probability' : numpy.linspace(0.0001,0.004,20),'l4_cortex_exc.params.density' : numpy.linspace(4000,204000,30)}).run_parameter_search()
