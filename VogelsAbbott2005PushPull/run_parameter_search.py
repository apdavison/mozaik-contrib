# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackend
import numpy
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=4,num_mpi=1),{'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : numpy.linspace(0.0,0.006,20),'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : numpy.linspace(0.0,0.05,20)}).run_parameter_search()
CombinationParameterSearch(SlurmSequentialBackend(num_threads=4,num_mpi=1),{'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : numpy.linspace(0.0,0.006,20),'l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.push_pull_ratio' : numpy.linspace(0.0,1.0,10)}).run_parameter_search()
