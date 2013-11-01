# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackend
import numpy

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=8),{'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : numpy.linspace(0.00065*2000,0.0008*2000,20),'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : numpy.linspace(0.0,0.004*2000,20)}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=16,num_mpi=1),{'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : numpy.linspace(0.0007,0.0011,5),'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : numpy.linspace(0.002,0.007,5)}).run_parameter_search()
CombinationParameterSearch(SlurmSequentialBackend(num_threads=8,num_mpi=1),{'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : numpy.linspace(0.0007,0.0011,5),'l4_cortex_exc.rand_struct_ratio' : numpy.linspace(0.5,1.0,10)}).run_parameter_search()


