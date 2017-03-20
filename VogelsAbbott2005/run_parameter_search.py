# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackend
import numpy
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=4,num_mpi=1),{'l4_cortex_exc.L4ExcL4ExcConnection.weights' : [0.002,0.004],'l4_cortex_inh.L4InhL4ExcConnection.weights' : [0.2,0.51]}).run_parameter_search()
CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=8),{'l4_cortex_exc.L4ExcL4ExcConnection.weights' : numpy.linspace(0.002,0.006,5),'l4_cortex_inh.L4InhL4ExcConnection.weights' : numpy.linspace(0.001,0.03,5),
'l4_cortex_exc.params.cell.params.tau_m' : [10,11,12]}).run_parameter_search()
