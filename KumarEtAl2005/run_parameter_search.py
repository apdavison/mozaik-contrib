# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackend
import numpy
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=4,num_mpi=1),{'l4_cortex_exc.L4ExcL4ExcConnection.weights' : [0.002,0.004],'l4_cortex_inh.L4InhL4ExcConnection.weights' : [0.2,0.51]}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=64),{'l4_cortex_exc.L4ExcL4ExcConnection.weights' : numpy.linspace(0.001,0.006,3),'l4_cortex_inh.L4InhL4ExcConnection.weights' : numpy.linspace(0.048,0.096,2)}).run_parameter_search()
CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=1),{'l4_cortex_exc.params.artificial_stimulators.background_act.params.exc_firing_rate' : numpy.linspace(0.0,100000,400)}).run_parameter_search()

