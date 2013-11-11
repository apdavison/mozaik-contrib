# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackend
import numpy
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=4,num_mpi=1),{'l4_cortex_exc.L4ExcL4ExcConnection.weights' : [0.002,0.004],'l4_cortex_inh.L4InhL4ExcConnection.weights' : [0.2,0.51]}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4ExcConnection.weights' : numpy.linspace(0.0,0.006,40),'l4_cortex_inh.L4InhL4ExcConnection.weights' : numpy.linspace(0.0,0.1,40)}).run_parameter_search()

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=8),{'l4_cortex_exc.L4ExcL4InhConnection.weights' : numpy.linspace(0.0006,0.0011,40),'l4_cortex_inh.L4InhL4ExcConnection.weights' : numpy.linspace(0.0,0.006,40)}).run_parameter_search()


#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=4),{'l4_cortex_exc.params.density' : numpy.linspace(5000,50000,10),'l4_cortex_inh.L4InhL4ExcConnection.weights' : numpy.linspace(0.0,0.1,20)}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=4),{'l4_cortex_exc.params.cell.params.tau_m' : numpy.linspace(5,50,10),'l4_cortex_inh.L4InhL4ExcConnection.weights' : numpy.linspace(0.0,0.1,20)}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=4),{'l4_cortex_exc.params.cell.params.tau_syn_E' : numpy.linspace(0.1,1.0,20),'l4_cortex_inh.L4InhL4ExcConnection.weights' : numpy.linspace(0.0,0.1,20)}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=4),{'l4_cortex_exc.L4ExcL4InhConnection.weights' : numpy.linspace(0.003,0.01,20),'l4_cortex_inh.L4InhL4ExcConnection.weights' : numpy.linspace(0.0,0.1,20)}).run_parameter_search()


CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=8),{'l4_cortex_exc.L4ExcL4InhConnection.weights' : numpy.linspace(0.0006,0.0008,20)}).run_parameter_search()