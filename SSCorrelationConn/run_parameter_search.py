# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackend
import numpy

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=8),{'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : numpy.linspace(0.00065*2000,0.0008*2000,20),'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : numpy.linspace(0.0,0.004*2000,20)}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=16,num_mpi=1),{'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : numpy.linspace(0.0007,0.0011,5),'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : numpy.linspace(0.002,0.007,5)}).run_parameter_search()

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : numpy.linspace(0.0004,0.0005,30),'l4_cortex_exc.rand_struct_ratio' : numpy.linspace(0.5,1.0,20)}).run_parameter_search()

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : numpy.linspace(0.00035,0.00055,6),'l4_cortex_exc.AfferentConnection.num_samples' : numpy.linspace(20,90,4)}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : numpy.linspace(0.00040,0.00045,6),'l4_cortex_exc.AfferentConnection.num_samples' : numpy.linspace(1,40,3),'l4_cortex_inh.aff_inh_gain' : numpy.linspace(1.0,1.5,5)}).run_parameter_search()

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.AfferentConnection.num_samples' : numpy.linspace(10,100,10)}).run_parameter_search()

#'l4_cortex_exc.AfferentConnection.base_weight' : numpy.linspace(0.0,0.0004,2)

CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : numpy.linspace(0.0005,0.0007,5),'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : numpy.linspace(0.002,0.004,4)}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : numpy.linspace(0.0007,0.0011,5),'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : numpy.linspace(0.0005,0.002,5)}).run_parameter_search()


#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : numpy.linspace(0.00042,0.00044,3),'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : numpy.linspace(0.0048,0.0058,4),'l4_cortex_exc.AfferentConnection.num_samples' : numpy.linspace(20,90,2)}).run_parameter_search()

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_inh.aff_inh_gain' : numpy.linspace(1.0,2.0,5),'l4_cortex_exc.AfferentConnection.num_samples' : numpy.linspace(50,90,2),'l4_cortex_exc.rand_struct_ratio' : [0.85,0.9,0.99]}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_inh.aff_inh_gain' : numpy.linspace(1.0,2.0,5),'l4_cortex_exc.AfferentConnection.num_samples' : numpy.linspace(50,90,2)}).run_parameter_search()

