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

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : numpy.linspace(0.0007,0.0009,5),'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : numpy.linspace(0.0007,0.002,5)}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : numpy.linspace(0.0007,0.0011,5),'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : numpy.linspace(0.0005,0.002,5)}).run_parameter_search()

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : numpy.linspace(0.0007,0.0009,2),'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : numpy.linspace(0.0015,0.002,2), 'l4_cortex_inh.ExcInhAfferentRatio' : numpy.linspace(1.8,1.2,3),'retina_lgn.params.gain' : (4.0,6.0)}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : numpy.linspace(0.0005,0.00045,2),'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : numpy.linspace(0.0007,0.00065,2),'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : numpy.linspace(0.0006,0.0007,2), 'l4_cortex_inh.ExcInhAfferentRatio' : numpy.linspace(1.2,1.5,2),'retina_lgn.params.gain' : (4.0,6.0,8.0)}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.00045],'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0007],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0008], 'l4_cortex_inh.ExcInhAfferentRatio' : [1.0],'retina_lgn.params.gain' : [10.0,15.0,20.0],'l4_cortex_exc.rand_struct_ratio' : [0.9,0.85,0.95], 'l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [0.7,0.8,0.9]}).run_parameter_search()


#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.00045],'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0007],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0008], 'l4_cortex_inh.ExcInhAfferentRatio' : [1.0],'retina_lgn.params.gain' : [10.0,15.0,20.0],'l4_cortex_exc.rand_struct_ratio' : [0.9,0.85,0.95], 'l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [0.7,0.8,0.9],'l4_cortex_exc.params.cell.params.v_reset' : [-51,-55,-60]}).run_parameter_search()


#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : numpy.linspace(0.00042,0.0005,5),'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : numpy.linspace(0.001,0.005,5),'l4_cortex_exc.AfferentConnection.num_samples' : numpy.linspace(40,60,2)}).run_parameter_search()

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_inh.aff_inh_gain' : numpy.linspace(1.0,2.0,5),'l4_cortex_exc.AfferentConnection.num_samples' : numpy.linspace(50,90,2),'l4_cortex_exc.rand_struct_ratio' : [0.85,0.9,0.99]}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_inh.aff_inh_gain' : numpy.linspace(1.0,2.0,5),'l4_cortex_exc.AfferentConnection.num_samples' : numpy.linspace(50,90,2)}).run_parameter_search()

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=20),{'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.00045],'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0007],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0008], 'l4_cortex_inh.ExcInhAfferentRatio' : [1.0],'retina_lgn.params.gain' : [2.0,6.0],'l4_cortex_exc.rand_struct_ratio' : [0.85,0.95], 'l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [0.5,0.9],'l4_cortex_exc.params.cell.params.v_reset' : [-51], 'l4_cortex_exc.AfferentConnection.num_samples' : [20,50]}).run_parameter_search()


CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.00045,0.0006],'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0007,0.0003,0.0005,0.0009],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0007,0.0008,0.0009], 'l4_cortex_inh.ExcInhAfferentRatio' : [1.0],'retina_lgn.params.gain' : [3.0],'l4_cortex_exc.rand_struct_ratio' : [0.85,0.95], 'l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [0.5,0.9],'l4_cortex_exc.params.cell.params.v_reset' : [-51], 'l4_cortex_exc.AfferentConnection.num_samples' : [50]}).run_parameter_search()