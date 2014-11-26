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

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : numpy.linspace(0.0007,0.0009,2),'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : numpy.linspace(0.0015,0.002,2), 'l4_cortex_inh.ExcInhAfferentRatio' : numpy.linspace(1.8,1.2,3),: (4.0,6.0)}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : numpy.linspace(0.0005,0.00045,2),'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : numpy.linspace(0.0007,0.00065,2),'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : numpy.linspace(0.0006,0.0007,2), 'l4_cortex_inh.ExcInhAfferentRatio' : numpy.linspace(1.2,1.5,2),'retina_lgn.params.gain' : (4.0,6.0,8.0)}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.00045],'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0007],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0008], 'l4_cortex_inh.ExcInhAfferentRatio' : [1.0],'retina_lgn.params.gain' : [10.0,15.0,20.0],'l4_cortex_exc.rand_struct_ratio' : [0.9,0.85,0.95], 'l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [0.7,0.8,0.9]}).run_parameter_search()


#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.00045],'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0007],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0008], 'l4_cortex_inh.ExcInhAfferentRatio' : [1.0],'retina_lgn.params.gain' : [10.0,15.0,20.0],'l4_cortex_exc.rand_struct_ratio' : [0.9,0.85,0.95], 'l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [0.7,0.8,0.9],'l4_cortex_exc.params.cell.params.v_reset' : [-51,-55,-60]}).run_parameter_search()


#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : numpy.linspace(0.00042,0.0005,5),'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : numpy.linspace(0.001,0.005,5),'l4_cortex_exc.AfferentConnection.num_samples' : numpy.linspace(40,60,2)}).run_parameter_search()

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_inh.aff_inh_gain' : numpy.linspace(1.0,2.0,5),'l4_cortex_exc.AfferentConnection.num_samples' : numpy.linspace(50,90,2),'l4_cortex_exc.rand_struct_ratio' : [0.85,0.9,0.99]}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_inh.aff_inh_gain' : numpy.linspace(1.0,2.0,5),'l4_cortex_exc.AfferentConnection.num_samples' : numpy.linspace(50,90,2)}).run_parameter_search()

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=20),{'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.00045],'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0007],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0008], 'l4_cortex_inh.ExcInhAfferentRatio' : [1.0],'retina_lgn.params.gain' : [2.0,6.0],'l4_cortex_exc.rand_struct_ratio' : [0.85,0.95], 'l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [0.5,0.9],'l4_cortex_exc.params.cell.params.v_reset' : [-51], 'l4_cortex_exc.AfferentConnection.num_samples' : [20,50]}).run_parameter_search()


#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.00045],'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0011,0.0013,0.00175],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0008,0.001,0.0012], 'l4_cortex_inh.ExcInhAfferentRatio' : [1.0,1.1,1.3],'retina_lgn.params.gain' : [2.0],'l4_cortex_exc.rand_struct_ratio' : [0.9], 'l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [0.9],'l4_cortex_exc.params.cell.params.v_reset' : [-51], 'l4_cortex_exc.AfferentConnection.num_samples' : [50]}).run_parameter_search()


#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.0004],'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.00005],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0002,0.0006], 'l4_cortex_inh.ExcInhAfferentRatio' : [1.5,2.0],'retina_lgn.params.gain' : [1.0,2.0],'l4_cortex_exc.rand_struct_ratio' : [0.85], 'l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [0.3,0.9],'l4_cortex_exc.params.cell.params.v_rest' : [-80],'l4_cortex_exc.params.cell.params.v_reset' : [-55],'l4_cortex_exc.AfferentConnection.base_weight' : [0.0004,0.0006]}).run_parameter_search()


#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.0005],'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0004,0.0005,0.0006],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0006,0.0005,0.0004,0.0003,0.0002], 'l4_cortex_inh.ExcInhAfferentRatio' : [2.0],'retina_lgn.params.gain' : [0.5],'l4_cortex_exc.rand_struct_ratio' : [0.85], 'l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [0.5],'l4_cortex_exc.params.cell.params.v_rest' : [-70],'l4_cortex_exc.params.cell.params.v_reset' : [-51],'l4_cortex_exc.AfferentConnection.base_weight' : [0.0012],'l4_cortex_exc.AfferentConnection.num_samples' : [25],'l4_cortex_exc.params.cell.params.cm' : [0.1],'l4_cortex_exc.params.cell.params.v_thresh' : [-50]}).run_parameter_search()

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=20),{'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0003,0.0004],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0005,0.0006],'l4_cortex_exc.rand_struct_ratio' : [0.8,0.7],'retina_lgn.params.gain' : [0.5,0.8],'l4_cortex_inh.L4InhL4InhConnection.base_weight' : [0.0006,0.0008]}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=20),{'retina_lgn.params.gain' : [0.3,0.5,0.7],'retina_lgn.params.noise.stdev' : [2.7,3.3,4.0], 'l4_cortex_exc.AfferentConnection.num_samples' : [50,40,30]}).run_parameter_search()

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.rand_struct_ratio' : [0.6], 'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0005,0.0004,0.0003],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0002,0.0003,0.0005], 'l4_cortex_inh.params.cell.params.v_thresh' : [-50,-53],'retina_lgn.params.gain' : [0.7],'retina_lgn.params.noise.stdev' : [2.7,3.3,4.0]}).run_parameter_search()

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.rand_struct_ratio' : [0.8], 'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0003,0.0001],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0005,0.0006], 'l4_cortex_inh.params.cell.params.v_thresh' : [-53],'retina_lgn.params.gain' : [0.7],'l4_cortex_inh.ExcInhAfferentRatio' : [1.7,2.0],'l4_cortex_exc.AfferentConnection.base_weight' : [0.0007],'retina_lgn.params.noise.stdev' : [2.7,3.5]}).run_parameter_search()


#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.rand_struct_ratio' : [0.6], 'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.0006,0.0007],'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0003],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0005,0.001], 'l4_cortex_inh.params.cell.params.v_thresh' : [-53],'retina_lgn.params.gain' : [0.9,1.3],'l4_cortex_inh.ExcInhAfferentRatio' : [1.3],'l4_cortex_exc.AfferentConnection.base_weight' : [0.0007],'retina_lgn.params.noise.stdev' : [3.5],'l4_cortex_inh.L4InhL4InhConnection.base_weight' : [0.0002,0.0003,0.0005]}).run_parameter_search()


#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.rand_struct_ratio' : [0.6,0.7], 'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.0006,0.00063,0.00065,0.00067],'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0003],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0005], 'l4_cortex_inh.params.cell.params.v_thresh' : [-53],'retina_lgn.params.gain' : [1.3],'l4_cortex_inh.ExcInhAfferentRatio' : [1.3],'l4_cortex_exc.AfferentConnection.base_weight' : [0.0025],'l4_cortex_exc.AfferentConnection.num_samples' : [10],'retina_lgn.params.noise.stdev' : [3.5,5.0],'l4_cortex_inh.L4InhL4InhConnection.base_weight' : [0.0002,0.00025]}).run_parameter_search()


#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.rand_struct_ratio' : [0.6], 'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.0],'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0005,0.0004,0.0003,0.0006], 'l4_cortex_inh.params.cell.params.v_thresh' : [-53],'retina_lgn.params.gain' : [1.3],'l4_cortex_inh.ExcInhAfferentRatio' : [1.3],'l4_cortex_exc.AfferentConnection.base_weight' : [0.0025],'l4_cortex_exc.AfferentConnection.num_samples' : [10,15],'retina_lgn.params.noise.stdev' : [3.5,4.0,5.0],'l4_cortex_inh.L4InhL4InhConnection.base_weight' : [0.0]}).run_parameter_search()

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.rand_struct_ratio' : [0.6,0.3], 'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.0],'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0003,0.0006], 'l4_cortex_inh.params.cell.params.v_thresh' : [-53,-54],'retina_lgn.params.gain' : [1.3],'l4_cortex_inh.ExcInhAfferentRatio' : [1.0,1.2],'l4_cortex_exc.AfferentConnection.base_weight' : [0.0025,0.0015],'l4_cortex_exc.AfferentConnection.num_samples' : [10],'retina_lgn.params.noise.stdev' : [2.7,3.5],'l4_cortex_inh.L4InhL4InhConnection.base_weight' : [0.0]}).run_parameter_search()

CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0015,0.002],'l4_cortex_inh.L4InhL4ExcConnection.short_term_plasticity.tau_rec' : [200,100],'l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma': [0.3],'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0006],'l4_cortex_exc.params.cell.params.tau_syn_I' : [15,17],'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.0002],'retina_lgn.params.gain' : [1.7,2.0,2.5],'l4_cortex_inh.params.cell.params.tau_refrac' : [1.0],'l4_cortex_exc.AfferentConnection.num_samples' : [25]}).run_parameter_search()


#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.params.cell.params.tau_syn_I' : [12,13,14],'retina_lgn.params.gain' : [1.1,1.3,1.4],'retina_lgn.params.noise.stdev' : [2.7,3.0],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.001,0.0008,0.0006]}).run_parameter_search()

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma': [0.4],'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0004,0.0006],'l4_cortex_exc.params.cell.params.tau_syn_I' : [15],'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.0002],'retina_lgn.params.gain' : [1.1,1.3],'l4_cortex_inh.params.cell.params.tau_refrac' : [1.0,2.0],'l4_cortex_exc.AfferentConnection.num_samples' : [20,30,40]}).run_parameter_search()



#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.0001,0.00015,0.0002,0.00025,0.0003],'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0002,0.00025,0.0003,0.00035,0.0004]}).run_parameter_search()



#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.params.density' : [12000],'l4_cortex_exc.rand_struct_ratio' : [0.8,0.7], 'l4_cortex_inh.ExcInhAfferentRatio' : [1.2],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.00014,0.0002,0.00025],'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0005,0.0006,0.0007]}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.00055],'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.001],'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0006], 'l4_cortex_inh.ExcInhAfferentRatio' : [1.0],'retina_lgn.params.gain' : [2.0],'l4_cortex_exc.rand_struct_ratio' : [0.6,0.8,0.95], 'l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [0.9],'l4_cortex_exc.params.cell.params.v_reset' : [-51],'l4_cortex_exc.params.cell.params.tau_w' : [100,200,400],'l4_cortex_exc.params.cell.params.b' : [0.4,0.6,1.0,2,4,10]}).run_parameter_search()