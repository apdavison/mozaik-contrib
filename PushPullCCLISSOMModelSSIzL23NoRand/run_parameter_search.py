# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackend
import numpy
import time



if False:
    CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{
									     'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0017],
#									     'l4_cortex_exc.inhibitory_connection_ratio' : [0.7],
#									     'l4_cortex_inh.params.cell.params.tau_m' : [10.0],
#									     'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.0014],
#									     'l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.tau_rec' : [50],
#									     'l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.U' : [0.75],
									     'l4_cortex_exc.params.sx' : [3500],
									     'l4_cortex_exc.params.cell.params.tau_syn_E' : [1.5,1.3,1.1,0.9],
									     'l4_cortex_exc.params.cell.params.tau_syn_I' : [2.7,2.5,2.3,2.1,1.9,1.7,1.5,1.3],
									     'l4_cortex_exc.K' : [1600],
									     'l4_cortex_exc.layer23_aff_ratio' : [0.5],
									     'l4_cortex_exc.params.density' : [2000],
#									     'l23_cortex_exc.params.density' : [2000],
									     'retina_lgn.params.noise.stdev' : [2.7],
									     'l4_cortex_exc.params.cell.params.v_rest' : [-72],
									     'retina_lgn.params.gain_control.gain' : [50],
#									     'l23_cortex_exc.L23ExcL4ExcConnection.weight_functions.f1.params.arborization_constant' : [30,50],
#									     'l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f3.params.arborization_scaler' : [4.0],
#									     'l23_cortex_exc.L4ExcL23ExcConnection.weight_functions.f2.params.sigma' : [0.1,0.3],
#									     'l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f1.params.sigma' : [0.3,0.5],
									     'l4_cortex_inh.L4InhL4ExcConnection.weight_functions.f1.params.sigma' : [3.0],
									     'l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [0.5],
#									     'l4_cortex_inh.L4InhL4InhConnection.base_weight' : [0.001,0.0006],
									     'l4_cortex_inh.AfferentConnection.delay' : [2.0],
									     'l4_cortex_exc.AfferentConnection.size' : [0.25],
#									     'feedback' : [False],
									     'l23' : [False],
#								     	     'l23_cortex_exc.L23ExcL4ExcConnection.num_samples' : [20,30],
									     }).run_parameter_search()






if False:
    CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{
									     'sheets.l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0020,0.0021,0.0022],
									     'sheets.l4_cortex_exc.inhibitory_connection_ratio' : [0.75],
									     'sheets.l4_cortex_inh.params.cell.params.tau_m' : [10.0],
#									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.0014],
									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.tau_rec' : [50],
#									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.U' : [0.75],
#									     'sheets.l4_cortex_exc.params.sx' : [3500],
#									     'sheets.l4_cortex_exc.params.cell.params.tau_syn_E' : [1.1],
									     'sheets.l4_cortex_exc.params.cell.params.tau_syn_I' : [2.5],
#									     'sheets.l4_cortex_exc.params.cell.params.delta_T' : [2.0],
									     'sheets.l4_cortex_exc.K' : [1000],
									     'sheets.l4_cortex_exc.layer23_aff_ratio' : [0.35],
#									     'sheets.l23_cortex_exc.feedback_conn_ratio' : [0.13],
#									     'sheets.l4_cortex_exc.params.density' : [3000],
#									     'sheets.l23_cortex_exc.params.density' : [2000],
									     'sheets.retina_lgn.params.noise.stdev' : [2.5],
#									     'sheets.l4_cortex_inh.params.cell.params.v_thresh' : [-52],
#									     'sheets.l4_cortex_exc.params.cell.params.v_rest' : [-72],
									     'sheets.l4_cortex_exc.params.cell.params.v_thresh' : [-56],
#									     'sheets.l4_cortex_inh.params.cell.params.v_thresh' : [-52],
#									     'sheets.l23_cortex_exc.params.cell.params.v_thresh' : [-58],
#									     'sheets.l23_cortex_exc.params.cell.params.v_thresh' : [-60,-62],
#									     'sheets.l4_cortex_inh.ExcInhAfferentRatio' : [1.1,1.15],
#									     'sheets.retina_lgn.params.gain_control.gain' : [50],
#									     'sheets.l23_cortex_exc.L23ExcL4InhConnection.weight_functions.f1.params.arborization_constant' : [100],
#									     'sheets.l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f3.params.arborization_scaler' : [4.0],
#									     'sheets.l23_cortex_exc.L4ExcL23ExcConnection.weight_functions.f2.params.sigma' : [0.075,0.3],
#									     'sheets.l23_cortex_inh.L4ExcL23InhConnection.weight_functions.f2.params.sigma' : [2.5],
									     'sheets.l23_cortex_inh.L4ExcL23InhConnection.num_samples' : [140],
#									     'sheets.l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f1.params.sigma' : [0.3,2.5],
#									     'sheets.l4_cortex_inh.L4InhL4ExcConnection.weight_functions.f1.params.sigma' : [2.5],
									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [1.4],
#									     'sheets.l4_cortex_exc.L4ExcL4InhConnection.weight_functions.f1.params.sigma' : [3.0],
#									     'sheets.l4_cortex_inh.L4InhL4InhConnection.base_weight' : [0.001,0.0006],
#									     'sheets.l4_cortex_inh.AfferentConnection.delay' : [1.2],
#									     'sheets.l4_cortex_exc.AfferentConnection.size' : [0.15,0.5],
									     'feedback' : [True],
									     'l23' : [True],
									     }).run_parameter_search()



if False:
    CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=64),{
									     'sheets.l23_cortex_exc.L23ExcL4InhConnection.base_weight' : [0.0016],
									     'sheets.l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0023],
#									     'sheets.l4_cortex_exc.inhibitory_connection_ratio' : [0.8],
#									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.tau_rec' : [50,100],
									     'sheets.l4_cortex_exc.params.cell.params.tau_syn_I' : [2.4,2.5],
#									     'sheets.l4_cortex_exc.K' : [1480],
#									     'l23' : [False],
									     'sheets.l4_cortex_exc.params.density' : [2000],
#									     'sheets.l4_cortex_exc.params.sx' : [3500],
#									     'sheets.l23_cortex_exc.feedback_conn_ratio' : [0.13],
									     'sheets.l4_cortex_exc.layer23_aff_ratio' : [0.3],
									     'sheets.retina_lgn.params.noise.stdev' : [2.3],
									     'sheets.retina_lgn.params.gain_control.gain' : [20,30,40],
#									     'sheets.l4_cortex_exc.AfferentConnection.short_term_plasticity.tau_rec' : [100,150],
				        				     'sheets.retina_lgn.params.gain_control.non_linear_gain.contrast_scaler' : [0.1],
				        				     'sheets.retina_lgn.params.gain_control.non_linear_gain.luminance_gain' : [0.0],
#									     'sheets.l4_cortex_exc.AfferentConnection.size' : [0.15,0.2],
#									     'sheets.l23_cortex_exc.feedback_conn_ratio' : [0.13],
									     'sheets.l4_cortex_exc.params.cell.params.v_thresh' : [-57],
									     'sheets.l23_cortex_inh.L4ExcL23InhConnection.num_samples' : [125],
									     'sheets.l4_cortex_inh.ExcInhAfferentRatio' : [1.0,1.05,1.1],
#									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [1.2,1.4],
									     'sheets.l4_cortex_exc.L4ExcL4InhConnection.weight_functions.f1.params.sigma' : [2.0,3.0],
#									     'sheets.l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f3.params.arborization_scaler' : [3.0],
#									     'sheets.l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f1.params.sigma' : [1.0,1.4],
#									     'sheets.l23_cortex_inh.L4ExcL23InhConnection.weight_expression' : ["\\\'f1\\\'","\\\'f1*f2\\\'"],
#									     'sheets.l23_cortex_exc.L23ExcL4InhConnection.weight_functions.f1.params.arborization_constant' : [100],
									     }).run_parameter_search()



if True:
    CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=32),{
									     'sheets.l23_cortex_exc.L23ExcL4InhConnection.base_weight' : [0.002],
									     'sheets.l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0023,0.0024],
									     'sheets.l4_cortex_exc.inhibitory_connection_ratio' : [0.75],
									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.tau_rec' : [100,150],
#									     'sheets.l4_cortex_exc.params.cell.params.tau_syn_I' : [1.9],
#									     'sheets.l4_cortex_exc.K' : [1480],
#									     'sheets.l4_cortex_exc.params.density' : [2000],
#									     'sheets.l4_cortex_exc.params.sx' : [3500],
#									     'sheets.l23_cortex_exc.feedback_conn_ratio' : [0.13],
									     'sheets.l4_cortex_exc.layer23_aff_ratio' : [0.3],
									     'sheets.retina_lgn.params.noise.stdev' : [2.4],
									     'sheets.retina_lgn.params.gain_control.gain' : [40,50],
									     'sheets.l4_cortex_exc.AfferentConnection.short_term_plasticity.tau_rec' : [200,250],
#				        				     'sheets.retina_lgn.params.gain_control.non_linear_gain.contrast_scaler' : [0,1,10],
#				        				     'sheets.retina_lgn.params.gain_control.non_linear_gain.luminance_gain' : [0.0],
#									     'sheets.l4_cortex_exc.AfferentConnection.size' : [0.15,0.2],
#									     'sheets.l23_cortex_exc.feedback_conn_ratio' : [0.13],
									     'sheets.l4_cortex_inh.ExcInhAfferentRatio' : [1.0,1.1],
									     'sheets.l4_cortex_exc.params.cell.params.v_thresh' : [-57],
									     'sheets.l23_cortex_inh.L4ExcL23InhConnection.num_samples' : [115,125],
#									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [1.2,1.4],
#									     'sheets.l4_cortex_exc.L4ExcL4InhConnection.weight_functions.f1.params.sigma' : [2.0,1.4],
#									     'sheets.l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f3.params.arborization_scaler' : [3.0],
#									     'sheets.l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f1.params.sigma' : [1.0,1.4],
#									     'sheets.l23_cortex_inh.L4ExcL23InhConnection.weight_expression' : ["\\\'f1\\\'","\\\'f1*f2\\\'"],
#									     'sheets.l23_cortex_exc.L23ExcL4InhConnection.weight_functions.f1.params.arborization_constant' : [100],
									     }).run_parameter_search()


#import time
#time.sleep(5) # delays for 5 seconds

#if True:
#    CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{
#									     'sheets.l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0020,0.0021,0.0022,0.0023],
#									     'sheets.l4_cortex_exc.inhibitory_connection_ratio' : [0.75],
#									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.tau_rec' : [50],
#									     'sheets.l4_cortex_exc.params.cell.params.tau_syn_I' : [2.5],
#									     'sheets.l4_cortex_exc.K' : [1000],
#									     'sheets.l4_cortex_exc.layer23_aff_ratio' : [0.35],
#									     'sheets.retina_lgn.params.noise.stdev' : [2.3,2.5],
#									     'sheets.l4_cortex_exc.params.cell.params.v_thresh' : [-56],
#									     'sheets.l23_cortex_inh.L4ExcL23InhConnection.num_samples' : [140],
#									     }).run_parameter_search()








