# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackend,SlurmSequentialBackendIoV
import numpy
import time




if True:
    CombinationParameterSearch(SlurmSequentialBackendIoV(num_threads=45,num_mpi=1),{
									     'sheets.l23_cortex_exc.L23ExcL4InhConnection.base_weight' : [0.0017],
									     'sheets.l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0026],
#									     'sheets.l4_cortex_exc.inhibitory_connection_ratio' : [0.7],
									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.tau_rec' : [100],
									     'sheets.l4_cortex_exc.params.cell.params.tau_syn_I' : [1.7],
#									     'sheets.l4_cortex_exc.K' : [1480],
#									     'sheets.l4_cortex_exc.params.density' : [1400,2000],
									     'sheets.l4_cortex_exc.params.sx' : [3200,3500],
#									     'sheets.l23_cortex_exc.feedback_conn_ratio' : [0.13],
#									     'sheets.l4_cortex_exc.layer23_aff_ratio' : [0.35],
									     'sheets.retina_lgn.params.noise.stdev' : [2.7],
									     'sheets.retina_lgn.params.gain_control.gain' : [30],
#									     'sheets.l4_cortex_exc.AfferentConnection.short_term_plasticity.tau_rec' : [250],
									     'sheets.l4_cortex_exc.AfferentConnection.base_weight' : [0.0016],
#				        				     'sheets.retina_lgn.params.gain_control.non_linear_gain.contrast_scaler' : [0,1,10],
#				        				     'sheets.retina_lgn.params.gain_control.non_linear_gain.luminance_gain' : [0.0],
#									     'sheets.l4_cortex_exc.AfferentConnection.size' : [0.15],
#									     'sheets.l23_cortex_exc.feedback_conn_ratio' : [0.13],
									     'sheets.l4_cortex_inh.ExcInhAfferentRatio' : [1.1],
#									     'sheets.l4_cortex_exc.params.cell.params.v_thresh' : [-57],
#									     'sheets.l4_cortex_exc.params.cell.params.v_rest' : [-71],
									     'sheets.l23_cortex_inh.L4ExcL23InhConnection.num_samples' : [110,115,120,125],
									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [1.4],
#									     'sheets.l4_cortex_exc.L4ExcL4InhConnection.weight_functions.f1.params.sigma' : [3.0],
#									     'sheets.l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f3.params.arborization_scaler' : [1.0,4.0],
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



