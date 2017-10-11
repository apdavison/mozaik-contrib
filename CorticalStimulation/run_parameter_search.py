# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackendIoV
import numpy
import time

if True:
    CombinationParameterSearch(SlurmSequentialBackendIoV(num_threads=12,num_mpi=1),{
#									     'sheets.l23_cortex_exc.L23ExcL4InhConnection.base_weight' : [0.002],
#									     'sheets.l4_cortex_exc.layer23_aff_ratio' : [0.3,0.35],
#									     'sheets.l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0024],
#									     'sheets.l23_cortex_inh.L4ExcL23InhConnection.num_samples' : [125],
#									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.tau_rec' : [100,200],
#									     'feedback' : [False],
#									     'sheets.l4_cortex_exc.params.cell.params.tau_syn_I' : [1.9,2.4,3.0],
									     'with_cortical_conn' : [True,False],
#									     'sheets.l4_cortex_exc.K' : [1000,1480],
#									     'sheets.retina_lgn.params.noise.stdev' : [2.1],
									     }).run_parameter_search()

if False:
    CombinationParameterSearch(SlurmSequentialBackendIoV(num_threads=10,num_mpi=1),{
									     'sheets.l4_cortex_exc.K' : [1480],
									     'sheets.l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0021,0.0022,0.0023,0.0024],
									     'sheets.l4_cortex_exc.params.cell.params.tau_syn_I' : [2.1],
									     'sheets.retina_lgn.params.noise.stdev' : [2.8],
									     'sheets.l4_cortex_exc.layer23_aff_ratio' : [0.3],
									     'sheets.l4_cortex_exc.inhibitory_connection_ratio' : [0.6,0.7],
									     'sheets.l23_cortex_inh.L4ExcL23InhConnection.num_samples' : [135],
									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.tau_rec' : [100],
									     'sheets.l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f3.params.arborization_scaler' : [1.0,3.0],
									     'sheets.l23_cortex_inh.L4ExcL23InhConnection.weight_expression' : ["\\\'f1\\\'","\\\'f1*f2\\\'"],
									     'sheets.l4_cortex_exc.params.density' : [2000],

									     }).run_parameter_search()






if False:
    CombinationParameterSearch(SlurmSequentialBackendIoV(num_threads=12,num_mpi=1),{
				         'sheets.l23_cortex_exc.L23ExcL4InhConnection.base_weight' : [0.002],
				         'sheets.l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0023],
				         'sheets.l4_cortex_exc.inhibitory_connection_ratio' : [0.75],
				         'sheets.l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.tau_rec' : [100],
				         'sheets.l4_cortex_exc.params.cell.params.tau_syn_I' : [1.9],
				         'sheets.l4_cortex_exc.layer23_aff_ratio' : [0.3],
				         'sheets.retina_lgn.params.noise.stdev' : [2.4],
				         'sheets.retina_lgn.params.gain_control.gain' : [40],
				         'sheets.l4_cortex_exc.AfferentConnection.short_term_plasticity.tau_rec' : [200],
				         'sheets.l4_cortex_exc.AfferentConnection.size' : [0.2],
				         'sheets.l4_cortex_inh.ExcInhAfferentRatio' : [1.1],
				         'sheets.l4_cortex_exc.params.cell.params.v_thresh' : [-57],
				         'sheets.l23_cortex_inh.L4ExcL23InhConnection.num_samples' : [115],
				         }).run_parameter_search()
