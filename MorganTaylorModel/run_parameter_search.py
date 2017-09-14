# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackendIoV
import numpy
import time

if True:
    CombinationParameterSearch(SlurmSequentialBackendIoV(num_threads=6,num_mpi=1),{
									     'sheets.l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0021],
									     'sheets.l4_cortex_exc.inhibitory_connection_ratio' : [0.75],
									     'sheets.l4_cortex_exc.params.cell.params.tau_syn_I' : [1.9],
									     'sheets.l4_cortex_exc.params.cell.params.v_thresh' : [-56],
									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.tau_rec' : [30],
									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [1.6],
#									     'sheets.l4_cortex_exc.L4ExcL4InhConnection.weight_functions.f1.params.sigma' : [3.0],
#									     'sheets.l4_cortex_exc.AfferentConnection.off_bias' : [1.0],
									     'sheets.retina_lgn.params.gain_control.gain' : [20],
#									     'sheets.retina_lgn.params.gain_control.non_linear_gain' : [None],
#									     'sheets.retina_lgn.params.receptive_field.func_params.sigma_s' : [1.0,0.9,0.8,0.7],
#									     'sheets.retina_lgn.params.receptive_field.func_params.sigma_c' : [0.4,0.3,0.2],
									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [0.4,0.6,0.8,1.0,1.2,1.6,2.0,2.4,2.8,3.2,3.6,4.0],
									     'sheets.l4_cortex_exc.L4ExcL4InhConnection.weight_functions.f1.params.sigma' : [0.4,0.6,0.8,1.0,1.2,1.6,2.0,2.4,2.8,3.2,3.6,4.0],
									     'sheets.retina_lgn.params.gain_control.non_linear_gain.contrast_scaler' : [200000],
									     'sheets.retina_lgn.params.gain_control.non_linear_gain.luminance_gain' : [0.0],
									     'sheets.l4_cortex_exc.AfferentConnection.size' : [0.15],
#									     'only_afferent' : [True,False],
#									     'sheets.retina_lgn.params.cached' : [True,False],
#									     'sheets.retina_lgn.params.gain_control.non_linear_gain.luminance_scaler' : [0.1,0.5],
									     'sheets.retina_lgn.params.noise.stdev' : [2.2],
									     'sheets.l4_cortex_exc.params.density' : [4000],
									     }).run_parameter_search()


if False:
    CombinationParameterSearch(SlurmSequentialBackendIoV(num_threads=8,num_mpi=1),{
									     'sheets.l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0024],
									     'sheets.l4_cortex_exc.params.cell.params.tau_syn_I' : [3.0],
									     'sheets.l4_cortex_exc.params.cell.params.v_thresh' : [-56],
									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.tau_rec' : [70],
									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [0.4,0.8,1.2,1.4,1.6,2.0,2.4,2.8,3.2,3.6,4.0],
									     'sheets.l4_cortex_exc.L4ExcL4InhConnection.weight_functions.f1.params.sigma' : [0.4,0.8,1.2,1.4,1.6,2.0,2.4,2.8,3.2,3.6,4.0],
									     }).run_parameter_search()



if False:
    CombinationParameterSearch(SlurmSequentialBackendIoV(num_threads=8,num_mpi=1),{
									     'sheets.l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0025],
									     'sheets.l4_cortex_exc.params.cell.params.tau_syn_I' : [3.0],
									     'sheets.l4_cortex_exc.params.cell.params.v_thresh' : [-56],
									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.tau_rec' : [70],
									     'input_space.update_interval' : [1.0,15.0],
#									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.weight_functions.f1.params.sigma' : [2.0,2.4,2.8,3.2,3.6,4.0],
#									     'sheets.l4_cortex_exc.L4ExcL4InhConnection.weight_functions.f1.params.sigma' : [4.0],
									     }).run_parameter_search()






