# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackendIoV
import numpy
import time

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






