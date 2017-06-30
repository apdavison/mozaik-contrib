# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackendIoV
import numpy
import time

CombinationParameterSearch(SlurmSequentialBackendIoV(num_threads=10,num_mpi=1),{
									     'sheets.l4_cortex_exc.K' : [1000],
									     'sheets.l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0019,0.002],
									     'sheets.l4_cortex_exc.params.cell.params.tau_syn_I' : [1.9,2.1],
									     'sheets.retina_lgn.params.noise.stdev' : [3.1],
									     'sheets.l4_cortex_exc.layer23_aff_ratio' : [0.3,0.4],
									     'sheets.l4_cortex_exc.inhibitory_connection_ratio' : [0.6,0.7],
									     'sheets.l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.tau_rec' : [50,70],
									     'sheets.l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f3.params.arborization_scaler' : [4.0],
									     }).run_parameter_search()






