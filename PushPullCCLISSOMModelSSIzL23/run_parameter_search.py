# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackend
import numpy


CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{'l4_cortex_exc.rand_struct_ratio' : [0.7,0.9],'l23_cortex_exc.L4ExcL23ExcConnection.base_weight' : [0.0007,0.0001],'l23_cortex_inh.L4ExcL23InhConnection.base_weight' : [0.0003,0.0007,0.001]}).run_parameter_search()


