# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackend
import numpy
CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=1),{'l4_cortex_exc.params.artificial_stimulators.background_act.params.exc_firing_rate' : numpy.linspace(0.0,100000,400)}).run_parameter_search()

