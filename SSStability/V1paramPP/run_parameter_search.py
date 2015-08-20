# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackend
import numpy

# VogelsAbbott run
CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=2),{'exc1.Exc1Exc1Connection.weights' : numpy.linspace(0.0,0.002,20),'inh1.Inh1Exc1Connection.weights' : numpy.linspace(0.0,0.02,20)}).run_parameter_search()