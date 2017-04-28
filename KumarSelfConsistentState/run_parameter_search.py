# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackend
import numpy
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=1),{'l4_cortex_exc.params.artificial_stimulators.background_act.params.exc_firing_rate' : numpy.linspace(0.0,5000,400)}).run_parameter_search()
#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=1),{'l4_cortex_exc.params.artificial_stimulators.background_act.params.exc_firing_rate' : numpy.linspace(0.0,100000,100),'l4_cortex_exc.params.artificial_stimulators.background_act.params.exc_weight' : numpy.linspace(0.0,0.04,10),'l4_cortex_exc.params.artificial_stimulators.background_act.params.inh_weight' : numpy.linspace(0.0,0.2,10)}).run_parameter_search()


#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=1),{'l4_cortex_exc.params.artificial_stimulators.background_act.params.exc_firing_rate' : numpy.linspace(0.0,10000,400),'l4_cortex_exc.params.artificial_stimulators.background_act.params.exc_weight' : [0.0006,0.00062,0.00064,0.00066],'l4_cortex_exc.params.artificial_stimulators.background_act.params.inh_weight' : [0.0011]}).run_parameter_search()

#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=1),{'l4_cortex_exc.params.artificial_stimulators.background_act.params.exc_firing_rate' : numpy.linspace(0.0,500000,100),'l4_cortex_exc.params.artificial_stimulators.background_act.params.inh_firing_rate' : numpy.linspace(0.0,400000,100)}).run_parameter_search()


#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=1),{'l4_cortex_exc.params.artificial_stimulators.background_act.params.exc_firing_rate' : numpy.linspace(0.0,100000,50),'l4_cortex_exc.params.artificial_stimulators.background_act.params.inh_firing_rate' : numpy.linspace(0.0,100000,50),'l4_cortex_exc.params.artificial_stimulators.background_act.params.exc_weight' : numpy.linspace(0.0006,0.0008,10)}).run_parameter_search()



#CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=1),{'l4_cortex_exc.params.artificial_stimulators.background_act.params.exc_firing_rate' : numpy.linspace(0.0,10000,30),'l4_cortex_exc.params.cell.params.tau_syn_E' : [2.5,3.0,3.5,4.0], 'l4_cortex_exc.params.cell.params.tau_syn_I':  [4.0,5.0,6.0,7.0] }).run_parameter_search()

CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=1),{'l4_cortex_exc.params.artificial_stimulators.background_act.params.exc_firing_rate' : numpy.linspace(0.0,1000,30),'l4_cortex_exc.params.cell.params.cm' : [0.1,0.05], 'l4_cortex_exc.params.cell.params.tau_m':  [5.0,10.0,20.0,30.0] }).run_parameter_search()