# -*- coding: utf-8 -*-
"""

"""
import matplotlib
matplotlib.use('Agg')
from mpi4py import MPI 
#from pyNN import nest
import sys
import mozaik.controller
from mozaik.controller import run_workflow, setup_logging
import mozaik
from experiments import create_experiments_short,create_experiments_old,create_experiments,create_experiments_tmp
from model import SelfSustainedPushPull
from mozaik.storage.datastore import Hdf5DataStore,PickledDataStore
from analysis_and_visualization import perform_analysis_and_visualization
from parameters import ParameterSet

mpi_comm = MPI.COMM_WORLD

if True:
    data_store,model = run_workflow('CorticalStimulationModel',SelfSustainedPushPull,create_experiments_tmp)
    data_store.save() 
else: 
    setup_logging()
    data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':'CorticalStimulationModel_visual_stimulation_____base_weight:0.0022_inhibitory_connection_ratio:0.5_layer23_aff_ratio:0.4_stdev:2.7','store_stimuli' : False}),replace=True)

if mpi_comm.rank == 0:
   print "Starting visualization" 
   perform_analysis_and_visualization(data_store,gratings=True,cort_stim=False,nat_stim=False,tp=0)
