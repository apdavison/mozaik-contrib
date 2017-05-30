# -*- coding: utf-8 -*-
"""

"""
from mpi4py import MPI 
from pyNN import nest
import sys
import mozaik.controller
from mozaik.controller import run_workflow, setup_logging
import mozaik
from experiments import create_experiments
from model import SelfSustainedPushPull
from mozaik.storage.datastore import Hdf5DataStore,PickledDataStore
from analysis_and_visualization import perform_analysis_and_visualization
from parameters import ParameterSet


mpi_comm = MPI.COMM_WORLD

if True:
    data_store,model = run_workflow('TestLGN',SelfSustainedPushPull,create_experiments)
    data_store.save() 
else: 
    #setup_logging()
    data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':'TestLGN_rup=28;rfr=7_____','store_stimuli' : False}),replace=True)

if mpi_comm.rank == 0:
   print "Starting visualization" 
   perform_analysis_and_visualization(data_store)
#   data_store.save() 
