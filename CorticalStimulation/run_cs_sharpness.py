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
from experiments import create_experiments_cortical_stimulation_sharpness
from model import SelfSustainedPushPull
from mozaik.storage.datastore import Hdf5DataStore,PickledDataStore
from analysis_and_visualization import perform_analysis_and_visualization
from parameters import ParameterSet


mpi_comm = MPI.COMM_WORLD

if True:
    data_store,model = run_workflow('CorticalStimulationModel',SelfSustainedPushPull,create_experiments_cortical_stimulation_sharpness)
    data_store.save() 
else: 
    setup_logging()
    data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':'CorticalStimulationModel_visual_stimulation_full_protocol_____','store_stimuli' : False}),replace=True)

if mpi_comm.rank == 0:
   print "Starting visualization" 
   perform_analysis_and_visualization(data_store,gratings=False,cort_stim=True,nat_stim=False,tp=2)
