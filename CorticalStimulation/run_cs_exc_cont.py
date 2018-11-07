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
from experiments import create_experiments_cortical_stimulation_cont_exc
from model import SelfSustainedPushPull
from mozaik.storage.datastore import Hdf5DataStore,PickledDataStore
from analysis_and_visualization import perform_analysis_and_visualization
from parameters import ParameterSet


mpi_comm = MPI.COMM_WORLD

data_store,model = run_workflow('CorticalStimulationModel',SelfSustainedPushPull,create_experiments_cortical_stimulation_cont_exc)
data_store.save() 

if mpi_comm.rank == 0:
   print "Starting visualization" 
   perform_analysis_and_visualization(data_store,gratings=False,cort_stim=True,nat_stim=False,tp=1)
