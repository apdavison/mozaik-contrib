# -*- coding: utf-8 -*-
"""

"""
import matplotlib
matplotlib.use('Agg')
from mpi4py import MPI 
from pyNN import nest
import sys
import mozaik.controller
from mozaik.controller import run_workflow, setup_logging
import mozaik
from experiments import create_experiments,create_experiments_bar,create_experiments_short,create_experiments_old,create_experiments_old_short
from model import SelfSustainedPushPull
from mozaik.storage.datastore import Hdf5DataStore,PickledDataStore
from analysis_and_visualization import perform_analysis_and_visualization
from parameters import ParameterSet


mpi_comm = MPI.COMM_WORLD

data_store,model = run_workflow('MorganTaylorModel',SelfSustainedPushPull,create_experiments_old)
data_store.save() 

if mpi_comm.rank == 0:
   print "Starting visualization" 
   perform_analysis_and_visualization(data_store,gratings=True,bars=False,nat_movies=True)

