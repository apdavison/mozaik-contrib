# -*- coding: utf-8 -*-
"""
This is implementation of model of self-sustained activitity in balanced networks from: 
Vogels, T. P., & Abbott, L. F. (2005). 
Signal propagation and logic gating in networks of integrate-and-fire neurons. 
The Journal of neuroscience : the official journal of the Society for Neuroscience, 25(46), 10786–95. 
"""
from pyNN import nest
import sys
import mozaik.controller
from mozaik.controller import run_workflow, setup_logging
import mozaik
from experiments import create_experiments_spont
from model import SelfSustainedPushPull
from mozaik.storage.datastore import Hdf5DataStore,PickledDataStore
from analysis_and_visualization import perform_analysis_and_visualization_spont
from parameters import ParameterSet

from mpi4py import MPI 
mpi_comm = MPI.COMM_WORLD

data_store,model = run_workflow('SelfSustainedPushPull',SelfSustainedPushPull,create_experiments_spont)
data_store.save() 

perform_analysis_and_visualization_spont(data_store)