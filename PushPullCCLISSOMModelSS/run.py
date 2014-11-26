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
from experiments import create_experiments
from model import SelfSustainedPushPull
from mozaik.storage.datastore import Hdf5DataStore,PickledDataStore
from analysis_and_visualization import perform_analysis_and_visualization
from parameters import ParameterSet

from mpi4py import MPI 
mpi_comm = MPI.COMM_WORLD


if True:
    data_store,model = run_workflow('SelfSustainedPushPull',SelfSustainedPushPull,create_experiments)
    model.connectors['V1AffConnectionOn'].store_connections(data_store)    
    model.connectors['V1AffConnectionOff'].store_connections(data_store)    
    model.connectors['V1AffInhConnectionOn'].store_connections(data_store)    
    model.connectors['V1AffInhConnectionOff'].store_connections(data_store)    
    model.connectors['V1L4ExcL4ExcConnection'].store_connections(data_store)    
    model.connectors['V1L4ExcL4InhConnection'].store_connections(data_store)    
    model.connectors['V1L4InhL4ExcConnection'].store_connections(data_store)    
    model.connectors['V1L4InhL4InhConnection'].store_connections(data_store)    
    model.connectors['V1L4ExcL4ExcConnectionRand'].store_connections(data_store)    
    model.connectors['V1L4ExcL4InhConnectionRand'].store_connections(data_store)    
    model.connectors['V1L4InhL4ExcConnectionRand'].store_connections(data_store)    
    model.connectors['V1L4InhL4InhConnectionRand'].store_connections(data_store)    
    data_store.save() 
    
else: 
    setup_logging()
    data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':'SelfSustainedPushPull_BIG_OR____','store_stimuli' : False}),replace=True)

if mpi_comm.rank == 0:
   print "Starting visualization" 
   perform_analysis_and_visualization(data_store)
#   data_store.save() 
