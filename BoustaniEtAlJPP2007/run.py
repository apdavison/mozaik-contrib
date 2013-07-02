# -*- coding: utf-8 -*-
"""
This is implementation of model of self-sustained activitity in balanced networks from: 

El Boustani, S., Pospischil, M., Rudolph-Lilith, M., & Destexhe, A. (n.d.). 
Activated cortical states: experiments, analyses and models. 
Journal of physiology, Paris, 101(1-3), 99–109.
"""
from pyNN import nest
import sys
from mozaik.controller import run_workflow, setup_logging
import mozaik
from experiments import create_experiments
from model import Boustani2007
from mozaik.storage.datastore import Hdf5DataStore,PickledDataStore
from analysis_and_visualization import perform_analysis_and_visualization

from mpi4py import MPI 
mpi_comm = MPI.COMM_WORLD


if True:
    logger = mozaik.getMozaikLogger()
    data_store,model = run_workflow('Boustani2007',Boustani2007,create_experiments)
else: 
    setup_logging()
    data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':'A'}),replace=True)
    logger.info('Loaded data store')

if mpi_comm.rank == 0:
   print "Starting visualization" 
   perform_analysis_and_visualization(data_store)
   data_store.save() 
