# -*- coding: utf-8 -*-
"""
"""
import sys
from pyNN import nest
from mozaik.controller import run_workflow, setup_logging
import mozaik
from model import TestStimuliModel
from experiments import create_experiments
from mozaik.storage.datastore import Hdf5DataStore,PickledDataStore
from analysis_and_visualization import perform_analysis_and_visualization
from parameters import ParameterSet


try:
    from mpi4py import MPI
except ImportError:
    MPI = None
if MPI:
    mpi_comm = MPI.COMM_WORLD
MPI_ROOT = 0

logger = mozaik.getMozaikLogger()

if True:
    data_store,model = run_workflow('TestStimuliModel',TestStimuliModel,create_experiments)
    data_store.save()
    
else: 
    setup_logging()
    data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':'TestStimuliModel_test_____', 'store_stimuli' : False}),replace=True)
    logger.info('Loaded data store')
    data_store.save()

if mpi_comm.rank == MPI_ROOT:
    perform_analysis_and_visualization(data_store)
