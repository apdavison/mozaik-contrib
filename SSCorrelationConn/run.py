# -*- coding: utf-8 -*-
"""
"""
import sys
from pyNN import nest
import mozaik.controller
from mozaik.controller import run_workflow, setup_logging
import mozaik
from experiments import create_experiments
from model import SSCorrelationConnectivity
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

if True:
    data_store,model = run_workflow('SSCorrelationConnectivity',SSCorrelationConnectivity,create_experiments)
    data_store.save() 
else: 
    setup_logging()
    data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':'/home/jan/cluster/dev/pkg/mozaik/mozaik/contrib/SSCorrelationConn/20140313-113338[param_sd.defaults]CombinationParamSearch{7}/SSCorrelationConnectivity_ParameterSearch_____base_weight:0.00045_sigma:0.5_base_weight:0.0007_rand_struct_ratio:0.5_ExcInhAfferentRatio:1.0_base_weight:0.0007_gain:15.0', 'store_stimuli' : False}),replace=True)

if mpi_comm.rank == 0:
   print "Starting visualization" 
   perform_analysis_and_visualization(data_store)
#   data_store.save() 
