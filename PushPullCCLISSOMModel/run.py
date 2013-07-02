# -*- coding: utf-8 -*-
"""
This is implementation of model of push-pull connectvity: 
Jens Kremkow: Correlating Excitation and Inhibition in Visual Cortical Circuits: Functional Consequences and Biological Feasibility. PhD Thesis, 2009.
"""
import sys
from pyNN import nest
from mozaik.controller import run_workflow, setup_logging
import mozaik
from model import PushPullCCModel
from experiments import create_experiments
from mozaik.storage.datastore import Hdf5DataStore,PickledDataStore
from analysis_and_visualization import perform_analysis_and_visualization
from parameters import ParameterSet

logger = mozaik.getMozaikLogger()

if True:
    data_store,model = run_workflow('FFI',PushPullCCModel,create_experiments)
    #jens_model.connectors['ON_to_[V1_Exc_L4]'].store_connections(data_store)    
    #jens_model.connectors['OFF_to_[V1_Exc_L4]'].store_connections(data_store)    
    #jens_model.connectors['ON_to_[V1_Inh_L4]'].store_connections(data_store)    
    #jens_model.connectors['OFF_to_[V1_Inh_L4]'].store_connections(data_store)    
    #jens_model.connectors['V1L4ExcL4ExcConnection'].store_connections(data_store)    
    #jens_model.connectors['V1L4ExcL4InhConnection'].store_connections(data_store)    
    #jens_model.connectors['V1L4InhL4ExcConnection'].store_connections(data_store)    
    #jens_model.connectors['V1L4InhL4InhConnection'].store_connections(data_store)    
    #jens_model.connectors['V1ExcL23ExcL23Connection'].store_connections(data_store)    
    #jens_model.connectors['V1ExcL23InhL23Connection'].store_connections(data_store)    
    #jens_model.connectors['V1InhL23ExcL23Connection'].store_connections(data_store)    
    #jens_model.connectors['V1InhL23InhL23Connection'].store_connections(data_store)    
    #jens_model.connectors['V1ExcL4ExcL23Connection'].store_connections(data_store)    

else: 
    setup_logging()
    data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':'B'}),replace=True)
    logger.info('Loaded data store')

if mozaik.mpi_comm.rank == mozaik.MPI_ROOT:
    perform_analysis_and_visualization(data_store)
