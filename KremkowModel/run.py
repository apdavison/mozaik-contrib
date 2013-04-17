# -*- coding: utf-8 -*-
"""
This is implementation of model of push-pull connectvity: 
Jens Kremkow: Correlating Excitation and Inhibition in Visual Cortical Circuits: Functional Consequences and Biological Feasibility. PhD Thesis, 2009.
"""
import sys
sys.path.append('/home/jan/projects/mozaik0.8/')
import mozaik
from model import PushPullCCModel
from experiments import create_experiments
from mozaik.storage.datastore import Hdf5DataStore,PickledDataStore
from mozaik.framework.experiment_controller import run_workflow, setup_logging
from analysis_and_visualization import perform_analysis_and_visualization

if True:
    logger = mozaik.getMozaikLogger("Mozaik")
    data_store = run_workflow('FFI',PushPullCCModel,create_experiments)
    #jens_model.connectors['V1L4ExcL4ExcConnection'].store_connections(data_store)    
    #jens_model.connectors['V1L4ExcL4InhConnection'].store_connections(data_store)    
    #jens_model.connectors['V1L4InhL4ExcConnection'].store_connections(data_store)    
    #jens_model.connectors['V1L4InhL4InhConnection'].store_connections(data_store)    
else: 
    setup_logging()
    data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':'A'}),replace=True)
    logger.info('Loaded data store')

perform_analysis_and_visualization(data_store)
