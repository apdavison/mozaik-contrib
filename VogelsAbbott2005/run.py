# -*- coding: utf-8 -*-
"""
This is implementation of model of self-sustained activitity in balanced networks from: 
Vogels, T. P., & Abbott, L. F. (2005). 
Signal propagation and logic gating in networks of integrate-and-fire neurons. 
The Journal of neuroscience : the official journal of the Society for Neuroscience, 25(46), 10786–95. 
"""
from pyNN import nest
import sys
sys.path.append('/home/jan/projects/mozaik0.8/')
import mozaik.framework.experiment_controller
from mozaik.framework.experiment_controller import run_workflow, setup_logging
import mozaik
from experiments import create_experiments
from model import VogelsAbbott
from mozaik.storage.datastore import Hdf5DataStore,PickledDataStore
from analysis_and_visualization import perform_analysis_and_visualization

if True:
    logger = mozaik.getMozaikLogger("Mozaik")
    data_store,model = run_workflow('FFI',VogelsAbbott,create_experiments)
else: 
    setup_logging()
    data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':'A'}),replace=True)
    logger.info('Loaded data store')


perform_analysis_and_visualization(data_store)
