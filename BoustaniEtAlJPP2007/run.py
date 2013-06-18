# -*- coding: utf-8 -*-
"""
This is implementation of model of self-sustained activitity in balanced networks from: 

El Boustani, S., Pospischil, M., Rudolph-Lilith, M., & Destexhe, A. (n.d.). 
Activated cortical states: experiments, analyses and models. 
Journal of physiology, Paris, 101(1-3), 99–109.
"""
import sys
sys.path.append('/home/jan/projects/mozaik/')
import mozaik
from model import VogelsAbbott
from experiments import create_experiments
from mozaik.storage.datastore import Hdf5DataStore,PickledDataStore
from mozaik.framework.experiment_controller import run_workflow, setup_logging
from analysis_and_visualization import perform_analysis_and_visualization

if True:
    logger = mozaik.getMozaikLogger("Mozaik")
    data_store = run_workflow('FFI',VogelsAbbott,create_experiments) 
else: 
    setup_logging()
    data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':'A'}),replace=True)
    logger.info('Loaded data store')


perform_analysis_and_visualization(data_store)
