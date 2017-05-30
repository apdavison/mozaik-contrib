import mozaik
from mozaik.visualization.plotting import *
from mozaik.analysis.technical import NeuronAnnotationsToPerNeuronValues
from mozaik.analysis.analysis import *
from mozaik.analysis.vision import *
from mozaik.storage.queries import *
from mozaik.storage.datastore import PickledDataStore
from mozaik.controller import Global

import sys
sys.path.append('/home/antolikjan/projects/mozaik/contrib')
import Kremkow_plots
from Kremkow_plots import *
from lsv1m_paper import *


logger = mozaik.getMozaikLogger()

import os

def perform_analysis_and_visualization(data_store):
    
    
    NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
    
    analog_ids = param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids()
    spike_ids = param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_spike_train_ids()

    OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSX_ONAnalog1.png').plot()
    OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : analog_ids[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSX_ONAnalog2.png').plot()
    OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : analog_ids[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSX_ONAnalog3.png').plot()

    RasterPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})
    
