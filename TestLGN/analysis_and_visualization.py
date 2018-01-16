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
    spike_ids_off = param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_spike_train_ids()

    OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSX_ONAnalog1.png').plot()
    OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : analog_ids[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSX_ONAnalog2.png').plot()
    OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : analog_ids[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSX_ONAnalog3.png').plot()

    RasterPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})
    
def perform_analysis_and_visualization_bar(data_store):
    
    
    NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
    
    analog_ids = param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids()
    spike_ids = param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_spike_train_ids()


    OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSX_ONAnalog1.png').plot()
    OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : analog_ids[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSX_ONAnalog2.png').plot()
    OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : analog_ids[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSX_ONAnalog3.png').plot()

    RasterPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='ONRaster.png').plot({'SpikeRasterPlot.group_trials':True})
    RasterPlot(data_store,ParameterSet({'sheet_name' : 'X_OFF', 'neurons' : lgn_off_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='ONRaster.png').plot({'SpikeRasterPlot.group_trials':True})
    
    dsv = param_filter_query(data_store,st_name=['FlashedBar'],st_relative_luminance=1)    
    RasterPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neurons' : lgn_on_ids,'trial_averaged_histogram': True, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='ONBar_ONRaster.png').plot({'SpikeRasterPlot.group_trials':True})

    dsv = param_filter_query(data_store,st_name=['FlashedBar'],st_relative_luminance=1)    
    RasterPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neurons' : lgn_off_ids,'trial_averaged_histogram': True, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='ONBar_OFFRaster.png').plot({'SpikeRasterPlot.group_trials':True})

    dsv = param_filter_query(data_store,st_name=['FlashedBar'],st_relative_luminance=0)    
    RasterPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neurons' : lgn_on_ids,'trial_averaged_histogram': True, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='OFFBar_ONRaster.png').plot({'SpikeRasterPlot.group_trials':True})

    dsv = param_filter_query(data_store,st_name=['FlashedBar'],st_relative_luminance=0)    
    RasterPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neurons' : lgn_off_ids,'trial_averaged_histogram': True, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='OFFBar_OFFRaster.png').plot({'SpikeRasterPlot.group_trials':True})

def perform_analysis_and_visualization_contrast_sensitivity(data_store):

	TrialAveragedFiringRate(param_filter_query(data_store,sheet_name=['X_ON','X_OFF'],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({})).analyse()
        TrialAveragedFiringRate(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()

        lgn_on_ids = param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_spike_train_ids()
	lgn_off_ids = param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_spike_train_ids()

        
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate'])    
	PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'contrast', 'neurons': list(lgn_on_ids), 'sheet_name' : 'X_ON','centered'  : False,'mean' : True, 'polar' : False, 'pool'  : False}),plot_file_name='ContrastTuningON.png',fig_param={'dpi' : 100,'figsize': (32,7)}).plot()
        PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'contrast', 'neurons': list(lgn_off_ids), 'sheet_name' : 'X_OFF','centered'  : False,'mean' : True, 'polar' : False, 'pool'  : False}),plot_file_name='ContrastTuningOFF.png',fig_param={'dpi' : 100,'figsize': (32,7)}).plot()



def perform_analysis_and_visualization_small(data_store):

    dsv = queries.param_filter_query(data_store,st_name='FlashedBar')
    for ads in dsv.get_analysis_result():
            sid = MozaikParametrized.idd(ads.stimulus_id)
            sid.x=0
            ads.stimulus_id = str(sid)
    
    for seg in dsv.get_segments():    
            sid = MozaikParametrized.idd(seg.annotations['stimulus'])
            sid.x=0
            seg.annotations['stimulus'] = str(sid)


    PSTH(param_filter_query(data_store,sheet_name=['X_ON','X_OFF'],st_name='FlashedBar'),ParameterSet({'bin_length' : 10})).analyse()

    NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()

    dsv = param_filter_query(data_store,sheet_name=['X_ON','X_OFF'])
    ActionPotentialRemoval(dsv,ParameterSet({'window_length': 5.0})).analyse()

    TrialVariability(data_store,ParameterSet({'vm': False,  'cond_exc': False, 'cond_inh': False})).analyse()
    TrialMean(data_store,ParameterSet({'vm': False,  'cond_exc': False, 'cond_inh': False})).analyse()

    analog_ids_on = param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids()
    analog_ids_off = param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_esyn_ids()

    lgn_on_ids = param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_spike_train_ids()
    lgn_off_ids = param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_spike_train_ids()

    OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : analog_ids_on[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='ONAnalog.png').plot()
    OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : analog_ids_off[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='OFFAnalog.png').plot()
    
    dsv = queries.param_filter_query(data_store,st_name='FlashedBar',y_axis_name=['psth (bin=10) trial-to-trial mean','Vm (no AP) trial-to-trial mean'],sheet_name = 'X_ON',analysis_algorithm=['TrialMean','PSTH'],st_relative_luminance=1)
    PlotTemporalTuningCurve(dsv, ParameterSet({'parameter_name' : 'y', 'neurons': list(lgn_on_ids), 'sheet_name' : 'X_ON','mean' : False}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='AnalogSignalEvolutionOnBar_X_ON.png').plot({'*.y_label'  : 'offset','*.interpolation' : 'nearest'})

    dsv = queries.param_filter_query(data_store,st_name='FlashedBar',y_axis_name=['psth (bin=10) trial-to-trial mean','Vm (no AP) trial-to-trial mean'],sheet_name = 'X_OFF',analysis_algorithm=['TrialMean','PSTH'],st_relative_luminance=1)
    PlotTemporalTuningCurve(dsv, ParameterSet({'parameter_name' : 'y', 'neurons': list(lgn_off_ids), 'sheet_name' : 'X_OFF','mean' : False}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='AnalogSignalEvolutionOnBar_X_OFF.png').plot({'*.y_label'  : 'offset','*.interpolation' : 'nearest'})

    dsv = queries.param_filter_query(data_store,st_name='FlashedBar',y_axis_name=['psth (bin=10) trial-to-trial mean','Vm (no AP) trial-to-trial mean'],sheet_name = 'X_ON',analysis_algorithm=['TrialMean','PSTH'],st_relative_luminance=0)
    PlotTemporalTuningCurve(dsv, ParameterSet({'parameter_name' : 'y', 'neurons': list(lgn_on_ids), 'sheet_name' : 'X_ON','mean' : False}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='AnalogSignalEvolutionOffBar_X_ON.png').plot({'*.y_label'  : 'offset','*.interpolation' : 'nearest'})

    dsv = queries.param_filter_query(data_store,st_name='FlashedBar',y_axis_name=['psth (bin=10) trial-to-trial mean','Vm (no AP) trial-to-trial mean'],sheet_name = 'X_OFF',analysis_algorithm=['TrialMean','PSTH'],st_relative_luminance=0)
    PlotTemporalTuningCurve(dsv, ParameterSet({'parameter_name' : 'y', 'neurons': list(lgn_off_ids), 'sheet_name' : 'X_OFF','mean' : False}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='AnalogSignalEvolutionOffBar_X_OFF.png').plot({'*.y_label'  : 'offset','*.interpolation' : 'nearest'})

