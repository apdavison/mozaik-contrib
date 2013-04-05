#!/usr/local/bin/ipython -i 
# -*- coding: utf-8 -*-
"""
This is implementation of model of self-sustained activitity in balanced networks from: 
Vogels, T. P., & Abbott, L. F. (2005). 
Signal propagation and logic gating in networks of integrate-and-fire neurons. 
The Journal of neuroscience : the official journal of the Society for Neuroscience, 25(46), 10786–95. 
"""
import sys
sys.path.append('/home/jan/projects/mozaik/')
import matplotlib
import time
import pylab
from mozaik.framework.experiment import *
from pyNN import nest as sim
from model import VogelsAbbottGaussian
from mozaik.framework.experiment_controller import run_experiments, setup_experiments, setup_logging
from mozaik.framework.population_selector import RCRandomPercentage
from mozaik.visualization.plotting import *
from mozaik.analysis.technical import NeuronAnnotationsToPerNeuronValues
from mozaik.storage.datastore import Hdf5DataStore,PickledDataStore
from NeuroTools.parameters import ParameterSet
from mozaik.storage.queries import *
import mozaik
from mozaik.tools.circ_stat import circular_dist

logger = mozaik.getMozaikLogger("Mozaik")

if True:
    params = setup_experiments('FFI',sim)    
    model = VogelsAbbottGaussian(sim,params)
    l4exc_kick = RCRandomPercentage(model.sheets["V1_Exc_L4"],ParameterSet({'percentage': 20.0}))
    l4inh_kick = RCRandomPercentage(model.sheets["V1_Inh_L4"],ParameterSet({'percentage': 20.0}))

    experiment_list =   [
                           #Lets kick the network up into activation
                           PoissonNetworkKick(model,duration=100*7,sheet_list=["V1_Exc_L4","V1_Inh_L4"],recording_configuration_list=[l4exc_kick,l4inh_kick],lambda_list=[100,100]),
                           #Spontaneous Activity 
                           #MeasureSpontaneousActivity(model,duration=70*7,num_trials=1),
                           
                           # Evoked activity
                           #MeasureOrientationTuningFullfield(model,num_orientations=2,spatial_frequency=0.8,temporal_frequency=2,grating_duration=147*7,contrasts=[100],num_trials=1),
                           
                           #Spontaneous Activity 
                           MeasureSpontaneousActivity(model,duration=300*7,num_trials=1),
                        ]
    data_store = run_experiments(model,experiment_list)
    logger.info('Saving Datastore')
    data_store.save()
else:
    setup_logging()
    data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':'A'}),replace=True)
    logger.info('Loaded data store')

import resource
print "Current memory usage: %iMB" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024))

analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()
spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()

NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')
l4_exc_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Exc_L4')
l4_exc = analog_ids[numpy.argmin([circular_dist(o,numpy.pi/2,numpy.pi)  for (o,p) in zip(l4_exc_or[0].get_value_by_id(analog_ids),l4_exc_phase[0].get_value_by_id(analog_ids))])]
l4_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L4')
l4_inh_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Inh_L4')
l4_inh = analog_ids_inh[numpy.argmin([circular_dist(o,numpy.pi/2,numpy.pi)  for (o,p) in zip(l4_inh_or[0].get_value_by_id(analog_ids_inh),l4_inh_phase[0].get_value_by_id(analog_ids_inh))])]
l4_exc_or_many = numpy.array(l4_exc_or[0].ids)[numpy.nonzero(numpy.array([circular_dist(o,numpy.pi/2,numpy.pi)  for (o,p) in zip(l4_exc_or[0].values,l4_exc_phase[0].values)]) < 0.1)[0]]


if True: # PLOTTING
    activity_plot_param =    {
           'frame_rate' : 5,  
           'bin_width' : 5.0, 
           'scatter' :  True,
           'resolution' : 0
    }       
    dsv = param_filter_query(data_store,st_name=['InternalStimulus','Null'])    
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[1], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})    
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[2], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})
    
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[1], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[2], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})
    
    RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False})).plot({'SpikeRasterPlot.group_trials':True})
    RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False})).plot({'SpikeRasterPlot.group_trials':True})

    #dsv = param_filter_query(data_store,st_orientation=[0,numpy.pi/2],st_name='FullfieldDriftingSinusoidalGrating')    
    #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,500.0)})
    #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,500.0)})

    import pylab
    pylab.show()
    

