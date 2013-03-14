#!/usr/local/bin/ipython -i 
import sys
sys.path.append('/home/jan/projects/mozaik/')
import matplotlib
import time
import pylab
from mozaik.framework.experiment import *
from pyNN import nest as sim
from model import PushPullCCModel
from mozaik.framework.experiment_controller import run_experiments, setup_experiments, setup_logging
from mozaik.visualization.plotting import *
from mozaik.visualization.MRfig import *
from mozaik.analysis.analysis import *
from mozaik.analysis.technical import NeuronAnnotationsToPerNeuronValues
from mozaik.visualization.Kremkow_plots import *
from mozaik.storage.datastore import Hdf5DataStore,PickledDataStore
from NeuroTools.parameters import ParameterSet
from mozaik.storage.queries import *
from mozaik.tools.circ_stat import circular_dist
import mozaik
from contrib.misc import verify_connectivity, F0_F1table

logger = mozaik.getMozaikLogger("Mozaik")

try:
    from mpi4py import MPI
except ImportError:
    MPI = None
if MPI:
    mpi_comm = MPI.COMM_WORLD
MPI_ROOT = 0


if False:
    params = setup_experiments('FFI',sim)    
    jens_model = PushPullCCModel(sim,params)

    experiment_list =   [
                           #Size Tuning  
                           MeasureSizeTuning(jens_model,num_sizes=13,max_size=4.5,orientation=numpy.pi/2,spatial_frequency=0.8,temporal_frequency=2,grating_duration=147*7,contrasts=[100],num_trials=20),
    
                           #Spontaneous Activity 
                           #MeasureSpontaneousActivity(jens_model,duration=147*7,num_trials=5),
                    
                           #GRATINGS
                           #MeasureOrientationTuningFullfield(jens_model,num_orientations=2,spatial_frequency=0.8,temporal_frequency=2,grating_duration=147*7,contrasts=[100],num_trials=1),
                       
                           #IMAGES WITH EYEMOVEMENT
                           #MeasureNaturalImagesWithEyeMovement(jens_model,stimulus_duration=147*7,num_trials=15),

                           #GRATINGS WITH EYEMOVEMENT
                           #MeasureDriftingSineGratingWithEyeMovement(jens_model,spatial_frequency=0.8,temporal_frequency=2,stimulus_duration=147*7,num_trials=15,contrast=100),
                       
                           #SHORT ORIENTATION TUNING
                           #MeasureOrientationTuningFullfield(jens_model,num_orientations=12,spatial_frequency=0.8,temporal_frequency=2,grating_duration=147*7,contrasts=[10,30,60,100],num_trials=5),
                           
                        ]

    data_store = run_experiments(jens_model,experiment_list)
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

    logger.info('Saving Datastore')
    if (not MPI) or (mpi_comm.rank == MPI_ROOT):
        data_store.save()
else:
    setup_logging()
    data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':'D'}),replace=True)
    logger.info('Loaded data store')

import resource
print "Current memory usage: %iMB" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024))

pref_or = numpy.pi/2

#find neuron with preference closet to pref_or
l4_analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
l4_analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()
l23_analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_esyn_ids()
l23_analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_esyn_ids()
NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')
l4_exc = l4_analog_ids[numpy.argmin([circular_dist(o,pref_or,numpy.pi)  for o in l4_exc_or[0].get_value_by_id(l4_analog_ids)])]
l4_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L4')
l4_inh = l4_analog_ids_inh[numpy.argmin([circular_dist(o,pref_or,numpy.pi)  for o in l4_inh_or[0].get_value_by_id(l4_analog_ids_inh)])]
l23_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L2/3')
l23_exc = l23_analog_ids[numpy.argmin([circular_dist(o,pref_or,numpy.pi)  for o in l23_exc_or[0].get_value_by_id(l23_analog_ids)])]
l23_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L2/3')
l23_inh = l23_analog_ids_inh[numpy.argmin([circular_dist(o,pref_or,numpy.pi)  for o in l23_inh_or[0].get_value_by_id(l23_analog_ids_inh)])]
#l4_exc_or_many = numpy.array(l4_exc_or[0].ids)[numpy.nonzero(numpy.array([circular_dist(o,pref_or,numpy.pi)  for (o,p) in zip(l4_exc_or[0].values,l4_exc_phase[0].values)]) < 0.1)[0]]

# Find neurons for which spikes were recorded and whose orientation is close to pref_or
l4_spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
l4_spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()
l23_spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_spike_train_ids()
l23_spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_spike_train_ids()

l4_spike_exc_or_close_to_pi_half = numpy.array(l4_spike_ids)[numpy.nonzero(numpy.array([circular_dist(o,pref_or,numpy.pi)  for o in l4_exc_or[0].get_value_by_id(l4_spike_ids)]) < 0.1)[0]].tolist()
l4_spike_inh_or_close_to_pi_half = numpy.array(l4_spike_ids_inh)[numpy.nonzero(numpy.array([circular_dist(o,pref_or,numpy.pi)  for o in l4_inh_or[0].get_value_by_id(l4_spike_ids_inh)]) < 0.1)[0]].tolist()
l23_spike_exc_or_close_to_pi_half = numpy.array(l23_spike_ids)[numpy.nonzero(numpy.array([circular_dist(o,pref_or,numpy.pi)  for o in l23_exc_or[0].get_value_by_id(l23_spike_ids)]) < 0.1)[0]].tolist()
l23_spike_inh_or_close_to_pi_half = numpy.array(l23_spike_ids_inh)[numpy.nonzero(numpy.array([circular_dist(o,pref_or,numpy.pi)  for o in l23_inh_or[0].get_value_by_id(l23_spike_ids_inh)]) < 0.1)[0]].tolist()

print l4_spike_exc_or_close_to_pi_half


print "Prefered orientation of plotted exc neurons:"
print 'id: ' + str(l4_exc)
print "Prefered orientation of plotted inh neurons:"
print 'id: ' + str(l4_inh)


if False:
    TrialAveragedFiringRate(data_store,ParameterSet({'stimulus_type':"DriftingSinusoidalGratingDisk"})).analyse()
    #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate')    
    #GaussianTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
    #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'])   
    #Analog_F0andF1(dsv,ParameterSet({})).analyse()
    #TrialVariability(dsv,ParameterSet({'vm': True,  'cond_exc': False, 'cond_inh': False})).analyse()
    #TrialMean(dsv,ParameterSet({'vm': True,  'cond_exc': False, 'cond_inh': False})).analyse()
    
    #GSTA(param_filter_query(data_store,sheet_name='V1_Exc_L4'),ParameterSet({'neurons' : [l4_exc], 'length' : 250.0 }),tags=['GSTA']).analyse()
    #Precision(param_filter_query(data_store,sheet_name='V1_Exc_L4'),ParameterSet({'neurons' : [l4_exc], 'bin_length' : 10.0 })).analyse()
    
    #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'])  
    #PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
    data_store.save()
 
if True: # PLOTTING
    #verify_connectivity(data_store)
    
    activity_plot_param =    {
           'frame_rate' : 5,  
           'bin_width' : 20.0, 
           'scatter' :  True,
           'resolution' : 0
    }       
    
    dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'],value_name='LGNAfferentOrientation')   
    PerNeuronValuePlot(dsv,ParameterSet({})).plot()
    
    #dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'],value_name='orientation preference',analysis_algorithm='PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage',st_contrast=100)    
    #PerNeuronValuePlot(dsv,ParameterSet({})).plot()

   
    #dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],analysis_algorithm='TrialVariability',st_orientation = 0)    
    #PerNeuronValueScatterPlot(dsv,ParameterSet({})).plot({'ScatterPlot.identity_line' : True, 'ScatterPlot.mark_means' : True})
    
    #dsv = param_filter_query(data_store,value_name=['orientation HWHH'],sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'])    
    #PerNeuronValueScatterPlot(dsv,ParameterSet({})).plot({ 'ScatterPlot.x_lim' : (0,90), 'ScatterPlot.y_lim' : (0,90), 'ScatterPlot.identity_line' : True})

    dsv = param_filter_query(data_store,st_orientation=[0,numpy.pi/2],st_name='DriftingSinusoidalGratingDisk',st_contrast=100)    
    #dsv = param_filter_query(data_store,st_orientation=[0,numpy.pi/2],st_name='FullfieldDriftingSinusoidalGrating',st_contrast=100)    
    #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
    #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
    #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
    #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : l23_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
    
    
    #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
    #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
    #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : l23_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})

    
    dsv = param_filter_query(data_store,st_name='DriftingSinusoidalGratingDisk',analysis_algorithm=['TrialAveragedFiringRate'])    
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'radius', 'neurons': l4_spike_exc_or_close_to_pi_half, 'sheet_name' : 'V1_Exc_L4'})).plot({'TuningCurve_firing_rate.y_lim' : (0,50)})
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'radius', 'neurons': l4_spike_inh_or_close_to_pi_half, 'sheet_name' : 'V1_Inh_L4'})).plot({'TuningCurve_firing_rate.y_lim' : (0,50)})
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'radius', 'neurons': l23_spike_exc_or_close_to_pi_half, 'sheet_name' : 'V1_Exc_L2/3'})).plot({'TuningCurve_firing_rate.y_lim' : (0,50)})
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'radius', 'neurons': l23_spike_inh_or_close_to_pi_half, 'sheet_name' : 'V1_Inh_L2/3'})).plot({'TuningCurve_firing_rate.y_lim' : (0,50)})
    
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'radius', 'neurons': l23_spike_exc_or_close_to_pi_half[8:], 'sheet_name' : 'V1_Exc_L2/3'})).plot({'TuningCurve_firing_rate.y_lim' : (0,50)})
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'radius', 'neurons': l23_spike_exc_or_close_to_pi_half[16:], 'sheet_name' : 'V1_Exc_L2/3'})).plot({'TuningCurve_firing_rate.y_lim' : (0,50)})
    
    # tuninc curves
    #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate','Analog_F0andF1'])    
    #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(l4_analog_ids), 'sheet_name' : 'V1_Exc_L4'})).plot()
    #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate','Analog_F0andF1'])    
    #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(l23_analog_ids), 'sheet_name' : 'V1_Exc_L2/3'})).plot()

    #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialVariability'])    
    #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(l4_analog_ids), 'sheet_name' : 'V1_Exc_L4'})).plot()

    #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialVariability'])    
    #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(l4_analog_ids_inh), 'sheet_name' : 'V1_Inh_L4'})).plot()


    import pylab
    pylab.show()


