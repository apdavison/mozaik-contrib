#!/usr/bin/ipython -i
import sys
sys.path.append('/home/antolikjan/projects/mozaik/')
sys.path.append('/home/antolikjan/projects/topographica/')
sys.path.append('/home/antolikjan/projects/topographica/external/')
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


if True:
    params = setup_experiments('FFI',sim)    
    jens_model = PushPullCCModel(sim,params)

    experiment_list =   [
                           #Spontaneous Activity 
                           #MeasureSpontaneousActivity(jens_model,duration=147*7,num_trials=5),

                           #GRATINGS
                           #MeasureOrientationTuningFullfield(jens_model,num_orientations=2,spatial_frequency=0.8,temporal_frequency=2,grating_duration=147*7,contrasts=[5,10,20,30,40,50,60,70,80,90,100],num_trials=5),
                           MeasureOrientationTuningFullfield(jens_model,num_orientations=2,spatial_frequency=0.8,temporal_frequency=2,grating_duration=147*7,contrasts=[100],num_trials=2),
                       
                           #IMAGES WITH EYEMOVEMENT
                           #MeasureNaturalImagesWithEyeMovement(jens_model,stimulus_duration=147*7,num_trials=15),

                           #GRATINGS WITH EYEMOVEMENT
                           #MeasureDriftingSineGratingWithEyeMovement(jens_model,spatial_frequency=0.8,temporal_frequency=2,stimulus_duration=147*7,num_trials=15,contrast=100),
                       
                           #SHORT ORIENTATION TUNING
                           #MeasureOrientationTuningFullfield(jens_model,num_orientations=12,spatial_frequency=0.8,temporal_frequency=2,grating_duration=147*7,contrasts=[10,30,60,100],num_trials=5),
                           
                        ]

    data_store = run_experiments(jens_model,experiment_list)
    #jens_model.connectors['V1L4ExcL4ExcConnection'].store_connections(data_store)    
    #jens_model.connectors['V1L4ExcL4InhConnection'].store_connections(data_store)    
    #jens_model.connectors['V1L4InhL4ExcConnection'].store_connections(data_store)    
    #jens_model.connectors['V1L4InhL4InhConnection'].store_connections(data_store)    
    logger.info('Saving Datastore')
    if (not MPI) or (mpi_comm.rank == MPI_ROOT):
        data_store.save()
else:
    setup_logging()
    data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':'OR'}),replace=True)
    logger.info('Loaded data store')

import resource
print "Current memory usage: %iMB" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024))

analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()



#find neuron with preference closet to 0  
NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')
l4_exc_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Exc_L4')
l4_exc = analog_ids[numpy.argmin([circular_dist(o,numpy.pi/2,numpy.pi)  for (o,p) in zip(l4_exc_or[0].get_value_by_id(analog_ids),l4_exc_phase[0].get_value_by_id(analog_ids))])]
l4_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L4')
l4_inh_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Inh_L4')
l4_inh = analog_ids_inh[numpy.argmin([circular_dist(o,numpy.pi/2,numpy.pi)  for (o,p) in zip(l4_inh_or[0].get_value_by_id(analog_ids_inh),l4_inh_phase[0].get_value_by_id(analog_ids_inh))])]
print len(l4_exc_or[0].ids)
print len(l4_exc_or[0].values)
l4_exc_or_many = numpy.array(l4_exc_or[0].ids)[numpy.nonzero(numpy.array([circular_dist(o,numpy.pi/2,numpy.pi)  for (o,p) in zip(l4_exc_or[0].values,l4_exc_phase[0].values)]) < 0.1)[0]]

print "Prefered orientation of plotted exc neurons:"
print 'index ' + str(l4_exc)
print "Prefered phase of plotted exc neurons:"
print l4_exc_phase[0].get_value_by_id(l4_exc)
print "Prefered orientation of plotted inh neurons:"
print l4_inh_phase[0].get_value_by_id(l4_inh)
print 'index ' + str(l4_inh)
print "Prefered phase of plotted inh neurons:"
print l4_exc_phase[0].get_value_by_id(l4_exc)


if True:  #ANALYSIS
    dsv = param_filter_query(data_store,sheet_name='V1_Exc_L4')
    ActionPotentialRemoval(dsv,ParameterSet({'window_length' : 10.0})).analyse()
    TrialAveragedFiringRate(param_filter_query(data_store,sheet_name='V1_Exc_L4'),ParameterSet({'stimulus_type':"FullfieldDriftingSinusoidalGrating"})).analyse()
    TrialAveragedFiringRate(param_filter_query(data_store,sheet_name='V1_Inh_L4'),ParameterSet({'stimulus_type':"FullfieldDriftingSinusoidalGrating"})).analyse()
    
    #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate')    
    #GaussianTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name=['V1_Exc_L4','V1_Inh_L4'])   
    Analog_F0andF1(dsv,ParameterSet({})).analyse()
    TrialVariability(dsv,ParameterSet({'vm': True,  'cond_exc': False, 'cond_inh': False})).analyse()

    #GSTA(param_filter_query(data_store,sheet_name='V1_Exc_L4'),ParameterSet({'neurons' : [l4_exc], 'length' : 250.0 }),tags=['GSTA']).analyse()
    #Precision(param_filter_query(data_store,sheet_name='V1_Exc_L4'),ParameterSet({'neurons' : [l4_exc], 'bin_length' : 10.0 })).analyse()
    
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name=['V1_Exc_L4','V1_Inh_L4'])  
    PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
    data_store.save()
 
if True: # PLOTTING
    #dsv = param_filter_query(data_store,st_orientation=[0,numpy.pi/2],st_name='FullfieldDriftingSinusoidalGrating',st_trial=0)    
    #VmPlot(dsv,ParameterSet({'neuron' : l4_exc,'sheet_name' : 'V1_Exc_L4'})).plot({'Vm_plot.y_lim' : (-67,-56)})
    #dsv = param_filter_query(data_store,st_orientation=[0,numpy.pi/2],st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='ActionPotentialRemoval',st_trial=0)    
    #AnalogSignalListPlot(dsv,ParameterSet({'neurons' : [l4_exc],'sheet_name' : 'V1_Exc_L4'})).plot({'AnalogSignalPlot.y_lim' : (-67,-56)})
    #import pylab
    #pylab.show()

    F0_F1table(data_store,l4_exc)
    
    dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],value_name='LGNAfferentOrientation')   
    PerNeuronValuePlot(dsv,ParameterSet({})).plot()
    
    dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],value_name='orientation preference',
                              analysis_algorithm='PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage',st_contrast=100)    
    PerNeuronValuePlot(dsv,ParameterSet({})).plot()
    
    dsv = param_filter_query(data_store,st_orientation=[0,numpy.pi/2],st_name='FullfieldDriftingSinusoidalGrating')    
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot(
                                  {'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
    
    
    #dsv = param_filter_query(data_store,value_name=['orientation HWHH'],sheet_name=['V1_Exc_L4','V1_Inh_L4'])    
    #PerNeuronValueScatterPlot(dsv,ParameterSet({})).plot({ 'ScatterPlot.x_lim' : (0,90), 'ScatterPlot.y_lim' : (0,90), 'ScatterPlot.identity_line' : True})

    dsv = param_filter_query(data_store,st_orientation=[0,numpy.pi/2],st_name='FullfieldDriftingSinusoidalGrating')    
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
    
    # tuninc curves
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate','Analog_F0andF1'])    
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids), 'sheet_name' : 'V1_Exc_L4'})).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})

    #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialVariability'])    
    #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids), 'sheet_name' : 'V1_Exc_L4'})).plot()
    
    #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialVariability'])    
    #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh), 'sheet_name' : 'V1_Inh_L4'})).plot()

    import pylab
    pylab.show()


