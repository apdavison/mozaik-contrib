#!/usr/local/bin/ipython -i 
import sys
sys.path.append('/home/jan/projects/mozaik-recording-update/')
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
                           MeasureSpontaneousActivity(jens_model,duration=50*7,num_trials=2),

                           #GRATINGS
                           #MeasureOrientationTuningFullfield(jens_model,num_orientations=2,spatial_frequency=0.8,temporal_frequency=2,grating_duration=147*7,contrasts=[30,100],num_trials=2),
                       
                           #IMAGES WITH EYEMOVEMENT
                           #MeasureNaturalImagesWithEyeMovement(jens_model,stimulus_duration=147*7,num_trials=15),

                           #GRATINGS WITH EYEMOVEMENT
                           #MeasureDriftingSineGratingWithEyeMovement(jens_model,spatial_frequency=0.8,temporal_frequency=2,stimulus_duration=147*7,num_trials=15,contrast=100),
                       
                           #SHORT ORIENTATION TUNING
                           #MeasureOrientationTuningFullfield(jens_model,num_orientations=12,spatial_frequency=0.8,temporal_frequency=2,grating_duration=147*7,contrasts=[10,30,60,100],num_trials=5),
                           
                        ]

    data_store = run_experiments(jens_model,experiment_list)
    jens_model.connectors['V1L4ExcL4ExcConnection'].store_connections(data_store)    
    jens_model.connectors['V1L4ExcL4InhConnection'].store_connections(data_store)    
    jens_model.connectors['V1L4InhL4ExcConnection'].store_connections(data_store)    
    jens_model.connectors['V1L4InhL4InhConnection'].store_connections(data_store)    
    logger.info('Saving Datastore')
    if (not MPI) or (mpi_comm.rank == MPI_ROOT):
        data_store.save()
else:
    setup_logging()
    data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':'OR'}),replace=True)
    logger.info('Loaded data store')

import resource
print "Current memory usage: %iMB" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024))


#find neuron with preference closet to 0  
NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')
l4_exc_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Exc_L4')
l4_exc = numpy.argmin([circular_dist(o,numpy.pi/2,numpy.pi)  for (o,p) in zip(l4_exc_or[0].values[:20],l4_exc_phase[0].values[:20])])
l4_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L4')
l4_inh_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Inh_L4')
l4_inh = numpy.argmin([circular_dist(o,numpy.pi/2,numpy.pi)  for (o,p) in zip(l4_inh_or[0].values[:20],l4_inh_phase[0].values[:20])])
l4_exc_or_many = numpy.nonzero(numpy.array([circular_dist(o,numpy.pi/2,numpy.pi)  for (o,p) in zip(l4_exc_or[0].values,l4_exc_phase[0].values)]) < 0.1)[0]

print "Prefered orientation of plotted exc neurons:"
print l4_exc_or[0].values[l4_exc]
print 'index ' + str(l4_exc)
print "Prefered phase of plotted exc neurons:"
print l4_exc_phase[0].values[l4_exc]
print "Prefered orientation of plotted inh neurons:"
print l4_inh_or[0].values[l4_inh]
print 'index ' + str(l4_inh)
print "Prefered phase of plotted inh neurons:"
print l4_inh_phase[0].values[l4_inh]


if False:  #ANALYSIS
    
    TrialAveragedFiringRate(param_filter_query(data_store,sheet_name='V1_Exc_L4'),ParameterSet({'stimulus_type':"FullfieldDriftingSinusoidalGrating"})).analyse()
    TrialAveragedFiringRate(param_filter_query(data_store,sheet_name='V1_Inh_L4'),ParameterSet({'stimulus_type':"FullfieldDriftingSinusoidalGrating"})).analyse()

    #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate')    
    #GaussianTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
    #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name=['V1_Exc_L4','V1_Inh_L4'])   
    #Conductance_F0andF1(dsv,ParameterSet({})).analyse()
    #TrialVariability(data_store,ParameterSet({'vm': True,  'cond_exc': False, 'cond_inh': False})).analyse()
    #TrialMean(data_store,ParameterSet({'vm': True,  'cond_exc': False, 'cond_inh': False})).analyse()
    #GSTA(param_filter_query(data_store,sheet_name='V1_Exc_L4'),ParameterSet({'neurons' : [l4_exc], 'length' : 250.0 }),tags=['GSTA']).analyse()
    #Precision(param_filter_query(data_store,sheet_name='V1_Exc_L4'),ParameterSet({'neurons' : [l4_exc], 'bin_length' : 10.0 })).analyse()
    data_store.save()
 
if False: # PLOTTING

    #F0_F1table(data_store,l4_exc)
    
    #PerNeuronValueScatterPlot(param_filter_query(data_store,value_name=['F0_Exc_Cond'],sheet_name='V1_Exc_L4',st_orientation=[0,numpy.pi/2],st_contrast=100),ParameterSet({})).plot({'ScatterPlot.title' : 'HC,DC','ScatterPlot.identity_line' : True})
    #PerNeuronValueScatterPlot(param_filter_query(data_store,value_name=['F0_Exc_Cond'],sheet_name='V1_Exc_L4',st_orientation=[0,numpy.pi/2],st_contrast=30),ParameterSet({})).plot({'ScatterPlot.title' : 'LC,DC','ScatterPlot.identity_line' : True})
    #PerNeuronValueScatterPlot(param_filter_query(data_store,value_name=['F1_Exc_Cond'],sheet_name='V1_Exc_L4',st_orientation=[0,numpy.pi/2],st_contrast=100),ParameterSet({})).plot({'ScatterPlot.title' : 'HC,F1','ScatterPlot.identity_line' : True})
    #PerNeuronValueScatterPlot(param_filter_query(data_store,value_name=['F1_Exc_Cond'],sheet_name='V1_Exc_L4',st_orientation=[0,numpy.pi/2],st_contrast=30),ParameterSet({})).plot({'ScatterPlot.title' : 'LC,F1','ScatterPlot.identity_line' : True})
    
    #dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],analysis_algorithm='TrialVariability',st_orientation = 0)    
    #PerNeuronValueScatterPlot(dsv,ParameterSet({})).plot({'ScatterPlot.identity_line' : True, 'ScatterPlot.mark_means' : True})
    
    #dsv = param_filter_query(data_store,value_name=['orientation HWHH'],sheet_name=['V1_Exc_L4','V1_Inh_L4'])    
    #PerNeuronValueScatterPlot(dsv,ParameterSet({})).plot({'ScatterPlot.title' : 'Exc', 'ScatterPlot.x_lim' : (0,90), 'ScatterPlot.y_lim' : (0,90), 'ScatterPlot.identity_line' : True})
    
    dsv = param_filter_query(data_store,st_orientation=[0,numpy.pi/2],st_name='FullfieldDriftingSinusoidalGrating')    
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
    
    
    # tuninc curves
    #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate','Conductance_F0andF1'])    
    #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(numpy.arange(0,7,1)), 'sheet_name' : 'V1_Exc_L4'})).plot()

    import pylab
    pylab.show()


