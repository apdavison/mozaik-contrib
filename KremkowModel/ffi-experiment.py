import sys
sys.path.append('/home/jan/projects/mozaik-queryunification/')
import matplotlib
import time
import pylab
from mozaik.framework.experiment import MeasureOrientationTuningFullfield, MeasureSpontaneousActivity, MeasureNaturalImagesWithEyeMovement, MeasureSizeTuning
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
                           
                           #SHORT ORIENTATION TUNING
                           MeasureOrientationTuningFullfield(jens_model,num_orientations=12,spatial_frequency=0.8,temporal_frequency=2,grating_duration=147*7,contrasts=[10,30,60,100],num_trials=5),
                           
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
    data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':'A'}))
    #SingleStimulus
    #NIWEM
    logger.info('Loaded data store')

import resource
print "Current memory usage: %iMB" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024))

if (not MPI) or (mpi_comm.rank == MPI_ROOT):
    
    activity_plot_param =    {
           'frame_rate' : 5,  
           'bin_width' : 5.0, 
           'scatter' :  True,
           'resolution' : 0
    }       
    NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
    
    # find neuron with preference closet to 0  
    if True:  
        l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')
        l4_exc_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Exc_L4')
        l4_exc = numpy.argmin([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_exc_or[0].values[:20],l4_exc_phase[0].values[:20])])
        l4_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L4')
        l4_inh_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Inh_L4')
        l4_inh = numpy.argmin([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_inh_or[0].values[:20],l4_inh_phase[0].values[:20])])
        l4_exc_or_many = numpy.nonzero(numpy.array([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_exc_or[0].values,l4_exc_phase[0].values)]) < 0.1)[0]

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
        
    #l4_exc = 13
    #l4_inh = 10

    
    
    #ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_exc,'sheet_name' : 'V1_Exc_L4','reversed' : True}),analysis_data_structure_parameter_filter_query(data_store,identifier='PerNeuronValue',value_name='LGNAfferentOrientation')).plot()
    #ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_exc,'sheet_name' : 'V1_Exc_L4','reversed' : True}),analysis_data_structure_parameter_filter_query(data_store,identifier='PerNeuronValue',value_name='LGNAfferentPhase')).plot()
    #ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_inh,'sheet_name' : 'V1_Inh_L4','reversed' : True}),analysis_data_structure_parameter_filter_query(data_store,identifier='PerNeuronValue',value_name='LGNAfferentOrientation')).plot()
    #ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_inh,'sheet_name' : 'V1_Inh_L4','reversed' : True}),analysis_data_structure_parameter_filter_query(data_store,identifier='PerNeuronValue',value_name='LGNAfferentPhase')).plot()
    
    GSTA(param_filter_query(data_store,sheet_name='V1_Exc_L4'),ParameterSet({'neurons' : [l4_exc], 'length' : 250.0 }),tags=['GSTA']).analyse()
    Precision(param_filter_query(data_store,sheet_name='V1_Exc_L4'),ParameterSet({'neurons' : [l4_exc], 'bin_length' : 10.0 })).analyse()

    l4_exc_data = param_filter_query(data_store,sheet_name='V1_Exc_L4')
    l4_inh_data = param_filter_query(data_store,sheet_name='V1_Inh_L4')

    TrialAveragedFiringRate(l4_exc_data,ParameterSet({'stimulus_type':"FullfieldDriftingSinusoidalGrating"})).analyse()
    TrialAveragedFiringRate(l4_inh_data,ParameterSet({'stimulus_type':"FullfieldDriftingSinusoidalGrating"})).analyse()
    
    OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : 0, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot()
    OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : 0, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot()
    OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot()
    #pylab.show()
    Figure2Gratings(data_store,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc}),plot_file_name='FigureL4Gratings.png',fig_param={'dpi' : 100,'figsize': (18,12)}).plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate')    
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neuron': l4_exc, 'sheet_name' : 'V1_Exc_L4'})).plot()
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neuron': l4_inh, 'sheet_name' : 'V1_Inh_L4'})).plot()
    
    PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage(param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({'parameter_name':'orientation'})).analyse()
    PerNeuronValuePlot(param_filter_query(data_store,identifier='PerNeuronValue',value_name='LGNAfferentOrientation',sheet_name='V1_Exc_L4'),ParameterSet({}),plot_file_name='ORSET.png').plot()
    PerNeuronValuePlot(param_filter_query(data_store,identifier='PerNeuronValue',value_name='orientation preference',sheet_name='V1_Exc_L4',st_name='FullfieldDriftingSinusoidalGrating',st_max_luminance=90.0),ParameterSet({}),plot_file_name='ORPREFL4.png').plot()
    
    if False:
        l4_exc_measured_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'orientation preference', sheet_name = 'V1_Exc_L4')
        l4_inh_measured_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'orientation preference', sheet_name = 'V1_Inh_L4')
        pylab.figure()
        pylab.subplot(2,1,1)
        pylab.plot(l4_exc_measured_or[0].values,l4_exc_or[0].values,'bo')
        pylab.title('Exc')
        pylab.xlabel('measured')
        pylab.ylabel('set')
        pylab.subplot(2,1,2)
        pylab.plot(l4_inh_measured_or[0].values,l4_inh_or[0].values,'bo')
        pylab.title('inh')
        pylab.xlabel('measured')
        pylab.ylabel('set')
    
    if True:
       dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate')    
       for i in xrange(0,min(10,len(l4_exc_or_many))):
           print 'cc'
           print l4_exc_or_many[i]
           print l4_exc 
           PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neuron': l4_exc_or_many[i], 'sheet_name' : 'V1_Exc_L4'})).plot()
    
    import resource
    logger.info("Current memory usage: %iMB" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024)))
    
import pylab
pylab.show()
