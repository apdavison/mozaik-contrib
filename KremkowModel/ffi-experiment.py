import sys
sys.path.append('/home/jan/projects/mozaik/')
import matplotlib
import time
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
                           #MeasureSpontaneousActivity(jens_model,duration=10*7,num_trials=15),
                           
                           #SHORT ORIENTATION TUNING
                           MeasureOrientationTuningFullfield(jens_model,num_orientations=6,spatial_frequency=0.8,temporal_frequency=2,grating_duration=148*7,contrasts=[1.0],num_trials=3),
                        ]

    data_store = run_experiments(jens_model,experiment_list)
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
    data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':'B'}))
    #SingleStimulus
    #NIWEM
    logger.info('Loaded data store')

import resource
print "Current memory usage: %iMB" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024))

if (not MPI) or (mpi_comm.rank == MPI_ROOT):
    
    activity_plot_param =    {
           'frame_rate' : 20,  
           'bin_width' : 50.0, 
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
        
    l4_exc = 13
    l4_inh = 10
    
    
    #ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_exc,'sheet_name' : 'V1_Exc_L4','reversed' : True}),analysis_data_structure_parameter_filter_query(data_store,identifier='PerNeuronValue',value_name='LGNAfferentOrientation')).plot()
    #ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_exc,'sheet_name' : 'V1_Exc_L4','reversed' : True}),analysis_data_structure_parameter_filter_query(data_store,identifier='PerNeuronValue',value_name='LGNAfferentPhase')).plot()
    #ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_inh,'sheet_name' : 'V1_Inh_L4','reversed' : True}),analysis_data_structure_parameter_filter_query(data_store,identifier='PerNeuronValue',value_name='LGNAfferentOrientation')).plot()
    #ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_inh,'sheet_name' : 'V1_Inh_L4','reversed' : True}),analysis_data_structure_parameter_filter_query(data_store,identifier='PerNeuronValue',value_name='LGNAfferentPhase')).plot()
    #ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_exc,'sheet_name' : 'V1_Exc_L2/3','reversed' : True})).plot()
    
    GSTA(select_result_sheet_query(data_store,"V1_Exc_L4"),ParameterSet({'neurons' : [l4_exc], 'length' : 50.0 }),tags=['GSTA']).analyse()
    Precision(select_result_sheet_query(data_store,"V1_Exc_L4"),ParameterSet({'neurons' : [l4_exc], 'bin_length' : 10.0 })).analyse()

    l4_exc_data = select_result_sheet_query(data_store,'V1_Exc_L4')
    l4_inh_data = select_result_sheet_query(data_store,'V1_Inh_L4')

    TrialAveragedFiringRate(l4_exc_data,ParameterSet({'stimulus_type':"FullfieldDriftingSinusoidalGrating"})).analyse()
    TrialAveragedFiringRate(l4_inh_data,ParameterSet({'stimulus_type':"FullfieldDriftingSinusoidalGrating"})).analyse()
   
    OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot()
    OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot()
    OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : 0, 'sheet_activity' : activity_plot_param}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot()
    OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : 0, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot()
    Figure2Gratings(data_store,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc}),plot_file_name='FigureL4Gratings.png',fig_param={'dpi' : 100,'figsize': (18,12)}).plot()
    dsv = analysis_data_structure_stimulus_filter_query(data_store,'FullfieldDriftingSinusoidalGrating')
    dsv = analysis_data_structure_parameter_filter_query(dsv,analysis_algorithm='TrialAveragedFiringRate')
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neuron': l4_exc, 'sheet_name' : 'V1_Exc_L4'})).plot()
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neuron': l4_inh, 'sheet_name' : 'V1_Inh_L4'})).plot()

    dsv = analysis_data_structure_stimulus_filter_query(data_store,'FullfieldDriftingSinusoidalGrating')
    PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage(dsv,ParameterSet({'parameter_name':'orientation'})).analyse()
    
    PerNeuronValuePlot(analysis_data_structure_parameter_filter_query(data_store,identifier='PerNeuronValue',value_name='LGNAfferentOrientation',sheet_name='V1_Exc_L4'),ParameterSet({}),plot_file_name='ORSET.png').plot()
    dsv = analysis_data_structure_stimulus_filter_query(data_store,'FullfieldDriftingSinusoidalGrating',max_luminance=90)
    PerNeuronValuePlot(analysis_data_structure_parameter_filter_query(dsv,identifier='PerNeuronValue',value_name='orientation preference',sheet_name='V1_Exc_L4'),ParameterSet({}),plot_file_name='ORPREFL4.png').plot()
    
    


    if True:
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
    
    
    import resource
    logger.info("Current memory usage: %iMB" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024)))
    
import pylab
pylab.show()
