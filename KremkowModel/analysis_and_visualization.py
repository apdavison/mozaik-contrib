import numpy
import mozaik
import pylab
from mozaik.visualization.plotting import *
from mozaik.analysis.technical import NeuronAnnotationsToPerNeuronValues
from mozaik.analysis.analysis import *
from mozaik.analysis.vision import *
from mozaik.storage.queries import *
from mozaik.storage.datastore import PickledDataStore
from mozaik.tools.circ_stat import circular_dist

def perform_analysis_and_visualization(data_store):
    if True:
        
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
        #dsv = param_filter_query(data_store,sheet_name='V1_Exc_L4')
        #ActionPotentialRemoval(dsv,ParameterSet({'window_length' : 10.0})).analyse()
        TrialAveragedFiringRate(param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],st_name="FullfieldDriftingSinusoidalGrating"),ParameterSet({})).analyse()
                
        #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name=['V1_Exc_L4','V1_Inh_L4'])    
        #GaussianTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name=['V1_Exc_L4','V1_Inh_L4'])   
        Analog_F0andF1(dsv,ParameterSet({})).analyse()
        #TrialVariability(dsv,ParameterSet({'vm': True,  'cond_exc': False, 'cond_inh': False})).analyse()

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
        
        #RetinalInputMovie(data_store,ParameterSet({'frame_rate' : 10})).plot()
        
        
        dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],value_name='LGNAfferentOrientation')   
        PerNeuronValuePlot(dsv,ParameterSet({}),plot_file_name="LGNAfferentOrientation.png").plot()
        
        dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],value_name='orientation preference',analysis_algorithm='PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage',st_contrast=100)    
        PerNeuronValuePlot(dsv,ParameterSet({})).plot()
        
        #dsv = param_filter_query(data_store,value_name=['orientation HWHH'],sheet_name=['V1_Exc_L4','V1_Inh_L4'])    
        #PerNeuronValueScatterPlot(dsv,ParameterSet({})).plot({ 'ScatterPlot.x_lim' : (0,90), 'ScatterPlot.y_lim' : (0,90), 'ScatterPlot.identity_line' : True})

        #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=[0,numpy.pi/2])    
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name="Exc1.png").plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating')    
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name="Exc2.png").plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="Inh.png").plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0On.png").plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0Off.png").plot()
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids())[1], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN1On.png").plot()
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_esyn_ids())[1], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN1Off.png").plot()
        
        # tuninc curves
        #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate','Analog_F0andF1'])    
        #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids), 'sheet_name' : 'V1_Exc_L4'})).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})
        #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh), 'sheet_name' : 'V1_Inh_L4'})).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})
        
        #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialVariability'])    
        #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids), 'sheet_name' : 'V1_Exc_L4'})).plot()
        
        #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialVariability'])    
        #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh), 'sheet_name' : 'V1_Inh_L4'})).plot()

        import pylab
        pylab.show()





