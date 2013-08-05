import numpy
import mozaik
from mozaik.visualization.plotting import *
from mozaik.analysis.technical import NeuronAnnotationsToPerNeuronValues
from mozaik.analysis.analysis import *
from mozaik.storage.queries import *
from mozaik.storage.datastore import PickledDataStore


def perform_analysis_and_visualization(data_store):
    analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()
    spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_on = param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_off = param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_spike_train_ids()
    
    
    if True: # PLOTTING
            activity_plot_param =    {
                   'frame_rate' : 5,  
                   'bin_width' : 5.0, 
                   'scatter' :  True,
                   'resolution' : 0
            }       
            
            
            #PSTH(param_filter_query(data_store,st_direct_stimulation_name="None",st_name=['Null']),ParameterSet({'bin_length' : 5.0})).analyse()
            #TrialAveragedFiringRate(param_filter_query(data_store,st_direct_stimulation_name="None",st_name=['Null']),ParameterSet({})).analyse()
            #Irregularity(param_filter_query(data_store,st_direct_stimulation_name="None",st_name=['Null']),ParameterSet({})).analyse()
            #NeuronToNeuronAnalogSignalCorrelations(param_filter_query(data_store,analysis_algorithm='PSTH'),ParameterSet({'convert_nan_to_zero' : True})).analyse()
            #PopulationMean(data_store,ParameterSet({})).analyse()
            
            #data_store.print_content(full_ADS=True)
    
    
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


            
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcAnalog.png').plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[1], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})    
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[2], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[3], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})    
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[4], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[5], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})    
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[6], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[7], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})    
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[8], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[9], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})    
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[10], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})


            
            
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='InhAnalog.png').plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[1], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[2], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[3], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[4], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})
            
            
            
            
            RasterPlot(param_filter_query(data_store,st_name=['Null']),ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False}),fig_param={'dpi' : 100,'figsize': (17,5)},plot_file_name='ExcRaster.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(param_filter_query(data_store,st_name=['Null']),ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False}),fig_param={'dpi' : 100,'figsize': (17,5)},plot_file_name='InhRaster.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(param_filter_query(data_store,st_name=['Null']),ParameterSet({'sheet_name' : 'X_ON', 'neurons' : spike_ids_on,'trial_averaged_histogram': False}),fig_param={'dpi' : 100,'figsize': (17,5)},plot_file_name='ExcRaster.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(param_filter_query(data_store,st_name=['Null']),ParameterSet({'sheet_name' : 'X_OFF', 'neurons' : spike_ids_off,'trial_averaged_histogram': False}),fig_param={'dpi' : 100,'figsize': (17,5)},plot_file_name='InhRaster.png').plot({'SpikeRasterPlot.group_trials':True})
            
            import pylab
            pylab.show()
            
            
    
