import numpy
import mozaik
import pylab
from mozaik.visualization.plotting import *
from mozaik.analysis.technical import NeuronAnnotationsToPerNeuronValues
from mozaik.storage.queries import *
from mozaik.storage.datastore import PickledDataStore

def perform_analysis_and_visualization(data_store):
    analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()
    spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()
    
    if True: # PLOTTING
        activity_plot_param =    {
               'frame_rate' : 5,  
               'bin_width' : 5.0, 
               'scatter' :  True,
               'resolution' : 0
        }       
        dsv = param_filter_query(data_store,st_name=['InternalStimulus'])    
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,500.0)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[1], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,500.0)})    
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[2], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,500.0)})
        
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,500.0)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[1], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,500.0)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[2], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,500.0)})
        
        
        RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False})).plot({'SpikeRasterPlot.group_trials':True})
        RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False})).plot({'SpikeRasterPlot.group_trials':True})
        
        import pylab
        pylab.show()
    

