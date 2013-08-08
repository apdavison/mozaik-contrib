import numpy
import mozaik
from mozaik.visualization.plotting import *
from mozaik.analysis.technical import NeuronAnnotationsToPerNeuronValues
from mozaik.analysis.analysis import *
from mozaik.storage.queries import *
from mozaik.storage.datastore import PickledDataStore


def perform_analysis_and_visualization(data_store):
    analog_ids = sorted(param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids())
    spike_ids = sorted(param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids())
    
    if True: # PLOTTING
            activity_plot_param =    {
                   'frame_rate' : 5,  
                   'bin_width' : 5.0, 
                   'scatter' :  True,
                   'resolution' : 0
            }       
            
            
            TrialAveragedFiringRate(param_filter_query(data_store,st_direct_stimulation_name="None"),ParameterSet({})).analyse()
            PopulationMean(data_store,ParameterSet({})).analyse()
            
            data_store.print_content(full_ADS=True)
            
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcAnalog1.png').plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,500.0)})

            import pylab
            pylab.show()
            
            
    
