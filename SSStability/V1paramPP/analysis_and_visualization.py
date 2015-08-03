import numpy
import mozaik
from mozaik.visualization.plotting import *
from mozaik.analysis.technical import NeuronAnnotationsToPerNeuronValues
from mozaik.analysis.analysis import *
from mozaik.analysis.spontaneous_activity import *
from mozaik.storage.queries import *
from mozaik.storage.datastore import PickledDataStore


def perform_analysis_and_visualization(data_store):
    analog_ids1 = sorted(param_filter_query(data_store,sheet_name="Exc1").get_segments()[0].get_stored_esyn_ids())
    analog_ids_inh1 = sorted(param_filter_query(data_store,sheet_name="Inh1").get_segments()[0].get_stored_esyn_ids())
    spike_ids1 = sorted(param_filter_query(data_store,sheet_name="Exc1").get_segments()[0].get_stored_spike_train_ids())
    spike_ids_inh1 = sorted(param_filter_query(data_store,sheet_name="Inh1").get_segments()[0].get_stored_spike_train_ids())

    analog_ids2 = sorted(param_filter_query(data_store,sheet_name="Exc2").get_segments()[0].get_stored_esyn_ids())
    analog_ids_inh2 = sorted(param_filter_query(data_store,sheet_name="Inh2").get_segments()[0].get_stored_esyn_ids())
    spike_ids2 = sorted(param_filter_query(data_store,sheet_name="Exc2").get_segments()[0].get_stored_spike_train_ids())
    spike_ids_inh2 = sorted(param_filter_query(data_store,sheet_name="Inh2").get_segments()[0].get_stored_spike_train_ids())

    
    if True: # PLOTTING
            activity_plot_param =    {
                   'frame_rate' : 5,  
                   'bin_width' : 5.0, 
                   'scatter' :  True,
                   'resolution' : 0
            }       

            dsv = param_filter_query(data_store,analysis_algorithm=['SpontaneousActivityLength'])    
            dsv.remove_ads_from_datastore()

            
            PSTH(param_filter_query(data_store,st_direct_stimulation_name="None"),ParameterSet({'bin_length' : 5.0})).analyse()
            TrialAveragedFiringRate(param_filter_query(data_store,st_direct_stimulation_name="None"),ParameterSet({})).analyse()
            Irregularity(param_filter_query(data_store,st_direct_stimulation_name="None"),ParameterSet({})).analyse()
            Analog_MeanSTDAndFanoFactor(param_filter_query(data_store,st_direct_stimulation_name="None"),ParameterSet({})).analyse()
            TrialAveragedVarianceAndVarianceRatioOfConductances(param_filter_query(data_store,st_direct_stimulation_name="None"),ParameterSet({})).analyse()
            CrossCorrelationOfExcitatoryAndInhibitoryConductances(param_filter_query(data_store,st_direct_stimulation_name="None"),ParameterSet({})).analyse()
            #GSTA(param_filter_query(data_store,st_direct_stimulation_name="None"),ParameterSet({'length' : 50, 'neurons' : analog_ids})).analyse()
            NeuronToNeuronAnalogSignalCorrelations(param_filter_query(data_store,analysis_algorithm='PSTH'),ParameterSet({'convert_nan_to_zero' : True})).analyse()
            PopulationMean(data_store,ParameterSet({})).analyse()
            SpontaneousActivityLength(data_store,ParameterSet({})).analyse()
                        
            data_store.save()

            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'Exc1', 'neuron' : analog_ids1[0], 'sheet_activity' : {}, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='Exc1Analog1.png').plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'Exc1', 'neuron' : analog_ids1[1], 'sheet_activity' : {}, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='Exc1Analog2.png').plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})    
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'Exc1', 'neuron' : analog_ids1[2], 'sheet_activity' : {}, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='Exc1Analog3.png').plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})

            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'Exc2', 'neuron' : analog_ids2[0], 'sheet_activity' : {}, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='Exc2Analog1.png').plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'Exc2', 'neuron' : analog_ids2[1], 'sheet_activity' : {}, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='Exc2Analog2.png').plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})    
            OverviewPlot(data_store,ParameterSet({'sheet_name' : 'Exc2', 'neuron' : analog_ids2[2], 'sheet_activity' : {}, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='Exc2Analog3.png').plot({'Vm_plot.y_lim' : (-80,-50),'Conductance_plot.y_lim' : (0,100.0)})

            
            RasterPlot(data_store,ParameterSet({'sheet_name' : 'Exc1', 'neurons' : spike_ids1,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (17,9)},plot_file_name='Exc1Raster.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(data_store,ParameterSet({'sheet_name' : 'Inh1', 'neurons' : spike_ids_inh1,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (17,9)},plot_file_name='Inh1Raster.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(data_store,ParameterSet({'sheet_name' : 'Exc2', 'neurons' : spike_ids2,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (17,9)},plot_file_name='Exc2Raster.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(data_store,ParameterSet({'sheet_name' : 'Inh2', 'neurons' : spike_ids_inh2,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (17,9)},plot_file_name='Inh2Raster.png').plot({'SpikeRasterPlot.group_trials':True})

            
            #dsv = param_filter_query(data_store,y_axis_name='Conductance^2')    
            #AnalogSignalPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcExcCross.png').plot()
            #AnalogSignalPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcInhCross.png').plot()
            
            dsv = param_filter_query(data_store,value_name=['Mean(ECond)','Mean(ICond)','Mean(VM)','Mean(ECond)/Mean(ICond)','Firing rate'])   
            PerNeuronValuePlot(dsv,ParameterSet({"cortical_view" : False}),plot_file_name='Histograms.png',fig_param={'dpi' : 100,'figsize': (25,12)}).plot({'*.title' : None,'HistogramPlot.plot[3,1].x_scale' : 'log','HistogramPlot.plot[3,1].log' : True,'HistogramPlot.plot[3,0].x_scale' : 'log','HistogramPlot.plot[3,0].log' : True})

	    


            
            
    
