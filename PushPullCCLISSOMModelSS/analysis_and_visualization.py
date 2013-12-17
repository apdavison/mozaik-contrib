import numpy
import mozaik
from mozaik.visualization.plotting import *
from mozaik.analysis.technical import NeuronAnnotationsToPerNeuronValues
from mozaik.analysis.analysis import *
from mozaik.analysis.vision import *
from mozaik.storage.queries import *
from mozaik.storage.datastore import PickledDataStore

class OR(Plotting):
    required_parameters = ParameterSet({
        'l4_exc_neurons': list,
        'l4_inh_neurons': list,
    })

    def subplot(self, subplotspec):
        gs = gridspec.GridSpecFromSubplotSpec(6, 18, subplot_spec=subplotspec,
                                              hspace=1.0, wspace=1.0)
        
        return {
                'Layer4Exc' : (PlotTuningCurve(self.datastore,ParameterSet({'parameter_name' : 'orientation', 'neurons': self.parameters.l4_exc_neurons, 'sheet_name' : 'V1_Exc_L4'})),gs[0:3,:],{}),
                'Layer4Inh' : (PlotTuningCurve(self.datastore,ParameterSet({'parameter_name' : 'orientation', 'neurons': self.parameters.l4_inh_neurons, 'sheet_name' : 'V1_Inh_L4'})),gs[3:6,:],{}),
        }


def perform_analysis_and_visualization(data_store):
    pref_or = numpy.pi/2
    analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()
    spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()

    NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
    l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')
    l4_exc_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Exc_L4')
    l4_exc = analog_ids[numpy.argmin([circular_dist(o,numpy.pi/2,numpy.pi)  for (o,p) in zip(l4_exc_or[0].get_value_by_id(analog_ids),l4_exc_phase[0].get_value_by_id(analog_ids))])]
    l4_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L4')
    l4_inh_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Inh_L4')
    l4_inh = analog_ids_inh[numpy.argmin([circular_dist(o,numpy.pi/2,numpy.pi)  for (o,p) in zip(l4_inh_or[0].get_value_by_id(analog_ids_inh),l4_inh_phase[0].get_value_by_id(analog_ids_inh))])]
    l4_exc_or_many = numpy.array(l4_exc_or[0].ids)[numpy.nonzero(numpy.array([circular_dist(o,numpy.pi/2,numpy.pi)  for (o,p) in zip(l4_exc_or[0].values,l4_exc_phase[0].values)]) < 0.1)[0]]
    l4_spike_exc_or_close_to_pi_half = numpy.array(spike_ids)[numpy.nonzero(numpy.array([circular_dist(o,pref_or,numpy.pi)  for o in l4_exc_or[0].get_value_by_id(spike_ids)]) < 0.1)[0]].tolist()
    l4_spike_inh_or_close_to_pi_half = numpy.array(spike_ids_inh)[numpy.nonzero(numpy.array([circular_dist(o,pref_or,numpy.pi)  for o in l4_inh_or[0].get_value_by_id(spike_ids_inh)]) < 0.1)[0]].tolist()
    
    if True:  #ANALYSIS
        dsv = param_filter_query(data_store,sheet_name='V1_Exc_L4')
        #ActionPotentialRemoval(dsv,ParameterSet({'window_length' : 10.0})).analyse()
        TrialAveragedFiringRate(param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],st_name="FullfieldDriftingSinusoidalGrating"),ParameterSet({})).analyse()
                
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name=['V1_Exc_L4','V1_Inh_L4'])    
        GaussianTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name=['V1_Exc_L4','V1_Inh_L4'])   
        Analog_F0andF1(dsv,ParameterSet({})).analyse()
        #TrialVariability(dsv,ParameterSet({'vm': True,  'cond_exc': False, 'cond_inh': False})).analyse()
        #param_filter_query(data_store,sheet_name='V1_Exc_L4',st_name='InternalStimulus').print_content()
        #GSTA(param_filter_query(data_store,sheet_name='V1_Exc_L4',st_name='InternalStimulus'),ParameterSet({'neurons' : [l4_exc], 'length' : 250.0 }),tags=['GSTA']).analyse()
        #GSTA(param_filter_query(data_store,sheet_name='V1_Exc_L4',st_orientation=[0,numpy.pi/2]),ParameterSet({'neurons' : [l4_exc], 'length' : 250.0 }),tags=['GSTA']).analyse()
        #Precision(param_filter_query(data_store,sheet_name='V1_Exc_L4'),ParameterSet({'neurons' : [l4_exc], 'bin_length' : 10.0 })).analyse()
        
        #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name=['V1_Exc_L4','V1_Inh_L4'])  
        #PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
        data_store.save()

    if False:
       ConductanceSignalListPlot(param_filter_query(data_store,analysis_algorithm='GSTA'),ParameterSet({'normalize_individually': True})).plot()

                    
    if True: # PLOTTING
            activity_plot_param =    {
                   'frame_rate' : 5,  
                   'bin_width' : 5.0, 
                   'scatter' :  True,
                   'resolution' : 0
            }       

            data_store.print_content(full_ADS=True)
            
            # self sustained plotting
            if True:
                #dsv = param_filter_query(data_store,st_name='DriftingGratingWithEyeMovement')    
                #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[11], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (17,6)},plot_file_name='NIwEM.png').plot()
                    
                dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],st_name='InternalStimulus',st_direct_stimulation_name='None')   
                RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False}),fig_param={'dpi' : 100,'figsize': (17,5)},plot_file_name='SSExcRaster.png').plot({'SpikeRasterPlot.group_trials':True})
                RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False}),fig_param={'dpi' : 100,'figsize': (17,5)},plot_file_name='SSInhRaster.png').plot({'SpikeRasterPlot.group_trials':True})

            
                
                dsv = param_filter_query(data_store,st_name='NaturalImageWithEyeMovement')    
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[4], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (17,6)},plot_file_name='DGwEM.png').plot()    
                #import pylab
                #pylab.show()



            #dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],st_name='FullfieldDriftingSinusoidalGrating',st_orientation=[0,numpy.pi/2],st_contrast=100)   
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='SSExcAnalog.png').plot()
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='SSInhAnalog.png').plot()
            
            #RasterPlot(data_store,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False}),fig_param={'dpi' : 100,'figsize': (17,5)},plot_file_name='SSExcRaster.png').plot({'SpikeRasterPlot.group_trials':True})
            #RasterPlot(data_store,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False}),fig_param={'dpi' : 100,'figsize': (17,5)},plot_file_name='SSInhRaster.png').plot({'SpikeRasterPlot.group_trials':True})

            #OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0On.png").plot()        
            #OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0Off.png").plot()
	    
            
            # orientation tuning plotting
            #dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],value_name='LGNAfferentOrientation')   
            #PerNeuronValuePlot(dsv,ParameterSet({}),plot_file_name='ORSet.png').plot()
            
            #dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],value_name='orientation preference',analysis_algorithm='PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage',st_contrast=100)    
            #PerNeuronValuePlot(dsv,ParameterSet({}),plot_file_name='ORComputed.png').plot()
            
            #dsv = param_filter_query(data_store,value_name=['orientation HWHH'],sheet_name=['V1_Exc_L4','V1_Inh_L4'])    
            #PerNeuronValueScatterPlot(dsv,ParameterSet({}),plot_file_name='HWHH.png').plot({ 'ScatterPlot.x_lim' : (0,90), 'ScatterPlot.y_lim' : (0,90), 'ScatterPlot.identity_line' : True})

            #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating')    
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="Exc.png").plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="Inh.png").plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0On.png").plot()
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0Off.png").plot()
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids())[1], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN1On.png").plot()
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_esyn_ids())[1], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN1Off.png").plot()
            
            # tuninc curves
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate'])    
            OR(dsv,ParameterSet({'l4_exc_neurons' : numpy.array(l4_spike_exc_or_close_to_pi_half)[[2,3,4,5]].tolist(),'l4_inh_neurons' : numpy.array(l4_spike_inh_or_close_to_pi_half)[[0,1,5,6]].tolist()}),fig_param={'dpi' : 100,'figsize': (16,6)},plot_file_name="OR.png").plot({'*.title' : None})
            #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(l4_spike_exc_or_close_to_pi_half), 'sheet_name' : 'V1_Exc_L4'}),plot_file_name='TCsExc.png').plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})
            #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(l4_spike_inh_or_close_to_pi_half), 'sheet_name' : 'V1_Inh_L4'}),plot_file_name='TCsInh.png').plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})
            import pylab
            pylab.show()
            

                
