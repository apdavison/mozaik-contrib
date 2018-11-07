import mozaik
from mozaik.visualization.plotting import *
from mozaik.analysis.technical import NeuronAnnotationsToPerNeuronValues
from mozaik.analysis.analysis import *
from mozaik.analysis.technical import *
from mozaik.analysis.helper_functions import *
from mozaik.analysis.vision import *
from mozaik.storage.queries import *
from mozaik.storage.datastore import PickledDataStore
from mozaik.controller import Global

import sys
sys.path.append('/home/jantolik/projects/mozaik/contrib')

from  cortical_stimulation_visualization import *

logger = mozaik.getMozaikLogger()

import os


    

def analysis(data_store,analog_ids,or_ids,analog_ids_inh,gratings,cort_stim,scale):
	scales = [0.01,0.07,0.14,1.0];
	contrasts = [3,5,10,100]
        #sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3']))
        #exc_sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4','V1_Exc_L2/3']))
	sheets= 'V1_Exc_L2/3'
	exc_sheets= 'V1_Exc_L2/3'
        dsv = param_filter_query(data_store,sheet_name=exc_sheets)
        ActionPotentialRemoval(dsv,ParameterSet({'window_length': 5.0})).analyse()
        
        TrialAveragedFiringRate(param_filter_query(data_store,st_direct_stimulation_name=None,st_name='InternalStimulus'),ParameterSet({})).analyse()
        TrialAveragedFiringRate(param_filter_query(data_store,st_direct_stimulation_name="LocalStimulatorArray",st_name='InternalStimulus'),ParameterSet({})).analyse()

        Irregularity(param_filter_query(data_store,st_direct_stimulation_name=None,st_name='InternalStimulus'),ParameterSet({})).analyse()
        
        PSTH(param_filter_query(data_store),ParameterSet({'bin_length' : 10.0 })).analyse()

	if cort_stim:
	    dsv = param_filter_query(data_store,st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',analysis_algorithm='TrialAveragedFiringRate',value_name='Firing rate',sheet_name='V1_Exc_L2/3')
	    PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage(dsv,ParameterSet({'parameter_name' : 'stimulating_signal_parameters_orientation'})).analyse()

	    if scale:
		ddsv = param_filter_query(data_store,st_stimulating_signal_parameters_scale = scales)    
	    else:
		ddsv = param_filter_query(data_store,st_stimulating_signal_parameters_contrast = contrasts)    
	    
	    dsv = param_filter_query(ddsv,st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',analysis_algorithm='TrialAveragedFiringRate',value_name='Firing rate',sheet_name='V1_Exc_L2/3')    
	    GaussianTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'stimulating_signal_parameters_orientation'})).analyse()
	    #TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(ddsv,sheet_name='V1_Exc_L2/3',st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_contrast = 100,analysis_algorithm='ActionPotentialRemoval'),ParameterSet({'neurons' : list(analog_ids), 'window_min' : 400, 'window_max' : 800})).analyse()
	    #TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(ddsv,sheet_name='V1_Exc_L2/3',st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_contrast = 100,analysis_algorithm='PSTH'),ParameterSet({'neurons' : list(analog_ids), 'window_min' : 40, 'window_max' : 80})).analyse()
	    if scale:
		    dsv = param_filter_query(data_store,st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',analysis_algorithm='TrialAveragedFiringRate',value_name='Firing rate',sheet_name='V1_Exc_L2/3',st_stimulating_signal_parameters_orientation = 0,st_stimulating_signal_parameters_scale = [0.005,0.01,0.02,0.04,0.08,0.16,0.32])    
		    NakaRushtonTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'stimulating_signal_parameters_scale', 'neurons' : or_ids})).analyse()
	    else:
		    dsv = param_filter_query(data_store,st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',value_name='Firing rate',analysis_algorithm='TrialAveragedFiringRate',sheet_name='V1_Exc_L2/3',st_stimulating_signal_parameters_orientation = 0)
		    NakaRushtonTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'stimulating_signal_parameters_contrast', 'neurons' : or_ids})).analyse()



        NeuronToNeuronAnalogSignalCorrelations(param_filter_query(data_store,analysis_algorithm='PSTH',sheet_name='V1_Exc_L2/3'),ParameterSet({'convert_nan_to_zero' : True})).analyse()


        if gratings:

            #TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(data_store,sheet_name='V1_Exc_L2/3',st_name='FullfieldDriftingSinusoidalGratingA',analysis_algorithm='ActionPotentialRemoval',st_contrast=100),ParameterSet({'neurons' : list(analog_ids), 'window_min' : 400, 'window_max' : 800})).analyse()
            #TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(data_store,sheet_name='V1_Exc_L2/3',st_name='FullfieldDriftingSinusoidalGratingA',analysis_algorithm='PSTH',st_contrast=100),ParameterSet({'neurons' : list(analog_ids), 'window_min' : 40, 'window_max' : 80})).analyse()

            TrialAveragedFiringRate(param_filter_query(data_store,sheet_name=sheets,st_name='FullfieldDriftingSinusoidalGratingA'),ParameterSet({})).analyse()
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGratingA',analysis_algorithm='TrialAveragedFiringRate',value_name='Firing rate',sheet_name=sheets)    
            GaussianTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()

            #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGratingA',sheet_name=sheets)   
            #Analog_F0andF1(dsv,ParameterSet({})).analyse()

            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGratingA',analysis_algorithm='TrialAveragedFiringRate',value_name='Firing rate',sheet_name=sheets)  
            PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()

            #ModulationRatio(param_filter_query(data_store,sheet_name=exc_sheets,st_contrast=[100]),ParameterSet({})).analyse()

            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGratingA',analysis_algorithm='TrialAveragedFiringRate',value_name='Firing rate')    
            CircularVarianceOfTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()

        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGratingA',analysis_algorithm='TrialAveragedFiringRate',value_name='Firing rate',sheet_name='V1_Exc_L2/3',st_orientation = 0)
        NakaRushtonTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'contrast', 'neurons' : or_ids})).analyse()



	dsv = param_filter_query(data_store,sheet_name='V1_Exc_L2/3',st_name=['InternalStimulus','FullfieldDriftingSinusoidalGratingA'])
	mozaik.analysis.analysis.TrialMean(data_store,ParameterSet({'vm': True,  'cond_exc': True, 'cond_inh': True})).analyse()

	dsv = param_filter_query(data_store,sheet_name='V1_Exc_L2/3',analysis_algorithm=['NeuronToNeuronAnalogSignalCorrelations'])
	mozaik.analysis.analysis.TrialMean(dsv,ParameterSet({'vm': False,  'cond_exc': False, 'cond_inh': False})).analyse()

	dsv = param_filter_query(data_store,sheet_name='V1_Exc_L2/3',analysis_algorithm=['ActionPotentialRemoval'],y_axis_name= 'Vm (no AP)')
	mozaik.analysis.analysis.TrialVariability(dsv,ParameterSet({'vm': False,  'cond_exc': False, 'cond_inh': False})).analyse()

        if cort_stim:
	    if scale:
		dsv = param_filter_query(data_store,sheet_name='V1_Exc_L2/3',y_axis_name = ['Vm (no AP) trial-to-trial variance','Vm (no AP)'],st_stimulating_signal_parameters_scale = 0.14)
	    else:
		dsv = param_filter_query(data_store,sheet_name='V1_Exc_L2/3',y_axis_name = ['Vm (no AP) trial-to-trial variance','Vm (no AP)'],st_stimulating_signal_parameters_contrast = 100)
	if gratings:
	    dsv = param_filter_query(data_store,sheet_name='V1_Exc_L2/3',y_axis_name = ['Vm (no AP) trial-to-trial variance','Vm (no AP)'],st_contrast=100)

	AnalogSignal_PerNeuronMeanVar(dsv,ParameterSet({})).analyse()

	dsv = param_filter_query(data_store,sheet_name='V1_Exc_L2/3',analysis_algorithm=['AnalogSignal_PerNeuronMeanVar'],value_name= ['PerNeuronVar(Vm (no AP))','PerNeuronMean(Vm (no AP))'])
	mozaik.analysis.analysis.TrialMean(dsv,ParameterSet({'vm': False,  'cond_exc': False, 'cond_inh': False})).analyse()


	Analog_MeanSTDAndFanoFactor(data_store,ParameterSet({})).analyse()

	pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Exc_L2/3',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(VM)',st_direct_stimulation_name=None).get_analysis_result()[0]
	if cort_stim:
		dsv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Exc_L2/3',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(VM)',st_direct_stimulation_name='LocalStimulatorArray')
	else:
		dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGratingA',sheet_name='V1_Exc_L2/3',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(VM)')
	mozaik.analysis.analysis.SubtractPNVfromPNVS(pnv,dsv,ParameterSet({})).analyse()       
        
        PopulationMeanAndVar(param_filter_query(data_store,identifier='PerNeuronValue'),ParameterSet({})).analyse()

	dsv = param_filter_query(data_store,sheet_name='V1_Exc_L2/3',value_name=['contrast Mean Naka-Rushton c50 of Firing rate','contrast Mean Naka-Rushton exponent of Firing rate', 'contrast Mean Naka-Rushton scaler of Firing rate','stimulating_signal_parameters_scale Mean Naka-Rushton c50 of Firing rate', 'stimulating_signal_parameters_scale Mean Naka-Rushton scaler of Firing rate','stimulating_signal_parameters_scale Mean Naka-Rushton baseline of Firing rate','stimulating_signal_parameters_scale Mean Naka-Rushton exponent of Firing rate'])
	SummarizeSingleValues(dsv,ParameterSet({'file_name' : 'stats'})).analyse()

	data_store.save()


def perform_analysis_and_visualization(data_store,gratings,cort_stim,nat_stim,tp,scale=True,sharpness=False,lee_size=False):
    scales = [0.01,0.07,0.14,1.0];
    
    sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4','V1_Inh_L4']))
    exc_sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4']))
    
    NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
    
    analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()
    analog_ids23 = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh23 = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_esyn_ids()

    spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()
    spike_ids23 = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh23 = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_spike_train_ids()
    
    l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')
    l4_exc_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Exc_L4')
    l4_exc = analog_ids[numpy.argmin([circular_dist(o,numpy.pi/2,numpy.pi)  for (o,p) in zip(l4_exc_or[0].get_value_by_id(analog_ids),l4_exc_phase[0].get_value_by_id(analog_ids))])]
    l4_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L4')
    l4_inh_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Inh_L4')
    l4_inh = analog_ids_inh[numpy.argmin([circular_dist(o,numpy.pi/2,numpy.pi)  for (o,p) in zip(l4_inh_or[0].get_value_by_id(analog_ids_inh),l4_inh_phase[0].get_value_by_id(analog_ids_inh))])]
    l4_exc_or_many = numpy.array(l4_exc_or[0].ids)[numpy.nonzero(numpy.array([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_exc_or[0].values,l4_exc_phase[0].values)]) < 0.1)[0]]
    l4_exc_or_many = list(set(l4_exc_or_many) &  set(spike_ids))

    l23_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L2/3')[0]
    l23_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L2/3')[0]
    l23_exc_or_many = numpy.array(spike_ids23)[numpy.nonzero(numpy.array([circular_dist(l23_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids23]) < 0.1)[0]]
    l23_exc_or_many_analog = numpy.array(analog_ids23)[numpy.nonzero(numpy.array([circular_dist(l23_exc_or.get_value_by_id(i),0,numpy.pi)  for i in analog_ids23]) < 0.1)[0]]
    l23_exc_or_many_analog_inh = numpy.array(analog_ids_inh23)[numpy.nonzero(numpy.array([circular_dist(l23_inh_or.get_value_by_id(i),0,numpy.pi)  for i in analog_ids_inh23]) < 0.1)[0]]
    l23_exc_or_many = list(set(l23_exc_or_many) &  set(spike_ids23))


        
    orr = list(set([MozaikParametrized.idd(s).orientation for s in queries.param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGratingA',st_contrast=100).get_stimuli()]))                
        
    l4_exc_or_many_analog = numpy.array(analog_ids)[numpy.nonzero(numpy.array([circular_dist(l4_exc_or[0].get_value_by_id(i),0,numpy.pi)  for i in analog_ids]) < 0.1)[0]]
    l4_inh_or_many_analog = numpy.array(analog_ids_inh)[numpy.nonzero(numpy.array([circular_dist(l4_inh_or[0].get_value_by_id(i),0,numpy.pi)  for i in analog_ids_inh]) < 0.15)[0]]

    lgn_on_ids = param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_spike_train_ids()
    lgn_off_ids = param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_spike_train_ids()

    analysis(data_store,analog_ids23,l23_exc_or_many,analog_ids_inh23,gratings,cort_stim,scale)

    data_store.save()

    def overviews(dsv,name_prefix):
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L23ExcAnalog1.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L23InhAnalog1.png').plot()    
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L23ExcAnalog2.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L23InhAnalog2.png').plot()    
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L23ExcAnalog3.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L23InhAnalog3.png').plot()    

        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L4ExcAnalog1.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L4InhAnalog1.png').plot()    
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L4ExcAnalog2.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L4InhAnalog2.png').plot()    
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L4ExcAnalog3.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L4InhAnalog3.png').plot()    


        RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'SSExcRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})
        RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'SSInhRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})
    
        RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neurons' : spike_ids23,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'SSExcRasterL23.png').plot({'SpikeRasterPlot.group_trials':True})
        RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neurons' : spike_ids_inh23,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'SSInhRasterL23.png').plot({'SpikeRasterPlot.group_trials':True})

        RasterPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neurons' : lgn_on_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'XONRaster.png').plot({'SpikeRasterPlot.group_trials':True})
        RasterPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neurons' : lgn_off_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'XOFFRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})

    # spontaneous activity overview
    dsv = param_filter_query(data_store,st_name='InternalStimulus',st_direct_stimulation_name=None)   
    overviews(dsv,"SPONT_")

    x = data_store.get_neuron_postions()['V1_Exc_L2/3'][0]
    y = data_store.get_neuron_postions()['V1_Exc_L2/3'][1]
    depth = data_store.get_neuron_postions()['V1_Exc_L2/3'][2]
    ors = numpy.array(l23_exc_or.ids)[numpy.nonzero(numpy.array([circular_dist(o,0,numpy.pi) for o in l23_exc_or.values]) < 0.15)[0]]
    close_spikes = numpy.array([a for a in ors if (a in spike_ids23)])
    close_spikes_sheet_indexes = data_store.get_sheet_indexes(sheet_name='V1_Exc_L2/3',neuron_ids=close_spikes)
    close_spikes = close_spikes[numpy.sqrt(x[close_spikes_sheet_indexes]**2 + y[close_spikes_sheet_indexes]**2)<1.5]
    depth_of_close_spikes = depth[data_store.get_sheet_indexes(sheet_name='V1_Exc_L2/3',neuron_ids=close_spikes)]
    depth_of_close_spikes, close_spikes = zip(*sorted(zip(depth_of_close_spikes,close_spikes),reverse=True))
    close_analogs = list(set(close_spikes) & set(analog_ids23))

    if cort_stim:
	if sharpness:
	    OrientationTuningSharpness(data_store,ParameterSet({}),fig_param={'dpi' : 300,'figsize': (13,10)},plot_file_name='Sharpness.png').plot()
	    return;
	if lee_size:
	    OrientationTuningLEE_size(data_store,ParameterSet({}),fig_param={'dpi' : 300,'figsize': (13,10)},plot_file_name='LEESizeFigure.png').plot()
	    return;

	if scale:
		ContrastResponse(data_store,ParameterSet({'cortical_stimulation' : True}),fig_param={'dpi' : 100,'figsize': (9,3.6)},plot_file_name='ContrastResponse.png').plot()
	else:
		ContrastResponseTransformed(data_store,ParameterSet({'cortical_stimulation' : True, 'type' : tp }),fig_param={'dpi' : 100,'figsize': (9,3.6)},plot_file_name='ContrastResponseTransformed.png').plot()

	if scale:
            dsv = param_filter_query(data_store,st_name='InternalStimulus',st_direct_stimulation_name="LocalStimulatorArray")
	    LightStimulationOverview(data_store,ParameterSet({}),fig_param={'dpi' : 300,'figsize': (18,6.9)},plot_file_name='LightStimulationOverview.png').plot({'*.SpikeRasterPlot.group_trials':True,'*.SpikeRasterPlot.x_label':'neuron #','*.title' : None,'Conductance_plot.y_lim' : (0,170)})

	if scale:
        	dsv = queries.param_filter_query(data_store,st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_scale = [0.08],st_stimulating_signal_parameters_orientation = 0)
	else:
		dsv = queries.param_filter_query(data_store,st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_contrast = [100],st_stimulating_signal_parameters_orientation = 0)
	ResponseOverview(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : close_analogs[12] }),fig_param={'dpi' : 100,'figsize': (16,2.5)},plot_file_name='ResponseOverview.png').plot({ 'Conductance_plot.y_lim' : (0,150)})    



	overviews(dsv,"CS_")

	if scale:
		OrientationTuningSummaryFiringRates(data_store,ParameterSet({'cortical_stimulation' : True, 'scale' : True}),plot_file_name='OrientationTuning.png',fig_param={'dpi' : 100,'figsize': (18,3.6)}).plot()
	else:
		OrientationTuningSummaryFiringRates(data_store,ParameterSet({'cortical_stimulation' : True, 'scale' : False}),plot_file_name='OrientationTuning.png',fig_param={'dpi' : 100,'figsize': (18,3.6)}).plot()

	StatisticsOverview(data_store,ParameterSet({'cortical_stimulation' : True ,'type' : tp, 'window_length': 200}),plot_file_name='StatisticsOverview.png',fig_param={'dpi' : 100,'figsize': (12,8)}).plot()

    if gratings:
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGratingA',st_orientation=[0,numpy.pi/2],st_contrast=100)   
        overviews(dsv,"GR_")
	dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGratingA',st_orientation=[0],st_contrast=100)   
	ResponseOverview(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : close_analogs[12] }),fig_param={'dpi' : 100,'figsize': (16,2.5)},plot_file_name='ResponseOverview.png').plot({ 'Conductance_plot.y_lim' : (0,150)})
	ContrastResponseTransformed(data_store,ParameterSet({'cortical_stimulation' : False, 'type' : tp }),fig_param={'dpi' : 100,'figsize': (9,3.6)},plot_file_name='ContrastResponseTransformed.png').plot()

	OrientationTuningSummaryFiringRates(data_store,ParameterSet({'cortical_stimulation' : False, 'scale' : False}),plot_file_name='OrientationTuning.png',fig_param={'dpi' : 100,'figsize': (18,3.6)}).plot()
	StatisticsOverview(data_store,ParameterSet({'cortical_stimulation' : False ,'type' : tp, 'window_length': 50}),plot_file_name='StatisticsOverview.png',fig_param={'dpi' : 100,'figsize': (12,8)}).plot()

        if False:
                pref_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],value_name='Firing rate',st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
                ort_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],value_name='Firing rate',st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
                pref_fr_50 = 0
                ort_fr_50 = 0
                if len(param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[5],analysis_algorithm=['TrialAveragedFiringRate'],value_name='Firing rate',st_orientation=[0]).get_analysis_result()) != 0:    
                    pref_fr_50 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[5],analysis_algorithm=['TrialAveragedFiringRate'],value_name='Firing rate',st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
                    ort_fr_50 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[5],analysis_algorithm=['TrialAveragedFiringRate'],value_name='Firing rate',st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
                spont = param_filter_query(data_store,st_name='InternalStimulus',sheet_name=['V1_Exc_L4'],analysis_algorithm=['TrialAveragedFiringRate'],value_name='Firing rate',ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)

                pylab.figure()
                pylab.bar([1,2,3,4,5],[numpy.mean(pref_fr_100),numpy.mean(ort_fr_100),numpy.mean(pref_fr_50),numpy.mean(ort_fr_50),numpy.mean(spont)])
                pylab.xticks([1.4,2.4,3.4,4.4,5.4],['PREF100','ORT100','PREF50','ORT50','SPONT'])
                pylab.xlim([0.8,6.0])
                pylab.ylabel("Firing rate")
                pylab.savefig(Global.root_directory+"Orientation_responseL4.png")
                if True:
                    pref_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L2/3'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],value_name='Firing rate',st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)
                    ort_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L2/3'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],value_name='Firing rate',st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)
                    pref_fr_50 = 0
                    ort_fr_50 = 0
                    
                    if(len(param_filter_query(data_store,sheet_name=['V1_Exc_L2/3'],st_contrast=[5],analysis_algorithm=['TrialAveragedFiringRate'],value_name='Firing rate',st_orientation=[0]).get_analysis_result())!=0):
                        pref_fr_50 = param_filter_query(data_store,sheet_name=['V1_Exc_L2/3'],st_contrast=[5],analysis_algorithm=['TrialAveragedFiringRate'],value_name='Firing rate',st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)
                        ort_fr_50 = param_filter_query(data_store,sheet_name=['V1_Exc_L2/3'],st_contrast=[5],analysis_algorithm=['TrialAveragedFiringRate'],value_name='Firing rate',st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)
                    spont = param_filter_query(data_store,st_name='InternalStimulus',sheet_name=['V1_Exc_L2/3'],analysis_algorithm=['TrialAveragedFiringRate'],value_name='Firing rate',ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)

                    pylab.figure()
                    pylab.bar([1,2,3,4,5],[numpy.mean(pref_fr_100),numpy.mean(ort_fr_100),numpy.mean(pref_fr_50),numpy.mean(ort_fr_50),numpy.mean(spont)])
                    pylab.xticks([1.4,2.4,3.4,4.4,5.4],['PREF100','ORT100','PREF50','ORT50','SPONT'])
                    pylab.xlim([0.8,6.0])
                    pylab.ylabel("Firing rate")
                
                    pylab.savefig(Global.root_directory+"Orientation_responseL23.png")



    
        SpontStatisticsOverview(data_store,ParameterSet({}), fig_param={'dpi' : 200,'figsize': (12,8)},plot_file_name='SpontStatisticsOverview.png').plot()
        SpontActOverview(data_store,ParameterSet({'l4_exc_neuron' : analog_ids[0], 'l4_inh_neuron' : analog_ids_inh[0],'l23_exc_neuron' : analog_ids23[0],'l23_inh_neuron' : analog_ids_inh23[0]}),plot_file_name='SpontActOverview.png', fig_param={'dpi' : 200,'figsize': (18,14.5)}).plot()

        dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L4')    
        MeanVsVarainceOfVM(dsv,ParameterSet({'neurons' : list(analog_ids[:5])}),fig_param={'dpi' : 100,'figsize': (15,7.5)},plot_file_name='TrialToTrialMeanVsVarianceOfVM.png').plot()
        
        dsv = queries.param_filter_query(data_store,value_name=['orientation HWHH of Firing rate','orientation CV(Firing rate)'],sheet_name=["V1_Exc_L2/3"],st_contrast=100)
        PerNeuronValueScatterPlot(dsv, ParameterSet({'only_matching_units' : False,'ignore_nan' : True}),plot_file_name='CVvsHWHH.png').plot({'*.x_lim' : (0,90),'*.y_lim' : (0,1.0)})
    
        MRfig(param_filter_query(data_store,sheet_name=['V1_Exc_L2/3','V1_Exc_L4'],st_contrast=[100],st_name='FullfieldDriftingSinusoidalGratingA'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'V1_Exc_L2/3'}),plot_file_name='MR.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()
        MRfigReal(param_filter_query(data_store,sheet_name=['V1_Exc_L2/3','V1_Exc_L4'],st_contrast=[100],st_name='FullfieldDriftingSinusoidalGratingA'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'V1_Exc_L2/3'}),plot_file_name='MRReal.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()
	
	dsv = param_filter_query(data_store,st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_orientation=0,st_stimulating_signal_parameters_scale=1.0)   
	ActivityMovie(dsv,ParameterSet({'bin_width': 20.0,'scatter':  True,'resolution': 10,'sheet_name': 'V1_Exc_L2/3', 'exp_time_constant': 200}),fig_param={'dpi' : 100,'figsize': (12,6)},plot_file_name='cs_or_0_l23').plot({'*.title' : None})


    if nat_stim:        

        dsv = param_filter_query(data_store,st_name='NaturalImageWithEyeMovement')   
        overviews(dsv,"NI_")

        TrialToTrialVariabilityComparison(data_store,ParameterSet({'sheet_name1' : 'V1_Exc_L4','sheet_name2' : 'V1_Exc_L4','data_dg' : 0.93 , 'data_ni' : 1.19}),fig_param={'dpi' : 200,'figsize': (15,7.5)},plot_file_name='TrialToTrialVariabilityComparison.png').plot()

        dsv = param_filter_query(data_store,st_name='NaturalImageWithEyeMovement')    
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}, 'spontaneous' : True}),plot_file_name='NMExc.png',fig_param={'dpi' : 100,'figsize': (28,12)}).plot({'Vm_plot.y_lim' : (-70,-50),'Conductance_plot.y_lim' : (0,50.0)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}, 'spontaneous' : True}),plot_file_name='NMInh.png',fig_param={'dpi' : 100,'figsize': (28,12)}).plot({'Vm_plot.y_lim' : (-70,-50),'Conductance_plot.y_lim' : (0,50.0)})
                    
        TrialCrossCorrelationAnalysis(data_store,ParameterSet({'neurons1' : list(analog_ids),'sheet_name1' : 'V1_Exc_L4','neurons2' : list(analog_ids23),'sheet_name2' : 'V1_Exc_L2/3', 'window_length' : 250}),fig_param={"dpi" : 100,"figsize": (15,6.5)},plot_file_name="trial-to-trial-cross-correlation.png").plot({'*.Vm.title' : None, '*.fontsize' : 19})
