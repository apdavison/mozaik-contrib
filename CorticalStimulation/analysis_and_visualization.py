import mozaik
from mozaik.visualization.plotting import *
from mozaik.analysis.technical import NeuronAnnotationsToPerNeuronValues
from mozaik.analysis.analysis import *
from mozaik.analysis.vision import *
from mozaik.storage.queries import *
from mozaik.storage.datastore import PickledDataStore
from mozaik.controller import Global

import sys
sys.path.append('/home/jantolik/projects/mozaik/contrib')
import Kremkow_plots
from Kremkow_plots import *
from lsv1m_paper import *


logger = mozaik.getMozaikLogger()

import os


    

def analysis(data_store,analog_ids,analog_ids_inh,gratings):
        sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3']))
        exc_sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4','V1_Exc_L2/3']))
        
        TrialAveragedFiringRate(param_filter_query(data_store,st_direct_stimulation_name=None,st_name='InternalStimulus'),ParameterSet({})).analyse()
        TrialAveragedFiringRate(param_filter_query(data_store,st_direct_stimulation_name="LocalStimulatorArray",st_name='InternalStimulus'),ParameterSet({})).analyse()

        Irregularity(param_filter_query(data_store,st_direct_stimulation_name=None,st_name='InternalStimulus'),ParameterSet({})).analyse()
        
        PSTH(param_filter_query(data_store),ParameterSet({'bin_length' : 10.0 })).analyse()
        
        NeuronToNeuronAnalogSignalCorrelations(param_filter_query(data_store,analysis_algorithm='PSTH',st_direct_stimulation_name=None,st_name='InternalStimulus'),ParameterSet({'convert_nan_to_zero' : True})).analyse()
        PopulationMeanAndVar(param_filter_query(data_store,st_direct_stimulation_name=None,st_name='InternalStimulus'),ParameterSet({})).analyse()

        dsv = param_filter_query(data_store,sheet_name=exc_sheets)
        ActionPotentialRemoval(dsv,ParameterSet({'window_length': 5.0})).analyse()

        dsv = param_filter_query(data_store,analysis_algorithm='ActionPotentialRemoval')
        TrialVariability(data_store,ParameterSet({'vm': True,  'cond_exc': True, 'cond_inh': True})).analyse()
        TrialMean(data_store,ParameterSet({'vm': True,  'cond_exc': True, 'cond_inh': True})).analyse()

        dsv = param_filter_query(data_store,st_name='InternalStimulus',st_direct_stimulation_name=None)
        Analog_MeanSTDAndFanoFactor(dsv,ParameterSet({})).analyse()

        data_store.save()

        if gratings:
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating')    
            NeuronToNeuronAnalogSignalCorrelations(param_filter_query(dsv,analysis_algorithm='PSTH'),ParameterSet({'convert_nan_to_zero' : True})).analyse()

            TrialAveragedFiringRate(param_filter_query(data_store,sheet_name=sheets,st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({})).analyse()
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name=sheets)    
            GaussianTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()

            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name=sheets)   
            Analog_F0andF1(dsv,ParameterSet({})).analyse()

            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name=sheets)  
            PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()

            #ModulationRatio(param_filter_query(data_store,sheet_name=exc_sheets,st_contrast=[100]),ParameterSet({})).analyse()

            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate')    
            CircularVarianceOfTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
        print "saving"
        data_store.save()


def perform_analysis_and_visualization1(data_store,gratings,cort_stim,nat_stim):


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
        
    orr = list(set([MozaikParametrized.idd(s).orientation for s in queries.param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_contrast=100).get_stimuli()]))                
        
    l4_exc_or_many_analog = numpy.array(analog_ids)[numpy.nonzero(numpy.array([circular_dist(l4_exc_or[0].get_value_by_id(i),0,numpy.pi)  for i in analog_ids]) < 0.1)[0]]
    l4_inh_or_many_analog = numpy.array(analog_ids_inh)[numpy.nonzero(numpy.array([circular_dist(l4_inh_or[0].get_value_by_id(i),0,numpy.pi)  for i in analog_ids_inh]) < 0.15)[0]]

    lgn_on_ids = param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_spike_train_ids()
    lgn_off_ids = param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_spike_train_ids()

    #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=[0,numpy.pi/2],st_contrast=100)   
    #name_prefix='GR_'
    #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L23ExcAnalog1.png').plot()
    #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L23InhAnalog1.png').plot()    
    #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L23ExcAnalog2.png').plot()
    #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L23InhAnalog2.png').plot()    
    #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L23ExcAnalog3.png').plot()
    #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L23InhAnalog3.png').plot()    

    dsv = param_filter_query(data_store,sheet_name='V1_Exc_L2/3',st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_orientation=0,st_stimulating_signal_parameters_scale=0.1)
    colors = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L2/3')[0].values/numpy.pi
    ids = data_store.get_sheet_indexes('V1_Exc_L2/3',dsv.get_segments()[0].get_stored_spike_train_ids())
    colors = colors[ids]
    mozaik.visualization.plotting.ActivityMovie(dsv,ParameterSet({'bin_width': 1.0,'scatter':  True,'resolution': 10,'sheet_name': 'V1_Exc_L2/3', 'exp_time_constant': 40}),fig_param={'dpi' : 100,'figsize': (6,6)},plot_file_name='cs_or_0_l23').plot({'*.dot_size' : 10,'*.title' : None,'*.colors' : colors})
    dsv = param_filter_query(data_store,sheet_name='V1_Exc_L2/3',st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_orientation=1.57079632679,st_stimulating_signal_parameters_scale=0.1)
    mozaik.visualization.plotting.ActivityMovie(dsv,ParameterSet({'bin_width': 1.0,'scatter':  True,'resolution': 10,'sheet_name': 'V1_Exc_L2/3', 'exp_time_constant': 40}),fig_param={'dpi' : 100,'figsize': (6,6)},plot_file_name='cs_or_90_l23').plot({'*.dot_size' : 10,'*.title' : None,'*.colors' : colors})



def perform_analysis_and_visualization(data_store,gratings,cort_stim,nat_stim):

    
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


        
    orr = list(set([MozaikParametrized.idd(s).orientation for s in queries.param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_contrast=100).get_stimuli()]))                
        
    l4_exc_or_many_analog = numpy.array(analog_ids)[numpy.nonzero(numpy.array([circular_dist(l4_exc_or[0].get_value_by_id(i),0,numpy.pi)  for i in analog_ids]) < 0.1)[0]]
    l4_inh_or_many_analog = numpy.array(analog_ids_inh)[numpy.nonzero(numpy.array([circular_dist(l4_inh_or[0].get_value_by_id(i),0,numpy.pi)  for i in analog_ids_inh]) < 0.15)[0]]

    lgn_on_ids = param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_spike_train_ids()
    lgn_off_ids = param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_spike_train_ids()

    #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name=sheets)    
    #GaussianTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()

    analysis(data_store,analog_ids,analog_ids_inh,gratings)
    
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

    if cort_stim:
        dsv = param_filter_query(data_store,st_name='InternalStimulus',st_direct_stimulation_name="LocalStimulatorArray")
        overviews(dsv,"CS_")

    if gratings:
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=[0,numpy.pi/2],st_contrast=100)   
        overviews(dsv,"GR_")

        if False:
                pref_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
                ort_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
                pref_fr_50 = 0
                ort_fr_50 = 0
                if len(param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[5],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0]).get_analysis_result()) != 0:    
                    pref_fr_50 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[5],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
                    ort_fr_50 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[5],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
                spont = param_filter_query(data_store,st_name='InternalStimulus',sheet_name=['V1_Exc_L4'],analysis_algorithm=['TrialAveragedFiringRate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)

                pylab.figure()
                pylab.bar([1,2,3,4,5],[numpy.mean(pref_fr_100),numpy.mean(ort_fr_100),numpy.mean(pref_fr_50),numpy.mean(ort_fr_50),numpy.mean(spont)])
                pylab.xticks([1.4,2.4,3.4,4.4,5.4],['PREF100','ORT100','PREF50','ORT50','SPONT'])
                pylab.xlim([0.8,6.0])
                pylab.ylabel("Firing rate")
                pylab.savefig(Global.root_directory+"Orientation_responseL4.png")
                if True:
                    pref_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L2/3'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)
                    ort_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L2/3'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)
                    pref_fr_50 = 0
                    ort_fr_50 = 0
                    
                    if(len(param_filter_query(data_store,sheet_name=['V1_Exc_L2/3'],st_contrast=[5],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0]).get_analysis_result())!=0):
                        pref_fr_50 = param_filter_query(data_store,sheet_name=['V1_Exc_L2/3'],st_contrast=[5],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)
                        ort_fr_50 = param_filter_query(data_store,sheet_name=['V1_Exc_L2/3'],st_contrast=[5],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)
                    spont = param_filter_query(data_store,st_name='InternalStimulus',sheet_name=['V1_Exc_L2/3'],analysis_algorithm=['TrialAveragedFiringRate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)

                    pylab.figure()
                    pylab.bar([1,2,3,4,5],[numpy.mean(pref_fr_100),numpy.mean(ort_fr_100),numpy.mean(pref_fr_50),numpy.mean(ort_fr_50),numpy.mean(spont)])
                    pylab.xticks([1.4,2.4,3.4,4.4,5.4],['PREF100','ORT100','PREF50','ORT50','SPONT'])
                    pylab.xlim([0.8,6.0])
                    pylab.ylabel("Firing rate")
                
                    pylab.savefig(Global.root_directory+"Orientation_responseL23.png")


        OrientationTuningSummaryFiringRates(data_store,ParameterSet({'exc_sheet_name1': 'V1_Exc_L4','inh_sheet_name1': 'V1_Inh_L4','exc_sheet_name2': 'V1_Exc_L2/3','inh_sheet_name2': 'V1_Inh_L2/3'}),fig_param={'dpi' : 200,'figsize': (18,12)},plot_file_name='OrientationTuningSummary.png').plot({'*.fontsize' : 19})            
        OrientationTuningSummaryAnalogSignals(data_store,ParameterSet({'exc_sheet_name1': 'V1_Exc_L4','inh_sheet_name1': 'V1_Inh_L4','exc_sheet_name2': 'V1_Exc_L2/3','inh_sheet_name2': 'V1_Inh_L2/3'}),fig_param={'dpi' : 200,'figsize': (18,12)},plot_file_name='OrientationTuningSummaryAnalogSignals.png').plot({'*.fontsize' : 19,'*.y_lim' : (0,None)})            

    
        SpontStatisticsOverview(data_store,ParameterSet({}), fig_param={'dpi' : 200,'figsize': (18,12)},plot_file_name='SpontStatisticsOverview.png').plot()
        SpontActOverview(data_store,ParameterSet({'l4_exc_neuron' : analog_ids[0], 'l4_inh_neuron' : analog_ids_inh[0],'l23_exc_neuron' : analog_ids23[0],'l23_inh_neuron' : analog_ids_inh23[0]}),plot_file_name='SpontActOverview.png', fig_param={'dpi' : 200,'figsize': (18,14.5)}).plot()

        dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L4')    
        MeanVsVarainceOfVM(dsv,ParameterSet({'neurons' : list(analog_ids[:5])}),fig_param={'dpi' : 100,'figsize': (15,7.5)},plot_file_name='TrialToTrialMeanVsVarianceOfVM.png').plot()
        
        dsv = queries.param_filter_query(data_store,value_name=['orientation HWHH of Firing rate','orientation CV(Firing rate)'],sheet_name=["V1_Exc_L2/3"],st_contrast=100)
        PerNeuronValueScatterPlot(dsv, ParameterSet({'only_matching_units' : False,'ignore_nan' : True}),plot_file_name='CVvsHWHH.png').plot({'*.x_lim' : (0,90),'*.y_lim' : (0,1.0)})
    
        MRfig(param_filter_query(data_store,sheet_name=['V1_Exc_L2/3','V1_Exc_L4'],st_contrast=[100],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'V1_Exc_L2/3'}),plot_file_name='MR.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()
        MRfigReal(param_filter_query(data_store,sheet_name=['V1_Exc_L2/3','V1_Exc_L4'],st_contrast=[100],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'V1_Exc_L2/3'}),plot_file_name='MRReal.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()
	
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
