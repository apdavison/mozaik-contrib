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

def analysis(data_store,analog_ids,analog_ids_inh,gratings,bars):
    sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4','V1_Inh_L4']))
    exc_sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4']))

    TrialAveragedFiringRate(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()

    if bars:
        TrialAveragedFiringRate(param_filter_query(data_store,st_name="FlashedBar"),ParameterSet({})).analyse()

    if gratings:
        TrialAveragedFiringRate(param_filter_query(data_store,st_name="FullfieldDriftingSinusoidalGrating"),ParameterSet({})).analyse()

    Irregularity(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()

    PSTH(param_filter_query(data_store),ParameterSet({'bin_length' : 10.0 })).analyse()

    NeuronToNeuronAnalogSignalCorrelations(param_filter_query(data_store,analysis_algorithm='PSTH',st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({'convert_nan_to_zero' : True})).analyse()
    PopulationMeanAndVar(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()

    dsv = param_filter_query(data_store,sheet_name=exc_sheets)
    ActionPotentialRemoval(dsv,ParameterSet({'window_length': 5.0})).analyse()

    dsv = param_filter_query(data_store,analysis_algorithm='ActionPotentialRemoval')
    TrialVariability(data_store,ParameterSet({'vm': True,  'cond_exc': True, 'cond_inh': True})).analyse()
    TrialMean(data_store,ParameterSet({'vm': True,  'cond_exc': True, 'cond_inh': True})).analyse()

    dsv = param_filter_query(data_store,st_name='InternalStimulus',st_direct_stimulation_name="None")
    Analog_MeanSTDAndFanoFactor(dsv,ParameterSet({})).analyse()

    if gratings:

        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name=sheets)    
        GaussianTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name=sheets)   
        Analog_F0andF1(dsv,ParameterSet({})).analyse()

        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name=sheets)  
        PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()

        TrialAveragedFiringRate(param_filter_query(data_store,sheet_name=sheets,st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({})).analyse()

        ModulationRatio(param_filter_query(data_store,sheet_name=exc_sheets,st_contrast=[100]),ParameterSet({})).analyse()

        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate')    
        CircularVarianceOfTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()

    data_store.save()



def perform_analysis_and_visualization(data_store,gratings,bars):
    
    sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4','V1_Inh_L4']))
    exc_sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4']))
    
    NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
    
    analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()
    spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()
    
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
    
    if bars:
        dsv = queries.param_filter_query(data_store,st_name='FlashedBar')
        for ads in dsv.get_analysis_result():
            sid = MozaikParametrized.idd(ads.stimulus_id)
            sid.x=0
            ads.stimulus_id = str(sid)
        for seg in dsv.get_segments():    
            sid = MozaikParametrized.idd(seg.annotations['stimulus'])
            sid.x=0
            seg.annotations['stimulus'] = str(sid)

    #analysis(data_store,analog_ids,analog_ids_inh,gratings,bars)

    # self sustained plotting
    dsv = param_filter_query(data_store,st_name='InternalStimulus')   
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalog1.png').plot()
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhAnalog1.png').plot()    
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalog2.png').plot()
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhAnalog2.png').plot()    
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalog3.png').plot()
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhAnalog3.png').plot()    
    
    RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})
    RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})

    if bars:
        queries.param_filter_query(data_store,st_name='FlashedBar',sheet_name = 'V1_Exc_L4',analysis_algorithm='TrialMean',st_relative_luminance=1).print_conent(full_ADS=True)
        dsv = queries.param_filter_query(data_store,st_name='FlashedBar',y_axis_name=['Vm (no AP) trial-to-trial mean','exc. conductance trial-to-trial mean','inh. conductance trial-to-trial mean'],sheet_name = 'V1_Exc_L4',analysis_algorithm='TrialMean',st_relative_luminance=1)
        dsv.print_conent(full_ADS=True)
        PlotTemporalTuningCurve(dsv, ParameterSet({'parameter_name' : 'y', 'neurons': list(analog_ids[:4]), 'sheet_name' : 'V1_Exc_L4','mean' : False}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='AnalogSignalEvolutionOnBar1').plot({'*.x_lim'  :(0,200),'*.y_label'  : 'offset'})
        dsv = queries.param_filter_query(data_store,st_name='FlashedBar',y_axis_name=['Vm (no AP) trial-to-trial mean','exc. conductance trial-to-trial mean','inh. conductance trial-to-trial mean'],sheet_name = 'V1_Exc_L4',analysis_algorithm='TrialMean',st_relative_luminance=0)
        PlotTemporalTuningCurve(dsv, ParameterSet({'parameter_name' : 'y', 'neurons': list(analog_ids[:4]), 'sheet_name' : 'V1_Exc_L4','mean' : False}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='AnalogSignalEvolutionOffBar1').plot({'*.x_lim'  :(0,200),'*.y_label'  : 'offset'})


        dsv = queries.param_filter_query(data_store,st_name='FlashedBar',y_axis_name=['Vm (no AP) trial-to-trial mean','exc. conductance trial-to-trial mean','inh. conductance trial-to-trial mean'],sheet_name = 'V1_Exc_L4',analysis_algorithm='TrialMean',st_relative_luminance=1)
        PlotTemporalTuningCurve(dsv, ParameterSet({'parameter_name' : 'y', 'neurons': list(analog_ids[4:8]), 'sheet_name' : 'V1_Exc_L4','mean' : False}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='AnalogSignalEvolutionOnBar2').plot({'*.x_lim'  :(0,200),'*.y_label'  : 'offset'})
        dsv = queries.param_filter_query(data_store,st_name='FlashedBar',y_axis_name=['Vm (no AP) trial-to-trial mean','exc. conductance trial-to-trial mean','inh. conductance trial-to-trial mean'],sheet_name = 'V1_Exc_L4',analysis_algorithm='TrialMean',st_relative_luminance=0)
        PlotTemporalTuningCurve(dsv, ParameterSet({'parameter_name' : 'y', 'neurons': list(analog_ids[4:8]), 'sheet_name' : 'V1_Exc_L4','mean' : False}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='AnalogSignalEvolutionOffBar2').plot({'*.x_lim'  :(0,200),'*.y_label'  : 'offset'})

        dsv = queries.param_filter_query(data_store,st_name='FlashedBar',y_axis_name=['Vm (no AP) trial-to-trial mean','exc. conductance trial-to-trial mean','inh. conductance trial-to-trial mean'],sheet_name = 'V1_Exc_L4',analysis_algorithm='TrialMean',st_relative_luminance=1)
        PlotTemporalTuningCurve(dsv, ParameterSet({'parameter_name' : 'y', 'neurons': list(analog_ids[8:12]), 'sheet_name' : 'V1_Exc_L4','mean' : False}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='AnalogSignalEvolutionOnBar3').plot({'*.x_lim'  :(0,200),'*.y_label'  : 'offset'})
        dsv = queries.param_filter_query(data_store,st_name='FlashedBar',y_axis_name=['Vm (no AP) trial-to-trial mean','exc. conductance trial-to-trial mean','inh. conductance trial-to-trial mean'],sheet_name = 'V1_Exc_L4',analysis_algorithm='TrialMean',st_relative_luminance=0)
        PlotTemporalTuningCurve(dsv, ParameterSet({'parameter_name' : 'y', 'neurons': list(analog_ids[8:12]), 'sheet_name' : 'V1_Exc_L4','mean' : False}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='AnalogSignalEvolutionOffBar3').plot({'*.x_lim'  :(0,200),'*.y_label'  : 'offset'})
    
    SpontActOverview(data_store,ParameterSet({'l4_exc_neuron' : analog_ids[0], 'l4_inh_neuron' : analog_ids_inh[0],'l23_exc_neuron' : -1,'l23_inh_neuron' : -1}),plot_file_name='SpontActOverview.png', fig_param={'dpi' : 200,'figsize': (18,14.5)}).plot()

    SpontStatisticsOverview(data_store,ParameterSet({}), fig_param={'dpi' : 200,'figsize': (18,12)},plot_file_name='SpontStatisticsOverview.png').plot()

    if gratings:    
        Kremkow_plots.OrientationTuningSummary(data_store,ParameterSet({'exc_sheet_name': 'V1_Exc_L4','inh_sheet_name': 'V1_Inh_L4'}),fig_param={'dpi' : 100,'figsize': (15,9)},plot_file_name='OrientationTuningSummaryL4.png').plot()            

        OrientationTuningSummaryAnalogSignals(data_store,ParameterSet({'exc_sheet_name1': 'V1_Exc_L4','inh_sheet_name1': 'V1_Inh_L4','exc_sheet_name2': 'None','inh_sheet_name2': 'None'}),fig_param={'dpi' : 200,'figsize': (18,12)},plot_file_name='OrientationTuningSummaryAnalogSignals.png').plot({'*.fontsize' : 19,'*.y_lim' : (0,None)})            
        
        TrialToTrialVariabilityComparison(data_store,ParameterSet({'sheet_name1' : 'V1_Exc_L4','sheet_name2' : 'V1_Exc_L4','data_dg' : 0.93 , 'data_ni' : 1.19}),fig_param={'dpi' : 200,'figsize': (15,7.5)},plot_file_name='TrialToTrialVariabilityComparison.png').plot()

        dsv = param_filter_query(data_store,st_name='NaturalImageWithEyeMovement')    
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}, 'spontaneous' : True}),plot_file_name='NMExc.png',fig_param={'dpi' : 100,'figsize': (28,12)}).plot({'Vm_plot.y_lim' : (-70,-50),'Conductance_plot.y_lim' : (0,50.0)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}, 'spontaneous' : True}),plot_file_name='NMInh.png',fig_param={'dpi' : 100,'figsize': (28,12)}).plot({'Vm_plot.y_lim' : (-70,-50),'Conductance_plot.y_lim' : (0,50.0)})
                    
        
        dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L4')    
        MeanVsVarainceOfVM(dsv,ParameterSet({'neurons' : list(analog_ids[:5])}),fig_param={'dpi' : 100,'figsize': (15,7.5)},plot_file_name='TrialToTrialMeanVsVarianceOfVM.png').plot()
        

        # tuninc curves
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate','Analog_F0andF1'])    
        PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[:6]), 'sheet_name' : 'V1_Exc_L4','centered'  : True,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='TCsExcL4.png',fig_param={'dpi' : 100,'figsize': (15,7.5)}).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})
        PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh[:6]), 'sheet_name' : 'V1_Inh_L4','centered'  : True,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='TCsInhL4.png',fig_param={'dpi' : 100,'figsize': (15,7.5)}).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})
        if l23_flag:
            PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids23[:6]), 'sheet_name' : 'V1_Exc_L2/3','centered'  : True,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='TCsExcL23.png',fig_param={'dpi' : 100,'figsize': (15,7.5)}).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})
            PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh23[:6]), 'sheet_name' : 'V1_Inh_L2/3','centered'  : True,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='TCsInhL23.png',fig_param={'dpi' : 100,'figsize': (15,7.5)}).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})
        
        MRfig(param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'None'}),plot_file_name='MR.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()
        MRfigReal(param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'None'}),plot_file_name='MR.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()

        dsv = queries.param_filter_query(data_store,value_name=['orientation HWHH of Firing rate','orientation CV(Firing rate)'],sheet_name=["V1_Exc_L2/3"],st_contrast=100)
        PerNeuronValueScatterPlot(dsv, ParameterSet({'only_matching_units' : False,'ignore_nan' : True}),plot_file_name='CVvsHWHH.png').plot({'*.x_lim' : (0,90),'*.y_lim' : (0,1.0)})
        
        TrialCrossCorrelationAnalysis(data_store,ParameterSet({'neurons1' : list(analog_ids),'sheet_name1' : 'V1_Exc_L4','neurons2' : [],'sheet_name2' : 'V1_Exc_L2/3', 'window_length' : 250}),fig_param={"dpi" : 100,"figsize": (15,6.5)},plot_file_name="trial-to-trial-cross-correlation.png").plot({'*.Vm.title' : None, '*.fontsize' : 19})


