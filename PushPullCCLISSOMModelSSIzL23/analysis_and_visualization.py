import mozaik
from mozaik.visualization.plotting import *
from mozaik.analysis.technical import NeuronAnnotationsToPerNeuronValues
from mozaik.analysis.analysis import *
from mozaik.analysis.vision import *
from mozaik.storage.queries import *
from mozaik.storage.datastore import PickledDataStore
from mozaik.controller import Global

import sys
sys.path.append('/home/antolikjan/dev/pkg/mozaik/mozaik/contrib')
import Kremkow_plots
from Kremkow_plots import *




def analysis(data_store,analog_ids,analog_ids_inh,analog_ids23,analog_ids_inh23):
    if True:  #ANALYSIS
        TrialAveragedFiringRate(param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({})).analyse()
        TrialAveragedFiringRate(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()
        #Irregularity(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()
        PopulationMean(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()
                
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'])    
        GaussianTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'])   
        Analog_F0andF1(dsv,ParameterSet({})).analyse()
        
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'])  
        PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
        
        dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L2/3','V1_Exc_L4'])
        ActionPotentialRemoval(dsv,ParameterSet({'window_length': 5.0})).analyse()
        
        GSTA(param_filter_query(data_store,sheet_name='V1_Exc_L4'),ParameterSet({'neurons' : list(analog_ids), 'length' : 250.0 }),tags=['GSTA']).analyse()
        GSTA(param_filter_query(data_store,sheet_name='V1_Inh_L4'),ParameterSet({'neurons' : list(analog_ids_inh), 'length' : 250.0 }),tags=['GSTA']).analyse()
        GSTA(param_filter_query(data_store,sheet_name='V1_Exc_L2/3'),ParameterSet({'neurons' : list(analog_ids23), 'length' : 250.0 }),tags=['GSTA']).analyse()
        GSTA(param_filter_query(data_store,sheet_name='V1_Inh_L2/3'),ParameterSet({'neurons' : list(analog_ids_inh23), 'length' : 250.0 }),tags=['GSTA']).analyse()
        
        dsv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name=['V1_Exc_L2/3','V1_Exc_L4'])
        Analog_MeanSTDAndFanoFactor(dsv,ParameterSet({})).analyse()

        pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Exc_L4',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(ECond)',st_direct_stimulation_name='None').get_analysis_result()[0]
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Exc_L4',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Exc_Cond')
        SubtractPNVfromPNVS(pnv,dsv,ParameterSet({})).analyse()

        pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Exc_L4',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(ICond)',st_direct_stimulation_name='None').get_analysis_result()[0]
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Exc_L4',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Inh_Cond')
        SubtractPNVfromPNVS(pnv,dsv,ParameterSet({})).analyse()

        pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Exc_L4',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(VM)',st_direct_stimulation_name='None').get_analysis_result()[0]
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Exc_L4',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Vm')
        SubtractPNVfromPNVS(pnv,dsv,ParameterSet({})).analyse()
        
        PSTH(param_filter_query(data_store,sheet_name=['V1_Exc_L2/3','V1_Exc_L4']),ParameterSet({'bin_length' : 2.0 })).analyse()
        SpikeCount(param_filter_query(data_store,sheet_name=['V1_Exc_L2/3','V1_Exc_L4']),ParameterSet({'bin_length' : 13.0 })).analyse()
    
        dsv = queries.param_filter_query(data_store,y_axis_name='spike count (bin=13.0)')   
        mozaik.analysis.analysis.TrialToTrialFanoFactorOfAnalogSignal(dsv,ParameterSet({})).analyse()
        
        TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(data_store,sheet_name='V1_Exc_L4',st_name="NaturalImageWithEyeMovement",analysis_algorithm='ActionPotentialRemoval'),ParameterSet({'neurons' : list(analog_ids)})).analyse()
        TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(data_store,sheet_name='V1_Exc_L4',st_name="NaturalImageWithEyeMovement",analysis_algorithm='PSTH'),ParameterSet({'neurons' : list(analog_ids)})).analyse()
        TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(data_store,sheet_name='V1_Exc_L2/3',st_name="NaturalImageWithEyeMovement",analysis_algorithm='ActionPotentialRemoval'),ParameterSet({'neurons' : list(analog_ids23)})).analyse()
        TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(data_store,sheet_name='V1_Exc_L2/3',st_name="NaturalImageWithEyeMovement",analysis_algorithm='PSTH'),ParameterSet({'neurons' : list(analog_ids23)})).analyse()


        dsv = param_filter_query(data_store,analysis_algorithm='ActionPotentialRemoval')
        TrialVariability(data_store,ParameterSet({'vm': True,  'cond_exc': False, 'cond_inh': False})).analyse()
        TrialMean(data_store,ParameterSet({'vm': True,  'cond_exc': False, 'cond_inh': False})).analyse()
       
        data_store.save()

def perform_analysis_and_visualization_or(data_store):
    analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()
    spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()

    analog_ids23 = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh23 = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_esyn_ids()
    spike_ids23 = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh23 = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_spike_train_ids()

    
    NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
    l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')
    l4_exc_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Exc_L4')
    l4_exc = analog_ids[numpy.argmin([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_exc_or[0].get_value_by_id(analog_ids),l4_exc_phase[0].get_value_by_id(analog_ids))])]
    l4_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L4')
    l4_inh_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Inh_L4')
    l4_inh = analog_ids_inh[numpy.argmin([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_inh_or[0].get_value_by_id(analog_ids_inh),l4_inh_phase[0].get_value_by_id(analog_ids_inh))])]

    # self sustained plotting
    dsv = param_filter_query(data_store,st_name=['InternalStimulus'])    
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalog.png').plot()
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhAnalog.png').plot()    
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSL23ExcAnalog.png').plot()
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSL23InhAnalog.png').plot()    
    
    
    RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcRaster.png').plot({'SpikeRasterPlot.group_trials':True})
    RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhRaster.png').plot({'SpikeRasterPlot.group_trials':True})


    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating')
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc1.png").plot({'Vm_plot.y_lim' : (-90,-50)})
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc2.png").plot({'Vm_plot.y_lim' : (-90,-50)})
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc3.png").plot({'Vm_plot.y_lim' : (-90,-50)})
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[3], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc4.png").plot({'Vm_plot.y_lim' : (-90,-50)})
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[4], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc5.png").plot({'Vm_plot.y_lim' : (-90,-50)})
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[5], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc6.png").plot({'Vm_plot.y_lim' : (-90,-50)})

    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL231.png").plot({'Vm_plot.y_lim' : (-90,-50)})
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL232.png").plot({'Vm_plot.y_lim' : (-90,-50)})
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL233.png").plot({'Vm_plot.y_lim' : (-90,-50)})
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[3], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL234.png").plot({'Vm_plot.y_lim' : (-90,-50)})
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[4], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL235.png").plot({'Vm_plot.y_lim' : (-90,-50)})
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[5], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL236.png").plot({'Vm_plot.y_lim' : (-90,-50)})


    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh1.png").plot({'Vm_plot.y_lim' : (-90,-50)})
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh2.png").plot({'Vm_plot.y_lim' : (-90,-50)})
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="InhL231.png").plot({'Vm_plot.y_lim' : (-90,-50)})
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="InhL232.png").plot({'Vm_plot.y_lim' : (-90,-50)})


    RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedExcRaster.png').plot({'SpikeRasterPlot.group_trials':True})
    RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedInhRaster.png').plot({'SpikeRasterPlot.group_trials':True})

    RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : list(analog_ids),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedExcRasterWithAnalog.png').plot({'SpikeRasterPlot.group_trials':True})
    RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : list(analog_ids_inh),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedInhRasterWithAnalog.png').plot({'SpikeRasterPlot.group_trials':True})
    RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neurons' : list(analog_ids23),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedExcL23RasterWithAnalog.png').plot({'SpikeRasterPlot.group_trials':True})
    RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neurons' : list(analog_ids_inh23),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedInhL23RasterWithAnalog.png').plot({'SpikeRasterPlot.group_trials':True})
        
    if True:
        analysis(data_store,analog_ids,analog_ids_inh,analog_ids23,analog_ids_inh23)
    
    dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L4')    
    MeanVsVarainceOfVM(dsv,ParameterSet({'neurons' : list(analog_ids[:5])}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='TrialToTrialMeanVsVarianceOfVM.png').plot()
    #TrialCrossCorrelationAnalysis(data_store,ParameterSet({'neurons' : list(analog_ids),'window_length' : 83}),fig_param={"dpi" : 100,"figsize": (25,12)},plot_file_name="trial-to-trial-cross-correlation.png").plot()
    
    FanoFactor_Baudot_et_al(data_store,ParameterSet({}),plot_file_name='FanoFactor.png').plot()
    TrialToTrialVariabilityComparison(data_store,ParameterSet({}),plot_file_name='TtTVar.png').plot()
    
    pref_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(analog_ids)
    ort_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(analog_ids)
    pref_fr_50 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[50],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(analog_ids)
    ort_fr_50 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[50],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(analog_ids)

    pylab.figure()
    pylab.bar([1,2,3,4],[numpy.mean(pref_fr_100),numpy.mean(ort_fr_100),numpy.mean(pref_fr_50),numpy.mean(ort_fr_50)])
    pylab.xticks([1.4,2.4,3.4,4.4],['PREF100','ORT100','PREF50','ORT50'])
    pylab.xlim([0.8,5.0])
    pylab.ylabel("Firing rate")
    
    pylab.savefig(Global.root_directory+"Orientation_response.png")
    
    
    #SNRAnalysis(data_store,ParameterSet({"neuron" : l4_exc}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='SNR1.png').plot()                        
    
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : l4_exc,'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[0],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview1.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[1],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview2.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[2],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview3.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[3],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview4.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[4],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview5.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[5],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview6.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[6],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview7.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[7],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview8.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[8],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview9.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[9],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview10.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[10],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview11.png').plot()

    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_exc, 'reversed' : True,'sheet_name' : 'V1_Exc_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='ExcConnections1.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : analog_ids[0], 'reversed' : True,'sheet_name' : 'V1_Exc_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='ExcConnections2.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : analog_ids[1], 'reversed' : True,'sheet_name' : 'V1_Exc_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='ExcConnections3.png').plot()

    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_inh, 'reversed' : True,'sheet_name' : 'V1_Inh_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='InhConnections1.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : analog_ids_inh[0], 'reversed' : True,'sheet_name' : 'V1_Inh_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='InhConnections2.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : analog_ids_inh[1], 'reversed' : True,'sheet_name' : 'V1_Inh_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='InhConnections3.png').plot()

    
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : analog_ids23[0], 'reversed' : True,'sheet_name' : 'V1_Exc_L2/3'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='L23ExcConnections1.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : analog_ids23[1], 'reversed' : True,'sheet_name' : 'V1_Exc_L2/3'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='L23ExcConnections2.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : analog_ids23[2], 'reversed' : True,'sheet_name' : 'V1_Exc_L2/3'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='L23ExcConnections3.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : analog_ids23[3], 'reversed' : True,'sheet_name' : 'V1_Exc_L2/3'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='L23ExcConnections4.png').plot()



def perform_analysis_and_visualization(data_store):
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
    
    orr = list(set([MozaikParametrized.idd(s).orientation for s in queries.param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_contrast=100).get_stimuli()]))                
    
    if True:
       analysis(data_store,analog_ids,analog_ids_inh,analog_ids23,analog_ids_inh23)
    
    
    #param_filter_query(data_store,analysis_algorithm='TrialToTrialCrossCorrelationOfAnalogSignalList').print_content(full_ADS=True)
    #TrialCrossCorrelationAnalysis(data_store,ParameterSet({'neurons' : list(analog_ids),'window_length' : 83}),fig_param={"dpi" : 100,"figsize": (25,12)},plot_file_name="trial-to-trial-cross-correlation.png").plot()
       
    
    if True: # PLOTTING
            activity_plot_param =    {
                   'frame_rate' : 5,  
                   'bin_width' : 5.0, 
                   'scatter' :  True,
                   'resolution' : 0
            }   
            
            MeanVsVarainceOfVM(data_store,ParameterSet({'neurons' : list(analog_ids[:5])}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='TrialToTrialMeanVsVarianceOfVM.png').plot()
             
            Kremkow_plots.OrientationTuningSummary(data_store,ParameterSet({}),fig_param={'dpi' : 100,'figsize': (20,12)},plot_file_name='OrientationTuningSummary.png').plot()            
            
            Kremkow_plots.ConductanceAndVmTuningSummary(data_store,ParameterSet({'many' : True}),fig_param={'dpi' : 100,'figsize': (25,16)},plot_file_name='Cond&VMTuning.png').plot()

            # self sustained plotting
            dsv = param_filter_query(data_store,st_name=['InternalStimulus'])    
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalog.png').plot()
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhAnalog.png').plot()    
            
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[2], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalog2.png').plot()
            
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcRaster.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhRaster.png').plot({'SpikeRasterPlot.group_trials':True})
           
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=numpy.pi/2,st_contrast=100)   
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : l4_exc,'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview.png').plot()
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[0]),numpy.pi)  for o in orr])],st_contrast=100)   
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[0],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview1.png').plot()
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[1]),numpy.pi)  for o in orr])],st_contrast=100)   
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[1],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview2.png').plot()
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[2]),numpy.pi)  for o in orr])],st_contrast=100)   
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[2],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview3.png').plot()
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[3]),numpy.pi)  for o in orr])],st_contrast=100)   
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[3],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview4.png').plot()
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[4]),numpy.pi)  for o in orr])],st_contrast=100)   
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[4],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview5.png').plot()
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[5]),numpy.pi)  for o in orr])],st_contrast=100)   
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[5],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview6.png').plot()
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[6]),numpy.pi)  for o in orr])],st_contrast=100)   
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[6],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview7.png').plot()
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[7]),numpy.pi)  for o in orr])],st_contrast=100)   
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[7],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview8.png').plot()
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[9]),numpy.pi)  for o in orr])],st_contrast=100)   
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[8],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview9.png').plot()
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[9]),numpy.pi)  for o in orr])],st_contrast=100)   
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[9],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview10.png').plot()
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[10]),numpy.pi)  for o in orr])],st_contrast=100)   
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[10],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverview11.png').plot()
            
            #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=numpy.pi/2,st_contrast=100)   
            #KremkowOverviewFigure(dsv,ParameterSet({'neuron' : l4_exc,'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview.png').plot()

            if True:            
                MeanVsVarainceOfVM(data_store,ParameterSet({'neurons' : list(analog_ids[:5])}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='TrialToTrialMeanVsVarianceOfVM.png').plot()
            
                dsv = param_filter_query(data_store,st_name='NaturalImageWithEyeMovement')            
                KremkowOverviewFigure(dsv,ParameterSet({'neuron' : l4_exc,'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='NMOverview.png').plot()
                
                FanoFactor_Baudot_et_al(data_store,ParameterSet({}),plot_file_name='FanoFactor.png').plot()
                TrialToTrialVariabilityComparison(data_store,ParameterSet({}),plot_file_name='TtTVar.png').plot()
                
                SNRAnalysis(data_store,ParameterSet({"neuron" : analog_ids[0]}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='SNR1.png').plot()                        
                SNRAnalysis(data_store,ParameterSet({"neuron" : analog_ids[1]}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='SNR2.png').plot()                        
                SNRAnalysis(data_store,ParameterSet({"neuron" : analog_ids[2]}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='SNR3.png').plot()                        

                StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids[0],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparison1.png').plot()
                StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids[1],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparison2.png').plot()
                StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids[2],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparison3.png').plot()
                StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids[3],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparison4.png').plot()
                
                dsv = param_filter_query(data_store,identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation')
                ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_exc, 'reversed' : True,'sheet_name' : 'V1_Exc_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='ExcConnections1.png').plot()
                ConnectivityPlot(data_store,ParameterSet({'neuron' : analog_ids[0], 'reversed' : True,'sheet_name' : 'V1_Exc_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='ExcConnections2.png').plot()
                ConnectivityPlot(data_store,ParameterSet({'neuron' : analog_ids[1], 'reversed' : True,'sheet_name' : 'V1_Exc_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='ExcConnections3.png').plot()

                ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_inh, 'reversed' : True,'sheet_name' : 'V1_Inh_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='InhConnections1.png').plot()
                ConnectivityPlot(data_store,ParameterSet({'neuron' : analog_ids_inh[0], 'reversed' : True,'sheet_name' : 'V1_Inh_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='InhConnections2.png').plot()
                ConnectivityPlot(data_store,ParameterSet({'neuron' : analog_ids_inh[1], 'reversed' : True,'sheet_name' : 'V1_Inh_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='InhConnections3.png').plot()

            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating')    
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedExcRaster.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedInhRaster.png').plot({'SpikeRasterPlot.group_trials':True})

            #OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0On.png").plot()        
            #OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0Off.png").plot()
    
            
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating')    
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            
            
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc1.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc2.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc3.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[3], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc4.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[4], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc5.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0On.png").plot()
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0Off.png").plot()
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids())[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN1On.png").plot()
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_esyn_ids())[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN1Off.png").plot()
            
            # tuninc curves
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate','Analog_F0andF1'])    
            PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids)[6], 'sheet_name' : 'V1_Exc_L4','centered'  : True,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='TCsExc.png',fig_param={'dpi' : 100,'figsize': (32,15)}).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})
            PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh)[6], 'sheet_name' : 'V1_Inh_L4','centered'  : True,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='TCsInh.png',fig_param={'dpi' : 100,'figsize': (32,15)}).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})

            dsv = param_filter_query(data_store,st_name='NaturalImageWithEyeMovement')    
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}, 'spontaneous' : True}),plot_file_name='NMExc.png',fig_param={'dpi' : 100,'figsize': (28,12)}).plot({'Vm_plot.y_lim' : (-70,-50),'Conductance_plot.y_lim' : (0,50.0)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}, 'spontaneous' : True}),plot_file_name='NMInh.png',fig_param={'dpi' : 100,'figsize': (28,12)}).plot({'Vm_plot.y_lim' : (-70,-50),'Conductance_plot.y_lim' : (0,50.0)})

            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_inh_or[0].get_value_by_id(l4_inh),numpy.pi)  for o in orr])],st_contrast=100)   
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : l4_inh,'sheet_name' : 'V1_Inh_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverviewInh.png').plot()
    
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_inh_or[0].get_value_by_id(analog_ids_inh[0]),numpy.pi)  for o in orr])],st_contrast=100)   
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids_inh[0],'sheet_name' : 'V1_Inh_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='InhOverview1.png').plot()
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_inh_or[0].get_value_by_id(analog_ids_inh[1]),numpy.pi)  for o in orr])],st_contrast=100)   
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids_inh[1],'sheet_name' : 'V1_Inh_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='InhOverview2.png').plot()
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_inh_or[0].get_value_by_id(analog_ids_inh[2]),numpy.pi)  for o in orr])],st_contrast=100)   
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids_inh[2],'sheet_name' : 'V1_Inh_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='InhOverview3.png').plot()
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_inh_or[0].get_value_by_id(analog_ids_inh[3]),numpy.pi)  for o in orr])],st_contrast=100)   
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids_inh[3],'sheet_name' : 'V1_Inh_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='InhOverview4.png').plot()
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_inh_or[0].get_value_by_id(analog_ids_inh[4]),numpy.pi)  for o in orr])],st_contrast=100)   
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids_inh[4],'sheet_name' : 'V1_Inh_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='InhOverview5.png').plot()

            dsv = param_filter_query(data_store,st_name=['InternalStimulus'])    
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {},'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhAnalog.png').plot()
    
            # orientation tuning plotting
            #dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],value_name='LGNAfferentOrientation')   
            #PerNeuronValuePlot(dsv,ParameterSet({"cortical_view" : True}),plot_file_name='ORSet.png').plot()
            
            #dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],value_name='orientation preference of Firing rate',analysis_algorithm='PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage',st_contrast=100)    
            #PerNeuronValuePlot(dsv,ParameterSet({"cortical_view" : True}),plot_file_name='ORComputed.png').plot()
            
