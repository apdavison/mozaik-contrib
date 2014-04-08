import mozaik
from mozaik.visualization.plotting import *
from mozaik.analysis.technical import NeuronAnnotationsToPerNeuronValues
from mozaik.analysis.analysis import *
from mozaik.analysis.vision import *
from mozaik.storage.queries import *
from mozaik.storage.datastore import PickledDataStore

import sys
sys.path.append('/home/antolikjan/dev/pkg/mozaik/mozaik/contrib/')
import Kremkow_plots
from Kremkow_plots import *

def perform_analysis_and_visualization_or(data_store):
    analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()
    spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()
    
    NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
    l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')
    l4_exc_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Exc_L4')
    l4_exc = analog_ids[numpy.argmin([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_exc_or[0].get_value_by_id(analog_ids),l4_exc_phase[0].get_value_by_id(analog_ids))])]
    l4_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L4')
    l4_inh_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Inh_L4')
    l4_inh = analog_ids_inh[numpy.argmin([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_inh_or[0].get_value_by_id(analog_ids_inh),l4_inh_phase[0].get_value_by_id(analog_ids_inh))])]
    
    GSTA(param_filter_query(data_store,sheet_name='V1_Exc_L4'),ParameterSet({'neurons' : list(analog_ids), 'length' : 250.0 }),tags=['GSTA']).analyse()

    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : l4_exc,'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[0],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview1.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[1],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview2.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[2],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview3.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[3],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview4.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[4],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview5.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[5],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview6.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[6],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview7.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[7],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview8.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[8],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview9.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[9],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview10.png').plot()
    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[10],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview11.png').plot()



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
    
    if True:  #ANALYSIS
        TrialAveragedFiringRate(param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({})).analyse()
        
        TrialAveragedFiringRate(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()
        #Irregularity(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()
        PopulationMean(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()
                
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name=['V1_Exc_L4','V1_Inh_L4'])    
        GaussianTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name=['V1_Exc_L4','V1_Inh_L4'])   
        Analog_F0andF1(dsv,ParameterSet({})).analyse()
        #TrialVariability(dsv,ParameterSet({'vm': True,  'cond_exc': False, 'cond_inh': False})).analyse()

        GSTA(param_filter_query(data_store,sheet_name='V1_Exc_L4'),ParameterSet({'neurons' : list(analog_ids), 'length' : 250.0 }),tags=['GSTA']).analyse()
        #Precision(param_filter_query(data_store,sheet_name='V1_Exc_L4'),ParameterSet({'neurons' : [l4_exc], 'bin_length' : 10.0 })).analyse()
        
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name=['V1_Exc_L4','V1_Inh_L4'])  
        PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
        data_store.save()
            
    if True: # PLOTTING
            activity_plot_param =    {
                   'frame_rate' : 5,  
                   'bin_width' : 5.0, 
                   'scatter' :  True,
                   'resolution' : 0
            }    
            
            # self sustained plotting
            dsv = param_filter_query(data_store,st_name=['InternalStimulus'])    
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalog.png').plot()
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[1], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalog1.png').plot()    
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[2], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalog2.png').plot()
            
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcRaster.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhRaster.png').plot({'SpikeRasterPlot.group_trials':True})

            
            if True:
                orr = list(set([MozaikParametrized.idd(s).orientation for s in queries.param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_contrast=100).get_stimuli()]))                
                
                dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=numpy.pi/2,st_contrast=100)   
                KremkowOverviewFigure(dsv,ParameterSet({'neuron' : l4_exc,'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview.png').plot()
                dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[0]),numpy.pi)  for o in orr])],st_contrast=100)   
                KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[0],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview1.png').plot()
                dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[1]),numpy.pi)  for o in orr])],st_contrast=100)   
                KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[1],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview2.png').plot()
                dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[2]),numpy.pi)  for o in orr])],st_contrast=100)   
                KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[2],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview3.png').plot()
                dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[3]),numpy.pi)  for o in orr])],st_contrast=100)   
                KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[3],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview4.png').plot()
                dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[4]),numpy.pi)  for o in orr])],st_contrast=100)   
                KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[4],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview5.png').plot()
                dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[5]),numpy.pi)  for o in orr])],st_contrast=100)   
                KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[5],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview6.png').plot()
                dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[5]),numpy.pi)  for o in orr])],st_contrast=100)   
                KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[6],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview7.png').plot()
                dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[6]),numpy.pi)  for o in orr])],st_contrast=100)   
                KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[7],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview8.png').plot()
                dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[7]),numpy.pi)  for o in orr])],st_contrast=100)   
                KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[8],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview9.png').plot()
                dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[8]),numpy.pi)  for o in orr])],st_contrast=100)   
                KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[9],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview10.png').plot()
                dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(analog_ids[9]),numpy.pi)  for o in orr])],st_contrast=100)   
                KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids[10],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview11.png').plot()
            
            Kremkow_plots.ConductanceAndVmTuningSummary(data_store,ParameterSet({'many' : True}),fig_param={'dpi' : 100,'figsize': (25,16)},plot_file_name='Cond&VMTuning.png').plot()
            
            #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=numpy.pi/2,st_contrast=100)   
            #KremkowOverviewFigure(dsv,ParameterSet({'neuron' : l4_exc,'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview.png').plot()
            
            dsv = param_filter_query(data_store,st_name='NaturalImageWithEyeMovement')            
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : l4_exc,'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='NMOverview.png').plot()
            
            #Kremkow_plots.StimulusResponseComparison(data_store,ParameterSet({'neuron' : l4_exc,'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='StimulusResponseComparison.png').plot()
            Kremkow_plots.OrientationTuningSummary(data_store,ParameterSet({}),fig_param={'dpi' : 100,'figsize': (20,12)},plot_file_name='OrientationTuningSummary.png').plot()
            

            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating')    
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedExcRaster.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedInhRaster.png').plot({'SpikeRasterPlot.group_trials':True})

            #OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0On.png").plot()        
            #OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0Off.png").plot()
    
            # orientation tuning plotting
            dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],value_name='LGNAfferentOrientation')   
            PerNeuronValuePlot(dsv,ParameterSet({"cortical_view" : True}),plot_file_name='ORSet.png').plot()
            
            dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],value_name='orientation preference',analysis_algorithm='PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage',st_contrast=100)    
            PerNeuronValuePlot(dsv,ParameterSet({"cortical_view" : True}),plot_file_name='ORComputed.png').plot()
            
            dsv = param_filter_query(data_store,value_name=['orientation HWHH'],sheet_name=['V1_Exc_L4','V1_Inh_L4'])    
            PerNeuronValueScatterPlot(dsv,ParameterSet({}),plot_file_name='HWHH.png',fig_param={'dpi' : 100,'figsize': (28,12)}).plot({ 'ScatterPlot.x_lim' : (0,90), 'ScatterPlot.y_lim' : (0,90), 'ScatterPlot.identity_line' : True})
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating')    
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="Exc.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="Inh.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            
            
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="Exc1.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[1], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="Exc2.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[2], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="Exc3.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[3], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="Exc4.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[4], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="Exc5.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0On.png").plot()
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0Off.png").plot()
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids())[1], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN1On.png").plot()
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_esyn_ids())[1], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN1Off.png").plot()
            
            # tuninc curves
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate','Analog_F0andF1'])    
            PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids), 'sheet_name' : 'V1_Exc_L4','centered'  : True,'mean' : False}),plot_file_name='TCsExc.png',fig_param={'dpi' : 100,'figsize': (28,12)}).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})
            PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh), 'sheet_name' : 'V1_Inh_L4','centered'  : True,'mean' : False}),plot_file_name='TCsInh.png',fig_param={'dpi' : 100,'figsize': (28,12)}).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})

            dsv = param_filter_query(data_store,st_name='NaturalImageWithEyeMovement')    
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),plot_file_name='NMExc.png',fig_param={'dpi' : 100,'figsize': (28,12)}).plot({'Vm_plot.y_lim' : (-70,-50),'Conductance_plot.y_lim' : (0,50.0)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),plot_file_name='NMInh.png',fig_param={'dpi' : 100,'figsize': (28,12)}).plot({'Vm_plot.y_lim' : (-70,-50),'Conductance_plot.y_lim' : (0,50.0)})
            pylab.show()
