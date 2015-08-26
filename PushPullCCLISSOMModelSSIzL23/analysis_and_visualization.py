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
from lsv1m_paper import *


logger = mozaik.getMozaikLogger()

import psutil
import os
process = psutil.Process(os.getpid())


def memory_usage_psutil():
    # return the memory usage in MB
    return process.memory_percent()
    

def analysis(data_store,analog_ids,analog_ids_inh,analog_ids23=None,analog_ids_inh23=None):
        sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3']))
        exc_sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4','V1_Exc_L2/3']))
        l23_flag = 'V1_Exc_L2/3' in set(sheets)
        
        logger.info('0: ' + str(memory_usage_psutil()))
        
        TrialAveragedFiringRate(param_filter_query(data_store,sheet_name=sheets,st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({})).analyse()
        TrialAveragedFiringRate(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()
        logger.info('1: ' + str(memory_usage_psutil()))
        Irregularity(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()
        PSTH(param_filter_query(data_store),ParameterSet({'bin_length' : 2.0 })).analyse()
        logger.info('2: ' + str(memory_usage_psutil()))
        #SpikeCount(param_filter_query(data_store,sheet_name=exc_sheets),ParameterSet({'bin_length' : 13.0 })).analyse()
        NeuronToNeuronAnalogSignalCorrelations(param_filter_query(data_store,analysis_algorithm='PSTH'),ParameterSet({'convert_nan_to_zero' : True})).analyse()
        
        logger.info('3: ' + str(memory_usage_psutil()))
        
        PopulationMean(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()
                
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name=sheets)    
        GaussianTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name=sheets)   
        Analog_F0andF1(dsv,ParameterSet({})).analyse()
        
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name=sheets)  
        PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
        logger.info('4: ' + str(memory_usage_psutil()))
        #data_store.save()
        
        dsv = param_filter_query(data_store,sheet_name=exc_sheets)
        ActionPotentialRemoval(dsv,ParameterSet({'window_length': 5.0})).analyse()
        
        logger.info('5: ' + str(memory_usage_psutil()))
        
        GSTA(param_filter_query(data_store,sheet_name='V1_Exc_L4',st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100),ParameterSet({'neurons' : list(analog_ids), 'length' : 250.0 }),tags=['GSTA']).analyse()
        GSTA(param_filter_query(data_store,sheet_name='V1_Inh_L4',st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100),ParameterSet({'neurons' : list(analog_ids_inh), 'length' : 250.0 }),tags=['GSTA']).analyse()
        
        logger.info('6: ' + str(memory_usage_psutil()))
        
        if l23_flag:
            GSTA(param_filter_query(data_store,sheet_name='V1_Exc_L2/3',st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100),ParameterSet({'neurons' : list(analog_ids23), 'length' : 250.0 }),tags=['GSTA']).analyse()
            GSTA(param_filter_query(data_store,sheet_name='V1_Inh_L2/3',st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100),ParameterSet({'neurons' : list(analog_ids_inh23), 'length' : 250.0 }),tags=['GSTA']).analyse()

        GSTA(param_filter_query(data_store,sheet_name='V1_Exc_L4',st_name='NaturalImageWithEyeMovement'),ParameterSet({'neurons' : list(analog_ids), 'length' : 250.0 }),tags=['GSTA']).analyse()
        GSTA(param_filter_query(data_store,sheet_name='V1_Inh_L4',st_name='NaturalImageWithEyeMovement'),ParameterSet({'neurons' : list(analog_ids_inh), 'length' : 250.0 }),tags=['GSTA']).analyse()

        if l23_flag:
            GSTA(param_filter_query(data_store,sheet_name='V1_Exc_L2/3',st_name='NaturalImageWithEyeMovement'),ParameterSet({'neurons' : list(analog_ids23), 'length' : 250.0 }),tags=['GSTA']).analyse()
            GSTA(param_filter_query(data_store,sheet_name='V1_Inh_L2/3',st_name='NaturalImageWithEyeMovement'),ParameterSet({'neurons' : list(analog_ids_inh23), 'length' : 250.0 }),tags=['GSTA']).analyse()
        
        logger.info('7: ' + str(memory_usage_psutil()))
        
        dsv = param_filter_query(data_store,st_name='InternalStimulus',st_direct_stimulation_name="None")
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

        pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Inh_L4',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(ECond)',st_direct_stimulation_name='None').get_analysis_result()[0]
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Inh_L4',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Exc_Cond')
        SubtractPNVfromPNVS(pnv,dsv,ParameterSet({})).analyse()

        pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Inh_L4',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(ICond)',st_direct_stimulation_name='None').get_analysis_result()[0]
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Inh_L4',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Inh_Cond')
        SubtractPNVfromPNVS(pnv,dsv,ParameterSet({})).analyse()

        pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Inh_L4',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(VM)',st_direct_stimulation_name='None').get_analysis_result()[0]
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Inh_L4',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Vm')
        SubtractPNVfromPNVS(pnv,dsv,ParameterSet({})).analyse()



        logger.info('8: ' + str(memory_usage_psutil()))
        if l23_flag:

            pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Exc_L2/3',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(ECond)',st_direct_stimulation_name='None').get_analysis_result()[0]
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Exc_L2/3',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Exc_Cond')
            SubtractPNVfromPNVS(pnv,dsv,ParameterSet({})).analyse()

            pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Exc_L2/3',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(ICond)',st_direct_stimulation_name='None').get_analysis_result()[0]
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Exc_L2/3',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Inh_Cond')
            SubtractPNVfromPNVS(pnv,dsv,ParameterSet({})).analyse()
            
            pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Exc_L2/3',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(VM)',st_direct_stimulation_name='None').get_analysis_result()[0]
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Exc_L2/3',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Vm')
            SubtractPNVfromPNVS(pnv,dsv,ParameterSet({})).analyse()
            
            pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Inh_L2/3',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(ECond)',st_direct_stimulation_name='None').get_analysis_result()[0]
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Inh_L2/3',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Exc_Cond')
            SubtractPNVfromPNVS(pnv,dsv,ParameterSet({})).analyse()

            pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Inh_L2/3',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(ICond)',st_direct_stimulation_name='None').get_analysis_result()[0]
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Inh_L2/3',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Inh_Cond')
            SubtractPNVfromPNVS(pnv,dsv,ParameterSet({})).analyse()

            pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Inh_L2/3',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(VM)',st_direct_stimulation_name='None').get_analysis_result()[0]
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Inh_L2/3',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Vm')
            SubtractPNVfromPNVS(pnv,dsv,ParameterSet({})).analyse()
            
            
        
        logger.info('9: ' + str(memory_usage_psutil()))
        dsv = queries.param_filter_query(data_store,y_axis_name='spike count (bin=13.0)')   
        mozaik.analysis.analysis.TrialToTrialFanoFactorOfAnalogSignal(dsv,ParameterSet({})).analyse()
        logger.info('10: ' + str(memory_usage_psutil()))
        TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(data_store,sheet_name='V1_Exc_L4',st_name="NaturalImageWithEyeMovement",analysis_algorithm='ActionPotentialRemoval'),ParameterSet({'neurons' : list(analog_ids)})).analyse()
        TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(data_store,sheet_name='V1_Exc_L4',st_name="NaturalImageWithEyeMovement",analysis_algorithm='PSTH'),ParameterSet({'neurons' : list(analog_ids)})).analyse()
        logger.info('11: ' + str(memory_usage_psutil()))
        if l23_flag:
            TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(data_store,sheet_name='V1_Exc_L2/3',st_name="NaturalImageWithEyeMovement",analysis_algorithm='ActionPotentialRemoval'),ParameterSet({'neurons' : list(analog_ids23)})).analyse()
            TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(data_store,sheet_name='V1_Exc_L2/3',st_name="NaturalImageWithEyeMovement",analysis_algorithm='PSTH'),ParameterSet({'neurons' : list(analog_ids23)})).analyse()
        logger.info('12: ' + str(memory_usage_psutil()))

        dsv = param_filter_query(data_store,analysis_algorithm='ActionPotentialRemoval')
        TrialVariability(data_store,ParameterSet({'vm': True,  'cond_exc': False, 'cond_inh': False})).analyse()
        TrialMean(data_store,ParameterSet({'vm': True,  'cond_exc': False, 'cond_inh': False})).analyse()
        logger.info('13: ' + str(memory_usage_psutil()))
        ModulationRatio(param_filter_query(data_store,sheet_name=exc_sheets,st_contrast=[100]),ParameterSet({})).analyse()
       
        logger.info('14: ' + str(memory_usage_psutil()))
        
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate')    
        CircularVarianceOfTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
        
        logger.info('15: ' + str(memory_usage_psutil()))
        
        data_store.save()

def perform_analysis_and_visualization_conn(data_store):
    analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()
    spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()

    analog_ids23 = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh23 = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_esyn_ids()
    spike_ids23 = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh23 = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_spike_train_ids()

    NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
    l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')[0]
    l23_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L2/3')[0]
    l4_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L4')[0]
    l23_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L2/3')[0]


    l4_exc_or_many = numpy.array(spike_ids)[numpy.nonzero(numpy.array([circular_dist(l4_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids]) < 0.1)[0]]
    l23_exc_or_many = numpy.array(spike_ids23)[numpy.nonzero(numpy.array([circular_dist(l23_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids23]) < 0.1)[0]]
    l4_inh_or_many = numpy.array(spike_ids_inh)[numpy.nonzero(numpy.array([circular_dist(l4_inh_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids_inh]) < 0.1)[0]]
    l23_inh_or_many = numpy.array(spike_ids_inh23)[numpy.nonzero(numpy.array([circular_dist(l23_inh_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids_inh23]) < 0.1)[0]]

    dsv = param_filter_query(data_store,identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation')
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_exc_or_many[0], 'reversed' : True,'sheet_name' : 'V1_Exc_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='ExcConnections1.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_exc_or_many[1], 'reversed' : True,'sheet_name' : 'V1_Exc_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='ExcConnections2.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_exc_or_many[2], 'reversed' : True,'sheet_name' : 'V1_Exc_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='ExcConnections3.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_exc_or_many[3], 'reversed' : True,'sheet_name' : 'V1_Exc_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='ExcConnections4.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_exc_or_many[4], 'reversed' : True,'sheet_name' : 'V1_Exc_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='ExcConnections5.png').plot()

    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_inh_or_many[0], 'reversed' : True,'sheet_name' : 'V1_Inh_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='InhConnections1.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_inh_or_many[1], 'reversed' : True,'sheet_name' : 'V1_Inh_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='InhConnections2.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_inh_or_many[2], 'reversed' : True,'sheet_name' : 'V1_Inh_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='InhConnections3.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_inh_or_many[3], 'reversed' : True,'sheet_name' : 'V1_Inh_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='InhConnections4.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l4_inh_or_many[4], 'reversed' : True,'sheet_name' : 'V1_Inh_L4'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='InhConnections5.png').plot()
    
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l23_exc_or_many[0], 'reversed' : True,'sheet_name' : 'V1_Exc_L2/3'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='L23ExcConnections1.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l23_exc_or_many[1], 'reversed' : True,'sheet_name' : 'V1_Exc_L2/3'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='L23ExcConnections2.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l23_exc_or_many[2], 'reversed' : True,'sheet_name' : 'V1_Exc_L2/3'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='L23ExcConnections3.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l23_exc_or_many[3], 'reversed' : True,'sheet_name' : 'V1_Exc_L2/3'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='L23ExcConnections4.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l23_exc_or_many[4], 'reversed' : True,'sheet_name' : 'V1_Exc_L2/3'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='L23ExcConnections5.png').plot()    
    
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l23_inh_or_many[0], 'reversed' : True,'sheet_name' : 'V1_Inh_L2/3'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='L23InhConnections1.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l23_inh_or_many[1], 'reversed' : True,'sheet_name' : 'V1_Inh_L2/3'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='L23InhConnections2.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l23_inh_or_many[2], 'reversed' : True,'sheet_name' : 'V1_Inh_L2/3'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='L23InhConnections3.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l23_inh_or_many[3], 'reversed' : True,'sheet_name' : 'V1_Inh_L2/3'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='L23InhConnections4.png').plot()
    mozaik.visualization.plotting.ConnectivityPlot(data_store,ParameterSet({'neuron' : l23_inh_or_many[4], 'reversed' : True,'sheet_name' : 'V1_Inh_L2/3'}),pnv_dsv=dsv,fig_param={'dpi' : 100,'figsize': (34,12)},plot_file_name='L23InhConnections5.png').plot()



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
    l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')[0]
    l4_exc_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Exc_L4')
    l4_exc = analog_ids[numpy.argmin([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_exc_or.get_value_by_id(analog_ids),l4_exc_phase[0].get_value_by_id(analog_ids))])]
    l4_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L4')[0]
    l4_inh_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Inh_L4')
    l4_inh = analog_ids_inh[numpy.argmin([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_inh_or.get_value_by_id(analog_ids_inh),l4_inh_phase[0].get_value_by_id(analog_ids_inh))])]

    l23_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L2/3')[0]
    l23_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L2/3')[0]

    l4_exc_or_many = numpy.array(spike_ids)[numpy.nonzero(numpy.array([circular_dist(l4_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids]) < 0.1)[0]]
    l23_exc_or_many = numpy.array(spike_ids23)[numpy.nonzero(numpy.array([circular_dist(l23_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids23]) < 0.1)[0]]

    l4_exc_or_many_analog = numpy.array(analog_ids)[numpy.nonzero(numpy.array([circular_dist(l4_exc_or.get_value_by_id(i),0,numpy.pi)  for i in analog_ids]) < 0.1)[0]]
    l23_exc_or_many_analog = numpy.array(analog_ids23)[numpy.nonzero(numpy.array([circular_dist(l23_exc_or.get_value_by_id(i),0,numpy.pi)  for i in analog_ids23]) < 0.1)[0]]
    
    l4_exc_or_many_analog_inh = numpy.array(analog_ids_inh)[numpy.nonzero(numpy.array([circular_dist(l4_inh_or.get_value_by_id(i),0,numpy.pi)  for i in analog_ids_inh]) < 0.1)[0]]
    l23_exc_or_many_analog_inh = numpy.array(analog_ids_inh23)[numpy.nonzero(numpy.array([circular_dist(l23_inh_or.get_value_by_id(i),0,numpy.pi)  for i in analog_ids_inh23]) < 0.1)[0]]
    


    
    if True:
        if True:
           analysis(data_store,l4_exc_or_many_analog,l4_exc_or_many_analog_inh,l23_exc_or_many_analog,l23_exc_or_many_analog_inh)
        
        if True:
            # self sustained plotting
            dsv = param_filter_query(data_store,st_name=['InternalStimulus'])    
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalog.png').plot()
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhAnalog.png').plot()    
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSL23ExcAnalog.png').plot()
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSL23InhAnalog.png').plot()    
            
            
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neurons' : spike_ids23,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcRasterL23.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neurons' : spike_ids_inh23,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhRasterL23.png').plot({'SpikeRasterPlot.group_trials':True})

            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating')
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc_or_many_analog[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc1.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc_or_many_analog[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc2.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc_or_many_analog[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc3.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc_or_many_analog[3], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc4.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc_or_many_analog[4], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc5.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc_or_many_analog[5], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc6.png").plot({'Vm_plot.y_lim' : (-80,-50)})

            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc_or_many_analog[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL231.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc_or_many_analog[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL232.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc_or_many_analog[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL233.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc_or_many_analog[3], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL234.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc_or_many_analog[4], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL235.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc_or_many_analog[5], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL236.png").plot({'Vm_plot.y_lim' : (-80,-50)})


            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh1.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh2.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="InhL231.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="InhL232.png").plot({'Vm_plot.y_lim' : (-80,-50)})


            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedExcRaster.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedInhRaster.png').plot({'SpikeRasterPlot.group_trials':True})

            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : list(analog_ids),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedExcRasterWithAnalog.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : list(analog_ids_inh),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedInhRasterWithAnalog.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neurons' : list(analog_ids23),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedExcL23RasterWithAnalog.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neurons' : list(analog_ids_inh23),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedInhL23RasterWithAnalog.png').plot({'SpikeRasterPlot.group_trials':True})



        #dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L4',st_contrast=100)    
        #FanoFactor_Baudot_et_al(dsv,ParameterSet({}),plot_file_name='FanoFactorL4.png').plot()
        
        #dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L2/3',st_contrast=100)    
        #FanoFactor_Baudot_et_al(dsv,ParameterSet({}),plot_file_name='FanoFactorL23.png').plot()
        if True:
            pref_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
            ort_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
            pref_fr_50 = 0#param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[50],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
            ort_fr_50 = 0#param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[50],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
            spont = param_filter_query(data_store,st_name='InternalStimulus',sheet_name=['V1_Exc_L4'],analysis_algorithm=['TrialAveragedFiringRate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)

            pylab.figure()
            pylab.bar([1,2,3,4,5],[numpy.mean(pref_fr_100),numpy.mean(ort_fr_100),numpy.mean(pref_fr_50),numpy.mean(ort_fr_50),numpy.mean(spont)])
            pylab.xticks([1.4,2.4,3.4,4.4,5.4],['PREF100','ORT100','PREF50','ORT50','SPONT'])
            pylab.xlim([0.8,6.0])
            pylab.ylabel("Firing rate")

            pylab.savefig(Global.root_directory+"Orientation_responseL4.png")
            
            
            pref_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L2/3'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)
            ort_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L2/3'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)
            pref_fr_50 = 0#param_filter_query(data_store,sheet_name=['V1_Exc_L2/3'],st_contrast=[50],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)
            ort_fr_50 = 0#param_filter_query(data_store,sheet_name=['V1_Exc_L2/3'],st_contrast=[50],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)
            spont = param_filter_query(data_store,st_name='InternalStimulus',sheet_name=['V1_Exc_L2/3'],analysis_algorithm=['TrialAveragedFiringRate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)

            pylab.figure()
            pylab.bar([1,2,3,4,5],[numpy.mean(pref_fr_100),numpy.mean(ort_fr_100),numpy.mean(pref_fr_50),numpy.mean(ort_fr_50),numpy.mean(spont)])
            pylab.xticks([1.4,2.4,3.4,4.4,5.4],['PREF100','ORT100','PREF50','ORT50','SPONT'])
            pylab.xlim([0.8,6.0])
            pylab.ylabel("Firing rate")
            
            pylab.savefig(Global.root_directory+"Orientation_responseL23.png")
        
        TrialToTrialVariabilityComparison(data_store,ParameterSet({'sheet_name1' : 'V1_Exc_L4','sheet_name2' : 'V1_Exc_L2/3','data_dg' : 0.93 , 'data_ni' : 1.19}),fig_param={'dpi' : 200,'figsize': (10,5)}).plot()
        
        StimulusResponseComparison(data_store,ParameterSet({'neuron' : l4_exc_or_many_analog[0],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparison1.png').plot()
        StimulusResponseComparison(data_store,ParameterSet({'neuron' : l4_exc_or_many_analog[1],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparison2.png').plot()
        StimulusResponseComparison(data_store,ParameterSet({'neuron' : l4_exc_or_many_analog[2],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparison3.png').plot()
        StimulusResponseComparison(data_store,ParameterSet({'neuron' : l4_exc_or_many_analog[3],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparison4.png').plot()

        StimulusResponseComparison(data_store,ParameterSet({'neuron' : l23_exc_or_many_analog[0],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparisonL231.png').plot()
        StimulusResponseComparison(data_store,ParameterSet({'neuron' : l23_exc_or_many_analog[1],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparisonL232.png').plot()
        StimulusResponseComparison(data_store,ParameterSet({'neuron' : l23_exc_or_many_analog[2],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparisonL233.png').plot()
        StimulusResponseComparison(data_store,ParameterSet({'neuron' : l23_exc_or_many_analog[3],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparisonL234.png').plot()

        import MRfig
        MRfig.MRfig(param_filter_query(data_store,sheet_name=['V1_Exc_L2/3','V1_Exc_L4'],st_contrast=[100],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'V1_Exc_L2/3'}),plot_file_name='MR.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()
        
        SpontActOverview(data_store,ParameterSet({'l4_exc_neuron' : analog_ids[0], 'l4_inh_neuron' : analog_ids_inh[0],'l23_exc_neuron' : analog_ids23[0],'l23_inh_neuron' : analog_ids_inh23[0]}),plot_file_name='SpontActOverview.png', fig_param={'dpi' : 200,'figsize': (15,12)}).plot()
        SpontStatisticsOverview(data_store,ParameterSet({}), fig_param={'dpi' : 200,'figsize': (15,10)},plot_file_name='SpontStatisticsOverview.png').plot()
        
        TrialCrossCorrelationAnalysis(data_store,ParameterSet({'neurons' : list(l23_exc_or_many_analog),'window_length' : 83, 'sheet_name' : 'V1_Exc_L2/3'}),fig_param={"dpi" : 100,"figsize": (25,12)},plot_file_name="trial-to-trial-cross-correlation23.png").plot()
        TrialCrossCorrelationAnalysis(data_store,ParameterSet({'neurons' : list(l4_exc_or_many_analog),'window_length' : 83, 'sheet_name' : 'V1_Exc_L4'}),fig_param={"dpi" : 100,"figsize": (25,12)},plot_file_name="trial-to-trial-cross-correlation4.png").plot()

        
        dsv = param_filter_query(data_store,value_name=['Mean(ECond)','Mean(ICond)','Mean(VM)','Mean(ECond)/Mean(ICond)','Firing rate'])   
        PerNeuronValuePlot(dsv,ParameterSet({"cortical_view" : False}),plot_file_name='Histograms.png',fig_param={'dpi' : 100,'figsize': (25,12)}).plot({'*.title' : None,'HistogramPlot.plot[3,1].x_scale' : 'log','HistogramPlot.plot[3,1].log' : True,'HistogramPlot.plot[3,0].x_scale' : 'log','HistogramPlot.plot[3,0].log' : True})
        
        dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L4')    
        MeanVsVarainceOfVM(dsv,ParameterSet({'neurons' : list(l4_exc_or_many_analog[:5])}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='TrialToTrialMeanVsVarianceOfVML4.png').plot()

        dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L2/3')    
        MeanVsVarainceOfVM(dsv,ParameterSet({'neurons' : list(l23_exc_or_many_analog[:5])}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='TrialToTrialMeanVsVarianceOfVML23.png').plot()
        
        
        #SNRAnalysis(data_store,ParameterSet({"neuron" : l4_exc}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='SNR1.png').plot()                        
      
    dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'],value_name='LGNAfferentOrientation')   
    dsv.print_content(full_ADS=True)
    PerNeuronValuePlot(dsv,ParameterSet({'cortical_view' : True}),plot_file_name='Map.png').plot({'*.colorbar' : False,'*.x_axis' : None, '*.y_axis' : None})


def perform_analysis_and_visualization_stc(data_store):
    analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()
    spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()

    analog_ids23 = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh23 = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_esyn_ids()
    spike_ids23 = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh23 = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_spike_train_ids()

    NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
    l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')[0]
    l23_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L2/3')[0]

    l4_exc_or_many = numpy.array(spike_ids)[numpy.nonzero(numpy.array([circular_dist(l4_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids]) < 0.1)[0]]
    l23_exc_or_many = numpy.array(spike_ids23)[numpy.nonzero(numpy.array([circular_dist(l23_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids23]) < 0.1)[0]]

    idx23 = data_store.get_sheet_indexes(sheet_name='V1_Exc_L2/3',neuron_ids=l23_exc_or_many)
    idx4 = data_store.get_sheet_indexes(sheet_name='V1_Exc_L4',neuron_ids=l4_exc_or_many)

    x = data_store.get_neuron_postions()['V1_Exc_L4'][0][idx4]
    y = data_store.get_neuron_postions()['V1_Exc_L4'][1][idx4]
    center4 = l4_exc_or_many[numpy.nonzero(numpy.sqrt(numpy.multiply(x,x)+numpy.multiply(y,y)) < 0.1)[0]]
    x = data_store.get_neuron_postions()['V1_Exc_L2/3'][0][idx23]
    y = data_store.get_neuron_postions()['V1_Exc_L2/3'][1][idx23]
    center23 = l23_exc_or_many[numpy.nonzero(numpy.sqrt(numpy.multiply(x,x)+numpy.multiply(y,y)) < 0.1)[0]]
    
    analog_center4 = set(center4).intersection(analog_ids)
    analog_center23 = set(center23).intersection(analog_ids23)

    if False:
        l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')
        l4_exc_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Exc_L4')
        l4_exc = analog_ids[numpy.argmin([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_exc_or[0].get_value_by_id(analog_ids),l4_exc_phase[0].get_value_by_id(analog_ids))])]
        l4_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L4')
        l4_inh_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Inh_L4')
        l4_inh = analog_ids_inh[numpy.argmin([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_inh_or[0].get_value_by_id(analog_ids_inh),l4_inh_phase[0].get_value_by_id(analog_ids_inh))])]

        #analysis(data_store,analog_ids,analog_ids_inh,analog_ids23,analog_ids_inh23)

    TrialAveragedFiringRate(param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'],st_name='DriftingSinusoidalGratingDisk'),ParameterSet({})).analyse()
        
    print len(l4_exc_or_many)
    print len(l23_exc_or_many)
    
    dsv = param_filter_query(data_store,st_name='DriftingSinusoidalGratingDisk',analysis_algorithm=['TrialAveragedFiringRate'])    
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'radius', 'neurons': list(center4), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='SizeTuningExcL4.png',fig_param={'dpi' : 100,'figsize': (32,7)}).plot()
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'radius', 'neurons': list(center23), 'sheet_name' : 'V1_Exc_L2/3','centered'  : False,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='SizeTuningExcL23.png',fig_param={'dpi' : 100,'figsize': (32,7)}).plot()        
    
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'radius', 'neurons': list(center4), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : True, 'polar' : False, 'pool'  : False}),plot_file_name='SizeTuningExcL4M.png',fig_param={'dpi' : 100,'figsize': (32,7)}).plot()
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'radius', 'neurons': list(center23), 'sheet_name' : 'V1_Exc_L2/3','centered'  : False,'mean' : True, 'polar' : False, 'pool'  : False}),plot_file_name='SizeTuningExcL23M.png',fig_param={'dpi' : 100,'figsize': (32,7)}).plot()        
    data_store.save()
    
    if True:
        dsv = param_filter_query(data_store,st_name=['DriftingSinusoidalGratingDisk'])    
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_center4[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_ExcL4_1.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_center4[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_ExcL4_2.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_center4[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_ExcL4_3.png').plot()
    
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_center23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_ExcL23_1.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_center23[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_ExcL23_2.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_center23[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_ExcL23_3.png').plot()
    


def perform_analysis_and_visualization_octc(data_store):
    analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()
    spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()

    analog_ids23 = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh23 = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_esyn_ids()
    spike_ids23 = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh23 = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_spike_train_ids()

    NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
    l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')[0]
    l23_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L2/3')[0]
    l4_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L4')[0]
    l23_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L2/3')[0]


    l4_exc_or_many = numpy.array(spike_ids)[numpy.nonzero(numpy.array([circular_dist(l4_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids]) < 0.1)[0]]
    l23_exc_or_many = numpy.array(spike_ids23)[numpy.nonzero(numpy.array([circular_dist(l23_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids23]) < 0.1)[0]]
    l4_inh_or_many = numpy.array(spike_ids_inh)[numpy.nonzero(numpy.array([circular_dist(l4_inh_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids_inh]) < 0.1)[0]]
    l23_inh_or_many = numpy.array(spike_ids_inh23)[numpy.nonzero(numpy.array([circular_dist(l23_inh_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids_inh23]) < 0.1)[0]]


    idx23 = data_store.get_sheet_indexes(sheet_name='V1_Exc_L2/3',neuron_ids=l23_exc_or_many)
    idx4 = data_store.get_sheet_indexes(sheet_name='V1_Exc_L4',neuron_ids=l4_exc_or_many)
    idx_inh23 = data_store.get_sheet_indexes(sheet_name='V1_Inh_L2/3',neuron_ids=l23_inh_or_many)
    idx_inh4 = data_store.get_sheet_indexes(sheet_name='V1_Inh_L4',neuron_ids=l4_inh_or_many)

    x = data_store.get_neuron_postions()['V1_Exc_L4'][0][idx4]
    y = data_store.get_neuron_postions()['V1_Exc_L4'][1][idx4]
    center4 = l4_exc_or_many[numpy.nonzero(numpy.sqrt(numpy.multiply(x,x)+numpy.multiply(y,y)) < 0.1)[0]]
    x = data_store.get_neuron_postions()['V1_Exc_L2/3'][0][idx23]
    y = data_store.get_neuron_postions()['V1_Exc_L2/3'][1][idx23]
    center23 = l23_exc_or_many[numpy.nonzero(numpy.sqrt(numpy.multiply(x,x)+numpy.multiply(y,y)) < 0.1)[0]]

    
    analog_center4 = set(center4).intersection(analog_ids)
    analog_center23 = set(center23).intersection(analog_ids23)

    TrialAveragedFiringRate(param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'],st_name='DriftingSinusoidalGratingCenterSurroundStimulus'),ParameterSet({})).analyse()
    
    dsv = param_filter_query(data_store,st_name='DriftingSinusoidalGratingCenterSurroundStimulus',analysis_algorithm=['TrialAveragedFiringRate'])    
    
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'surround_orientation', 'neurons': list(center4), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : True, 'polar' : False, 'pool'  : False}),plot_file_name='OCTCExcL4M.png',fig_param={'dpi' : 100,'figsize': (32,15)}).plot()
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'surround_orientation', 'neurons': list(center23), 'sheet_name' : 'V1_Exc_L2/3','centered'  : False,'mean' : True, 'polar' : False, 'pool'  : False}),plot_file_name='OCTCExcL23M.png',fig_param={'dpi' : 100,'figsize': (32,15)}).plot()        
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'surround_orientation', 'neurons': list(center4), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='OCTCExcL4.png',fig_param={'dpi' : 100,'figsize': (32,15)}).plot()
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'surround_orientation', 'neurons': list(center23), 'sheet_name' : 'V1_Exc_L2/3','centered'  : False,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='OCTCExcL23.png',fig_param={'dpi' : 100,'figsize': (32,15)}).plot()        
   
   
    print data_store.get_neuron_postions()['V1_Inh_L4'][0]
    print idx_inh4
    x = data_store.get_neuron_postions()['V1_Inh_L4'][0][idx_inh4]
    y = data_store.get_neuron_postions()['V1_Inh_L4'][1][idx_inh4]
    center4_inh = l4_inh_or_many[numpy.nonzero(numpy.sqrt(numpy.multiply(x,x)+numpy.multiply(y,y)) < 0.1)[0]]
    x = data_store.get_neuron_postions()['V1_Inh_L2/3'][0][idx_inh23]
    y = data_store.get_neuron_postions()['V1_Inh_L2/3'][1][idx_inh23]
    center23_inh = l23_inh_or_many[numpy.nonzero(numpy.sqrt(numpy.multiply(x,x)+numpy.multiply(y,y)) < 0.1)[0]]
    
    if True:
        dsv = param_filter_query(data_store,st_name=['DriftingSinusoidalGratingCenterSurroundStimulus'],st_surround_orientation=[0,numpy.pi/2])    
    
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : list(center23)[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_ExcL23_1.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : list(center23)[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_ExcL23_2.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : list(center23)[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_ExcL23_3.png').plot()

        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : list(center23_inh)[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_InhL23_1.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : list(center23_inh)[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_InhL23_2.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : list(center23_inh)[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_InhL23_3.png').plot()
        


def perform_analysis_and_visualization(data_store):
    sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3']))
    exc_sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4','V1_Exc_L2/3']))
    l23_flag = 'V1_Exc_L2/3' in set(sheets)

    NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
    
    analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()
    spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()
    
    if l23_flag:
        analog_ids23 = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_esyn_ids()
        analog_ids_inh23 = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_esyn_ids()
        spike_ids23 = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_spike_train_ids()
        spike_ids_inh23 = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_spike_train_ids()
    else:
        analog_ids23 = None
        analog_ids_inh23 = None

    data_store.print_content(full_ADS=True)

    if l23_flag:
        l23_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L2/3')[0]
        l23_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L2/3')[0]    
        
    
    l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')
    l4_exc_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Exc_L4')
    l4_exc = analog_ids[numpy.argmin([circular_dist(o,numpy.pi/2,numpy.pi)  for (o,p) in zip(l4_exc_or[0].get_value_by_id(analog_ids),l4_exc_phase[0].get_value_by_id(analog_ids))])]
    l4_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L4')
    l4_inh_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Inh_L4')
    l4_inh = analog_ids_inh[numpy.argmin([circular_dist(o,numpy.pi/2,numpy.pi)  for (o,p) in zip(l4_inh_or[0].get_value_by_id(analog_ids_inh),l4_inh_phase[0].get_value_by_id(analog_ids_inh))])]
    l4_exc_or_many = numpy.array(l4_exc_or[0].ids)[numpy.nonzero(numpy.array([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_exc_or[0].values,l4_exc_phase[0].values)]) < 0.1)[0]]
    l23_exc_or_many = numpy.array(l23_exc_or.ids)[numpy.nonzero(numpy.array([circular_dist(o,0,numpy.pi)  for o in l23_exc_or.values]) < 0.1)[0]]
    
    l4_exc_or_many = list(set(l4_exc_or_many) &  set(spike_ids))
    l23_exc_or_many = list(set(l23_exc_or_many) &  set(spike_ids23))
    
    
    
    orr = list(set([MozaikParametrized.idd(s).orientation for s in queries.param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_contrast=100).get_stimuli()]))                
        
    l4_exc_or_many_analog = numpy.array(analog_ids)[numpy.nonzero(numpy.array([circular_dist(l4_exc_or[0].get_value_by_id(i),0,numpy.pi)  for i in analog_ids]) < 0.1)[0]]
    l4_inh_or_many_analog = numpy.array(analog_ids_inh)[numpy.nonzero(numpy.array([circular_dist(l4_inh_or[0].get_value_by_id(i),0,numpy.pi)  for i in analog_ids_inh]) < 0.1)[0]]
    
    if l23_flag:
        l23_inh_or_many_analog = numpy.array(analog_ids_inh23)[numpy.nonzero(numpy.array([circular_dist(l23_inh_or.get_value_by_id(i),0,numpy.pi)  for i in analog_ids_inh23]) < 0.1)[0]]
        l23_exc_or_many_analog = numpy.array(analog_ids23)[numpy.nonzero(numpy.array([circular_dist(l23_exc_or.get_value_by_id(i),0,numpy.pi)  for i in analog_ids23]) < 0.1)[0]]
            
    if True:
       analysis(data_store,analog_ids,analog_ids_inh,analog_ids23,analog_ids_inh23)
    
    a = l4_exc_or[0].ids
    b = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0],ads_unique=True).get_analysis_result()[0].ids
    
    logger.info(str(len(a)))
    logger.info(str(len(b)))
    
    logger.info(str(a))
    logger.info(str(b))
    
    
    if True:
            pref_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
            ort_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
            pref_fr_50 = 0#param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[50],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
            ort_fr_50 = 0#param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[50],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
            spont = param_filter_query(data_store,st_name='InternalStimulus',sheet_name=['V1_Exc_L4'],analysis_algorithm=['TrialAveragedFiringRate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)

            pylab.figure()
            pylab.bar([1,2,3,4,5],[numpy.mean(pref_fr_100),numpy.mean(ort_fr_100),numpy.mean(pref_fr_50),numpy.mean(ort_fr_50),numpy.mean(spont)])
            pylab.xticks([1.4,2.4,3.4,4.4,5.4],['PREF100','ORT100','PREF50','ORT50','SPONT'])
            pylab.xlim([0.8,6.0])
            pylab.ylabel("Firing rate")

            pylab.savefig(Global.root_directory+"Orientation_responseL4.png")
            
            
            pref_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L2/3'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)
            ort_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L2/3'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)
            pref_fr_50 = 0#param_filter_query(data_store,sheet_name=['V1_Exc_L2/3'],st_contrast=[50],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)
            ort_fr_50 = 0#param_filter_query(data_store,sheet_name=['V1_Exc_L2/3'],st_contrast=[50],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)
            spont = param_filter_query(data_store,st_name='InternalStimulus',sheet_name=['V1_Exc_L2/3'],analysis_algorithm=['TrialAveragedFiringRate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)

            pylab.figure()
            pylab.bar([1,2,3,4,5],[numpy.mean(pref_fr_100),numpy.mean(ort_fr_100),numpy.mean(pref_fr_50),numpy.mean(ort_fr_50),numpy.mean(spont)])
            pylab.xticks([1.4,2.4,3.4,4.4,5.4],['PREF100','ORT100','PREF50','ORT50','SPONT'])
            pylab.xlim([0.8,6.0])
            pylab.ylabel("Firing rate")
            
            pylab.savefig(Global.root_directory+"Orientation_responseL23.png")

    
    if l23_flag and True:
        LSV1MReponseOverview(data_store,ParameterSet({'l4_exc_neuron' : l4_exc_or_many_analog[0], 'l4_inh_neuron' : l4_inh_or_many_analog[0],'l23_exc_neuron' : l23_exc_or_many_analog[0], 'l23_inh_neuron' : l23_inh_or_many_analog[0]}),fig_param={'dpi' : 70,'figsize': (30,15)},plot_file_name='SingleCellOverview.png').plot({'*.x_ticks':[0,2.0], "*.title" : None})
        
        dsv = param_filter_query(data_store,st_name='NaturalImageWithEyeMovement')
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[16], 'sheet_activity' : {}, 'spontaneous' : False}), fig_param={'dpi' : 100,'figsize': (10,6)},plot_file_name="NatExcL4_n16.png").plot({'Vm_plot.y_lim' : (-80,-50)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[15], 'sheet_activity' : {}, 'spontaneous' : False}), fig_param={'dpi' : 100,'figsize': (10,6)},plot_file_name="NatExcL23_n15.png").plot({'Vm_plot.y_lim' : (-80,-50)})

        #SpontRasterOverview(data_store,ParameterSet({}),fig_param={'dpi' : 70,'figsize': (20,20)},plot_file_name='SpontRasterOverview.png').plot()
        #SpontAnalogOverview(data_store,ParameterSet({}),fig_param={'dpi' : 70,'figsize': (20,20)},plot_file_name='SpontAnalogOverview.png').plot()
      
    
    if True: # PLOTTING
            activity_plot_param =    {
                   'frame_rate' : 5,  
                   'bin_width' : 5.0, 
                   'scatter' :  True,
                   'resolution' : 0
            }   
            
            Kremkow_plots.OrientationTuningSummary(data_store,ParameterSet({'exc_sheet_name': 'V1_Exc_L4','inh_sheet_name': 'V1_Inh_L4'}),fig_param={'dpi' : 100,'figsize': (20,12)},plot_file_name='OrientationTuningSummaryL4.png').plot()            
            if l23_flag:
                Kremkow_plots.OrientationTuningSummary(data_store,ParameterSet({'exc_sheet_name': 'V1_Exc_L2/3','inh_sheet_name': 'V1_Inh_L2/3'}),fig_param={'dpi' : 100,'figsize': (20,12)},plot_file_name='OrientationTuningSummaryL23.png').plot()            
            
            #Kremkow_plots.ConductanceAndVmTuningSummary(data_store,ParameterSet({'many' : True}),fig_param={'dpi' : 100,'figsize': (25,16)},plot_file_name='Cond&VMTuning.png').plot()
            
            # self sustained plotting
            dsv = param_filter_query(data_store,st_name=['InternalStimulus'],st_direct_stimulation_name="None")    
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalog.png').plot()
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhAnalog.png').plot()    
            
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})
            if l23_flag:
                RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neurons' : spike_ids23,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcRasterL23.png').plot({'SpikeRasterPlot.group_trials':True})
                RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neurons' : spike_ids_inh23,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhRasterL23.png').plot({'SpikeRasterPlot.group_trials':True})
            
            if True:
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

                if l23_flag:
                    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
                    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids23[0],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='L23ExcOverview1.png').plot()
                    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
                    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids23[1],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='L23ExcOverview2.png').plot()
                    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
                    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids23[2],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='L23ExcOverview3.png').plot()
                    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
                    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids23[3],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='L23ExcOverview4.png').plot()
                    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
                    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids23[4],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='L23ExcOverview5.png').plot()
                    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
                    KremkowOverviewFigure(dsv,ParameterSet({'neuron' : analog_ids23[5],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='L23ExcOverview6.png').plot()


            #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=numpy.pi/2,st_contrast=100)   
            #KremkowOverviewFigure(dsv,ParameterSet({'neuron' : l4_exc,'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,12)},plot_file_name='ExcOverview.png').plot()
            
            TrialToTrialVariabilityComparison(data_store,ParameterSet({'sheet_name1' : 'V1_Exc_L4','sheet_name2' : 'V1_Exc_L2/3','data_dg' : 0.93 , 'data_ni' : 1.19}),fig_param={'dpi' : 200,'figsize': (10,5)}).plot()
            
            if True:            
                dsv = param_filter_query(data_store,st_name='NaturalImageWithEyeMovement')            
                KremkowOverviewFigure(dsv,ParameterSet({'neuron' : l4_exc,'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='NMOverview.png').plot()

                SNRAnalysis(data_store,ParameterSet({"neuron" : analog_ids[0]}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='SNR1.png').plot()                        
                SNRAnalysis(data_store,ParameterSet({"neuron" : analog_ids[1]}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='SNR2.png').plot()                        
                SNRAnalysis(data_store,ParameterSet({"neuron" : analog_ids[2]}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='SNR3.png').plot()                        

                StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids[0],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparison1.png').plot()
                StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids[1],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparison2.png').plot()
                StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids[2],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparison3.png').plot()
                StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids[3],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparison4.png').plot()
                if l23_flag:
                    StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids23[0],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparison1L23.png').plot()
                    StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids23[1],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparison2L23.png').plot()
                    StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids23[2],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparison3L23.png').plot()
                    StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids23[3],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (10,12)},plot_file_name='StimulusResponseComparison4L23.png').plot()

                #dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L4')    
                #FanoFactor_Baudot_et_al(dsv,ParameterSet({}),plot_file_name='FanoFactorL4.png').plot()
                
                #dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L2/3')    
                #FanoFactor_Baudot_et_al(dsv,ParameterSet({}),plot_file_name='FanoFactorL23.png').plot()
               
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating')    
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedExcRaster.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedInhRaster.png').plot({'SpikeRasterPlot.group_trials':True})

            #OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0On.png").plot()        
            #OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0Off.png").plot()
    
            
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=[0,numpy.pi/2])    
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            
            
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc1.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc2.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc3.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[3], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc4.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[4], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc5.png").plot({'Vm_plot.y_lim' : (-90,-50)})

            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh1.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh2.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh3.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[3], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh4.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[4], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh5.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            
            if l23_flag:
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL231.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL232.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL233.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[3], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL234.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[4], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL235.png").plot({'Vm_plot.y_lim' : (-80,-50)})

                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="InhL231.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="InhL232.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="InhL233.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[3], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="InhL234.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[4], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="InhL235.png").plot({'Vm_plot.y_lim' : (-80,-50)})


            
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0On.png").plot()
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0Off.png").plot()
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids())[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN1On.png").plot()
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_esyn_ids())[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN1Off.png").plot()
            
            # tuninc curves
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate','Analog_F0andF1'])    
            PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[:6]), 'sheet_name' : 'V1_Exc_L4','centered'  : True,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='TCsExcL4.png',fig_param={'dpi' : 100,'figsize': (32,15)}).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})
            PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh[:6]), 'sheet_name' : 'V1_Inh_L4','centered'  : True,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='TCsInhL4.png',fig_param={'dpi' : 100,'figsize': (32,15)}).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})
            if l23_flag:
                PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids23[:6]), 'sheet_name' : 'V1_Exc_L2/3','centered'  : True,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='TCsExcL23.png',fig_param={'dpi' : 100,'figsize': (32,15)}).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})
                PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh23[:6]), 'sheet_name' : 'V1_Inh_L2/3','centered'  : True,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='TCsInhL23.png',fig_param={'dpi' : 100,'figsize': (32,15)}).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})
            
            
            MRfig(param_filter_query(data_store,sheet_name=exc_sheets,st_contrast=[100],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'V1_Exc_L2/3'}),plot_file_name='MR.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()

            SpontActOverview(data_store,ParameterSet({'l4_exc_neuron' : analog_ids[0], 'l4_inh_neuron' : analog_ids_inh[0],'l23_exc_neuron' : analog_ids23[0],'l23_inh_neuron' : analog_ids_inh23[0]}),plot_file_name='SpontActOverview.png', fig_param={'dpi' : 200,'figsize': (15,12)}).plot()
            SpontStatisticsOverview(data_store,ParameterSet({}), fig_param={'dpi' : 200,'figsize': (15,10)},plot_file_name='SpontStatisticsOverview.png').plot()
            OrientationTuningSummaryFiringRates(data_store,ParameterSet({'exc_sheet_name1': 'V1_Exc_L4','inh_sheet_name1': 'V1_Inh_L4','exc_sheet_name2': 'V1_Exc_L2/3','inh_sheet_name2': 'V1_Inh_L2/3'}),fig_param={'dpi' : 200,'figsize': (22.5,13.5)},plot_file_name='OrientationTuningSummary.png').plot({'*.fontsize' : 19})            
            OrientationTuningSummaryAnalogSignals(data_store,ParameterSet({'exc_sheet_name1': 'V1_Exc_L4','inh_sheet_name1': 'V1_Inh_L4','exc_sheet_name2': 'V1_Exc_L2/3','inh_sheet_name2': 'V1_Inh_L2/3'}),fig_param={'dpi' : 200,'figsize': (22.5,12)},plot_file_name='OrientationTuningSummaryAnalogSignals.png').plot({'*.fontsize' : 19})            
            
            TrialCrossCorrelationAnalysis(data_store,ParameterSet({'neurons' : list(analog_ids),'window_length' : 83,'sheet_name' : 'V1_Exc_L4'}),fig_param={"dpi" : 100,"figsize": (25,12)},plot_file_name="trial-to-trial-cross-correlation.png").plot()
            if l23_flag:
                TrialCrossCorrelationAnalysis(data_store,ParameterSet({'neurons' : list(analog_ids23),'window_length' : 83,'sheet_name' : 'V1_Exc_L2/3'}),fig_param={"dpi" : 100,"figsize": (25,12)},plot_file_name="trial-to-trial-cross-correlation23.png").plot()

            
            dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L4')    
            MeanVsVarainceOfVM(dsv,ParameterSet({'neurons' : list(analog_ids[:5])}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='TrialToTrialMeanVsVarianceOfVM.png').plot()

            if l23_flag:
                dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L2/3')    
                MeanVsVarainceOfVM(dsv,ParameterSet({'neurons' : list(analog_ids23[:5])}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='TrialToTrialMeanVsVarianceOfVM.png').plot()


            dsv = param_filter_query(data_store,st_name='NaturalImageWithEyeMovement')    
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}, 'spontaneous' : True}),plot_file_name='NMExc.png',fig_param={'dpi' : 100,'figsize': (28,12)}).plot({'Vm_plot.y_lim' : (-70,-50),'Conductance_plot.y_lim' : (0,50.0)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}, 'spontaneous' : True}),plot_file_name='NMInh.png',fig_param={'dpi' : 100,'figsize': (28,12)}).plot({'Vm_plot.y_lim' : (-70,-50),'Conductance_plot.y_lim' : (0,50.0)})

            dsv = param_filter_query(data_store,st_name=['InternalStimulus'])    
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {},'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhAnalog.png').plot()
    
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

            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_inh_or[0].get_value_by_id(l4_inh),numpy.pi)  for o in orr])],st_contrast=100)   
            KremkowOverviewFigure(dsv,ParameterSet({'neuron' : l4_inh,'sheet_name' : 'V1_Inh_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='ExcOverviewInh.png').plot()

            
            # orientation tuning plotting
            #dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],value_name='LGNAfferentOrientation')   
            #PerNeuronValuePlot(dsv,ParameterSet({"cortical_view" : True}),plot_file_name='ORSet.png').plot()
            
            #dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],value_name='orientation preference of Firing rate',analysis_algorithm='PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage',st_contrast=100)    
            #PerNeuronValuePlot(dsv,ParameterSet({"cortical_view" : True}),plot_file_name='ORComputed.png').plot()
            
def perform_analysis_and_visualization_spont(data_store):
    analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()
    spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()

    analog_ids23 = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh23 = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_esyn_ids()
    spike_ids23 = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh23 = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_spike_train_ids()

    dsv = param_filter_query(data_store,st_direct_stimulation_name="None",st_name="InternalStimulus")   
    TrialAveragedFiringRate(dsv,ParameterSet({})).analyse()
    Analog_MeanSTDAndFanoFactor(param_filter_query(dsv,st_direct_stimulation_name="None"),ParameterSet({})).analyse()
    TrialAveragedVarianceAndVarianceRatioOfConductances(param_filter_query(dsv,st_direct_stimulation_name="None"),ParameterSet({})).analyse()
    PSTH(param_filter_query(dsv,st_direct_stimulation_name="None"),ParameterSet({'bin_length' : 5.0})).analyse()
    Irregularity(param_filter_query(data_store,st_direct_stimulation_name="None"),ParameterSet({})).analyse()
    CrossCorrelationOfExcitatoryAndInhibitoryConductances(param_filter_query(data_store,st_direct_stimulation_name="None"),ParameterSet({})).analyse()
    PopulationMean(data_store,ParameterSet({})).analyse()

    dsv = param_filter_query(data_store,value_name=['Mean(ECond)','Mean(ICond)','Mean(VM)','Mean(ECond)/Mean(ICond)'],sheet_name=["V1_Exc_L4","V1_Exc_L2/3"])   
    PerNeuronValuePlot(dsv,ParameterSet({"cortical_view" : False}),plot_file_name='ExcInhMeasure_hist.png',fig_param={'dpi' : 100,'figsize': (16,8)}).plot({'*.title' : None,'HistogramPlot.plot[0,1].x_label' : 'Inh Cond (micro siemens)','HistogramPlot.plot[1,1].x_label' : 'Exc Cond (micro siemens)','HistogramPlot.plot[2,1].x_label' : 'Exc/Inh Cond','HistogramPlot.plot[3,1].x_label' : 'Membrane potential (mV)'})
    
    dsv = param_filter_query(data_store,value_name=['Firing rate'],sheet_name=["V1_Exc_L4","V1_Exc_L2/3"])
    PerNeuronValuePlot(dsv,ParameterSet({"cortical_view" : False}),plot_file_name='Spont_firing_rate.png',fig_param={'dpi' : 100,'figsize': (4,8)}).plot({'*.title' : None,'HistogramPlot.plot[0,1].x_scale' : 'log','HistogramPlot.plot[0,1].log' : True,'HistogramPlot.plot[0,0].x_scale' : 'log','HistogramPlot.plot[0,0].log' : True})

    dsv = param_filter_query(data_store,value_name=['Firing rate'],sheet_name=["V1_Exc_L4","V1_Exc_L2/3"])
    PerNeuronValuePlot(dsv,ParameterSet({"cortical_view" : False}),plot_file_name='Spont_firing_rate_loglog.png',fig_param={'dpi' : 100,'figsize': (4,8)}).plot({'*.title' : None,'HistogramPlot.plot[0,1].x_scale' : 'log','HistogramPlot.plot[0,1].log' : True,'HistogramPlot.plot[0,0].x_scale' : 'log','HistogramPlot.plot[0,0].log' : True,'HistogramPlot.plot[0,1].y_scale' : 'log','HistogramPlot.plot[0,0].y_scale' : 'log','*.y_label' : None})

