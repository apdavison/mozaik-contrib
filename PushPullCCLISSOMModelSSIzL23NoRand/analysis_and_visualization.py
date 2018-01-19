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

import psutil
import os
process = psutil.Process(os.getpid())


def memory_usage_psutil():
    # return the memory usage in MB
    return process.memory_percent()
    

def analysis(data_store,analog_ids,analog_ids_inh,analog_ids23=None,analog_ids_inh23=None):
        sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3']))
        exc_sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4','V1_Exc_L2/3']))
        l23_flag = ('V1_Exc_L2/3' in set(sheets)) and (analog_ids23 != None)
        
        logger.info('0: ' + str(memory_usage_psutil()))
        
        TrialAveragedFiringRate(param_filter_query(data_store,sheet_name=sheets,st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({})).analyse()
        TrialAveragedFiringRate(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()
        logger.info('1: ' + str(memory_usage_psutil()))
        Irregularity(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()
        
        
        PSTH(param_filter_query(data_store),ParameterSet({'bin_length' : 10.0 })).analyse()
        
        logger.info('2: ' + str(memory_usage_psutil()))
        #SpikeCount(param_filter_query(data_store,sheet_name=exc_sheets),ParameterSet({'bin_length' : 13.0 })).analyse()
        NeuronToNeuronAnalogSignalCorrelations(param_filter_query(data_store,analysis_algorithm='PSTH'),ParameterSet({'convert_nan_to_zero' : True})).analyse()
        
        logger.info('3: ' + str(memory_usage_psutil()))
        
        PopulationMeanAndVar(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()
                
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
        OperationPNVfromPNVS(pnv,lambda x,y: -(x+y),'-(x+y)',dsv,ParameterSet({})).analyse()
        
        pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Inh_L4',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(ECond)',st_direct_stimulation_name='None').get_analysis_result()[0]
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Inh_L4',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Exc_Cond')
        SubtractPNVfromPNVS(pnv,dsv,ParameterSet({})).analyse()

        pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Inh_L4',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(ICond)',st_direct_stimulation_name='None').get_analysis_result()[0]
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Inh_L4',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Inh_Cond')
        SubtractPNVfromPNVS(pnv,dsv,ParameterSet({})).analyse()

        pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Inh_L4',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(VM)',st_direct_stimulation_name='None').get_analysis_result()[0]
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Inh_L4',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Vm')
        OperationPNVfromPNVS(pnv,lambda x,y: -(x+y),'-(x+y)',dsv,ParameterSet({})).analyse()



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
            OperationPNVfromPNVS(pnv,lambda x,y: -(x+y),'-(x+y)',dsv,ParameterSet({})).analyse()
            
            pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Inh_L2/3',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(ECond)',st_direct_stimulation_name='None').get_analysis_result()[0]
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Inh_L2/3',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Exc_Cond')
            SubtractPNVfromPNVS(pnv,dsv,ParameterSet({})).analyse()

            pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Inh_L2/3',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(ICond)',st_direct_stimulation_name='None').get_analysis_result()[0]
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Inh_L2/3',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Inh_Cond')
            SubtractPNVfromPNVS(pnv,dsv,ParameterSet({})).analyse()

            pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Inh_L2/3',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(VM)',st_direct_stimulation_name='None').get_analysis_result()[0]
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Inh_L2/3',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Vm')
            OperationPNVfromPNVS(pnv,lambda x,y: -(x+y),'-(x+y)',dsv,ParameterSet({})).analyse()
            
            
        
        logger.info('9: ' + str(memory_usage_psutil()))
        dsv = queries.param_filter_query(data_store,y_axis_name='spike count (bin=13.0)')   
        mozaik.analysis.analysis.TrialToTrialFanoFactorOfAnalogSignal(dsv,ParameterSet({})).analyse()
        
        logger.info('10: ' + str(memory_usage_psutil()))
        
        TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(data_store,sheet_name='V1_Exc_L4',st_name="NaturalImageWithEyeMovement",analysis_algorithm='ActionPotentialRemoval'),ParameterSet({'neurons' : list(analog_ids)})).analyse()
        TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(data_store,sheet_name='V1_Exc_L4',st_name="NaturalImageWithEyeMovement",analysis_algorithm='PSTH'),ParameterSet({'neurons' : list(analog_ids)})).analyse()
        TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(data_store,sheet_name='V1_Exc_L4',st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='ActionPotentialRemoval',st_contrast=100),ParameterSet({'neurons' : list(analog_ids)})).analyse()
        TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(data_store,sheet_name='V1_Exc_L4',st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='PSTH',st_contrast=100),ParameterSet({'neurons' : list(analog_ids)})).analyse()

        
        logger.info('11: ' + str(memory_usage_psutil()))
        if l23_flag:
            TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(data_store,sheet_name='V1_Exc_L2/3',st_name="NaturalImageWithEyeMovement",analysis_algorithm='ActionPotentialRemoval'),ParameterSet({'neurons' : list(analog_ids23)})).analyse()
            TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(data_store,sheet_name='V1_Exc_L2/3',st_name="NaturalImageWithEyeMovement",analysis_algorithm='PSTH'),ParameterSet({'neurons' : list(analog_ids23)})).analyse()
            TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(data_store,sheet_name='V1_Exc_L2/3',st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='ActionPotentialRemoval',st_contrast=100),ParameterSet({'neurons' : list(analog_ids23)})).analyse()
            TrialToTrialCrossCorrelationOfAnalogSignalList(param_filter_query(data_store,sheet_name='V1_Exc_L2/3',st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='PSTH',st_contrast=100),ParameterSet({'neurons' : list(analog_ids23)})).analyse()
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
    l23 = 'V1_Exc_L2/3' in set(data_store.sheets()) 

    analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()
    spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()

    if l23:
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


    l4_exc_or_many = numpy.array(spike_ids)[numpy.nonzero(numpy.array([circular_dist(l4_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids]) < 0.1)[0]]
    l4_exc_or_many_analog = numpy.array(analog_ids)[numpy.nonzero(numpy.array([circular_dist(l4_exc_or.get_value_by_id(i),0,numpy.pi)  for i in analog_ids]) < 0.1)[0]]
    l4_exc_or_many_analog_inh = numpy.array(analog_ids_inh)[numpy.nonzero(numpy.array([circular_dist(l4_inh_or.get_value_by_id(i),0,numpy.pi)  for i in analog_ids_inh]) < 0.1)[0]]
    if l23:
        l23_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L2/3')[0]
        l23_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L2/3')[0]
        l23_exc_or_many = numpy.array(spike_ids23)[numpy.nonzero(numpy.array([circular_dist(l23_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids23]) < 0.1)[0]]
        l23_exc_or_many_analog = numpy.array(analog_ids23)[numpy.nonzero(numpy.array([circular_dist(l23_exc_or.get_value_by_id(i),0,numpy.pi)  for i in analog_ids23]) < 0.1)[0]]
        l23_exc_or_many_analog_inh = numpy.array(analog_ids_inh23)[numpy.nonzero(numpy.array([circular_dist(l23_inh_or.get_value_by_id(i),0,numpy.pi)  for i in analog_ids_inh23]) < 0.1)[0]]


    if False:
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name=['V1_Exc_L2/3','V1_Exc_L4'])   
        Analog_F0andF1(dsv,ParameterSet({})).analyse()
        
        pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Exc_L4',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(VM)',st_direct_stimulation_name='None').get_analysis_result()[0]
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Exc_L4',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Vm')
        OperationPNVfromPNVS(pnv,lambda x,y: -(x+y),'-(x+y)',dsv,ParameterSet({})).analyse()

        pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Exc_L2/3',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(VM)',st_direct_stimulation_name='None').get_analysis_result()[0]
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Exc_L2/3',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Vm')
        OperationPNVfromPNVS(pnv,lambda x,y: -(x+y),'-(x+y)',dsv,ParameterSet({})).analyse()

        if l23:
            MRfig(param_filter_query(data_store,sheet_name=['V1_Exc_L2/3','V1_Exc_L4'],st_contrast=[100],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'V1_Exc_L2/3'}),plot_file_name='MR.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()
            MRfigReal(param_filter_query(data_store,sheet_name=['V1_Exc_L2/3','V1_Exc_L4'],st_contrast=[100],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'V1_Exc_L2/3'}),plot_file_name='MRReal.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()
        else:
            MRfig(param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'None'}),plot_file_name='MR.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()
        0/0

        
    if True:
        if l23:
           analysis(data_store,l4_exc_or_many_analog,l4_exc_or_many_analog_inh,l23_exc_or_many_analog,l23_exc_or_many_analog_inh)
        else:
           analysis(data_store,l4_exc_or_many_analog,l4_exc_or_many_analog_inh)


    if True:
            # self sustained plotting
            dsv = param_filter_query(data_store,st_name=['InternalStimulus'])    
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})

            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalog.png').plot()
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhAnalog.png').plot()    
            if l23:
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSL23ExcAnalog.png').plot()
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSL23InhAnalog.png').plot()    
                RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neurons' : spike_ids23,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcRasterL23.png').plot({'SpikeRasterPlot.group_trials':True})
                RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neurons' : spike_ids_inh23,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhRasterL23.png').plot({'SpikeRasterPlot.group_trials':True})

            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating')
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc_or_many_analog[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc1.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc_or_many_analog[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc2.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc_or_many_analog[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc3.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc_or_many_analog[3], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc4.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc_or_many_analog[4], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc5.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc_or_many_analog[5], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc6.png").plot({'Vm_plot.y_lim' : (-80,-50)})

    	    if l23:
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc_or_many_analog[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL231.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc_or_many_analog[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL232.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc_or_many_analog[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL233.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc_or_many_analog[3], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL234.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc_or_many_analog[4], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL235.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc_or_many_analog[5], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL236.png").plot({'Vm_plot.y_lim' : (-80,-50)})


            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh1.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh2.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            if l23:
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="InhL231.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="InhL232.png").plot({'Vm_plot.y_lim' : (-80,-50)})


            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedExcRaster.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedInhRaster.png').plot({'SpikeRasterPlot.group_trials':True})

            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : list(analog_ids),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedExcRasterWithAnalog.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : list(analog_ids_inh),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedInhRasterWithAnalog.png').plot({'SpikeRasterPlot.group_trials':True})
            if l23:
                RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neurons' : list(analog_ids23),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedExcL23RasterWithAnalog.png').plot({'SpikeRasterPlot.group_trials':True})
                RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neurons' : list(analog_ids_inh23),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedInhL23RasterWithAnalog.png').plot({'SpikeRasterPlot.group_trials':True})



            #dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L4',st_contrast=100)    
            #FanoFactor_Baudot_et_al(dsv,ParameterSet({}),plot_file_name='FanoFactorL4.png').plot()
            
            #dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L2/3',st_contrast=100)    
            #FanoFactor_Baudot_et_al(dsv,ParameterSet({}),plot_file_name='FanoFactorL23.png').plot()
            if True:
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
                
                if l23:
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


            StimulusResponseComparison(data_store,ParameterSet({'neuron' : l4_exc_or_many_analog[0],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,5)},plot_file_name='StimulusResponseComparison1.png').plot()
            StimulusResponseComparison(data_store,ParameterSet({'neuron' : l4_exc_or_many_analog[1],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,5)},plot_file_name='StimulusResponseComparison2.png').plot()
            StimulusResponseComparison(data_store,ParameterSet({'neuron' : l4_exc_or_many_analog[2],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,5)},plot_file_name='StimulusResponseComparison3.png').plot()
            #StimulusResponseComparison(data_store,ParameterSet({'neuron' : l4_exc_or_many_analog[3],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (19,5)},plot_file_name='StimulusResponseComparison4.png').plot()
            if l23:
                StimulusResponseComparison(data_store,ParameterSet({'neuron' : l23_exc_or_many_analog[0],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (19,5)},plot_file_name='StimulusResponseComparisonL231.png').plot()
                StimulusResponseComparison(data_store,ParameterSet({'neuron' : l23_exc_or_many_analog[1],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (19,5)},plot_file_name='StimulusResponseComparisonL232.png').plot()
                StimulusResponseComparison(data_store,ParameterSet({'neuron' : l23_exc_or_many_analog[2],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (19,5)},plot_file_name='StimulusResponseComparisonL233.png').plot()
                #StimulusResponseComparison(data_store,ParameterSet({'neuron' : l23_exc_or_many_analog[3],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (19,5)},plot_file_name='StimulusResponseComparisonL234.png').plot()

            
            if l23:
                TrialCrossCorrelationAnalysis(data_store,ParameterSet({'neurons' : list(l23_exc_or_many_analog),'window_length' : 83, 'sheet_name' : 'V1_Exc_L2/3'}),fig_param={"dpi" : 100,"figsize": (25,12)},plot_file_name="trial-to-trial-cross-correlation23.png").plot()
            TrialCrossCorrelationAnalysis(data_store,ParameterSet({'neurons' : list(l4_exc_or_many_analog),'window_length' : 83, 'sheet_name' : 'V1_Exc_L4'}),fig_param={"dpi" : 100,"figsize": (25,12)},plot_file_name="trial-to-trial-cross-correlation4.png").plot()

            
            dsv = param_filter_query(data_store,value_name=['Mean(ECond)','Mean(ICond)','Mean(VM)','Mean(ECond)/Mean(ICond)','Firing rate'])   
            PerNeuronValuePlot(dsv,ParameterSet({"cortical_view" : False}),plot_file_name='Histograms.png',fig_param={'dpi' : 100,'figsize': (25,12)}).plot({'*.title' : None,'HistogramPlot.plot[3,1].x_scale' : 'log','HistogramPlot.plot[3,1].log' : True,'HistogramPlot.plot[3,0].x_scale' : 'log','HistogramPlot.plot[3,0].log' : True})
            
            dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L4')    
            MeanVsVarainceOfVM(dsv,ParameterSet({'neurons' : list(l4_exc_or_many_analog[:5])}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='TrialToTrialMeanVsVarianceOfVML4.png').plot()
            if l23:
                dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L2/3')    
                MeanVsVarainceOfVM(dsv,ParameterSet({'neurons' : list(l23_exc_or_many_analog[:5])}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='TrialToTrialMeanVsVarianceOfVML23.png').plot()
            
                
            #SNRAnalysis(data_store,ParameterSet({"neuron" : l4_exc}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='SNR1.png').plot()                        
            if l23:      
                dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'],value_name='LGNAfferentOrientation')   
                PerNeuronValuePlot(dsv,ParameterSet({'cortical_view' : True}),plot_file_name='Map.png').plot({'*.colorbar' : False,'*.x_axis' : None, '*.y_axis' : None})


def perform_analysis_and_visualization_stc(data_store):
    l23 = 'V1_Exc_L2/3' in set(data_store.sheets()) 
    analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()
    spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()

    NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()


    if l23:
        analog_ids23 = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_esyn_ids()
        analog_ids_inh23 = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_esyn_ids()
        spike_ids23 = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_spike_train_ids()
        spike_ids_inh23 = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_spike_train_ids()
        l23_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L2/3')[0]
        l23_exc_or_many = numpy.array(spike_ids23)[numpy.nonzero(numpy.array([circular_dist(l23_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids23]) < 0.1)[0]]    
        idx23 = data_store.get_sheet_indexes(sheet_name='V1_Exc_L2/3',neuron_ids=l23_exc_or_many)

    l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')[0]
    l4_exc_or_many = numpy.array(spike_ids)[numpy.nonzero(numpy.array([circular_dist(l4_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids]) < 0.1)[0]]
    idx4 = data_store.get_sheet_indexes(sheet_name='V1_Exc_L4',neuron_ids=l4_exc_or_many)

    x = data_store.get_neuron_postions()['V1_Exc_L4'][0][idx4]
    y = data_store.get_neuron_postions()['V1_Exc_L4'][1][idx4]
    center4 = l4_exc_or_many[numpy.nonzero(numpy.sqrt(numpy.multiply(x,x)+numpy.multiply(y,y)) < 0.2)[0]]
    analog_center4 = set(center4).intersection(analog_ids)
    logger.info(str(analog_center4))
    
    if l23:
        x = data_store.get_neuron_postions()['V1_Exc_L2/3'][0][idx23]
        y = data_store.get_neuron_postions()['V1_Exc_L2/3'][1][idx23]
        center23 = l23_exc_or_many[numpy.nonzero(numpy.sqrt(numpy.multiply(x,x)+numpy.multiply(y,y)) < 0.2)[0]]
        analog_center23 = set(center23).intersection(analog_ids23)
        logger.info(str(analog_center23))    
    
    
    TrialAveragedFiringRate(param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Exc_L2/3'],st_name='DriftingSinusoidalGratingDisk'),ParameterSet({})).analyse()
    
    dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Exc_L2/3'],st_name='DriftingSinusoidalGratingDisk')
    Analog_F0andF1(dsv,ParameterSet({})).analyse()

    dsv = param_filter_query(data_store,st_name='InternalStimulus',st_direct_stimulation_name="None")
    Analog_MeanSTDAndFanoFactor(dsv,ParameterSet({})).analyse()

    pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Exc_L4',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(VM)',st_direct_stimulation_name='None').get_analysis_result()[0]
    dsv = param_filter_query(data_store,st_name='DriftingSinusoidalGratingDisk',sheet_name='V1_Exc_L4',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Vm')
    OperationPNVfromPNVS(pnv,lambda x,y: -(x+y),'-(x+y)',dsv,ParameterSet({})).analyse()

    if l23:
        pnv = param_filter_query(data_store,st_name='InternalStimulus',sheet_name='V1_Exc_L2/3',analysis_algorithm=['Analog_MeanSTDAndFanoFactor'],value_name='Mean(VM)',st_direct_stimulation_name='None').get_analysis_result()[0]
        dsv = param_filter_query(data_store,st_name='DriftingSinusoidalGratingDisk',sheet_name='V1_Exc_L2/3',analysis_algorithm=['Analog_F0andF1'],value_name='F0_Vm')
        OperationPNVfromPNVS(pnv,lambda x,y: -(x+y),'-(x+y)',dsv,ParameterSet({})).analyse()
    
    
    dsv = param_filter_query(data_store,st_name='DriftingSinusoidalGratingDisk',analysis_algorithm='TrialAveragedFiringRate')    
    SizeTuningAnalysis(dsv,ParameterSet({'neurons' : center4.tolist(), 'sheet_name' : 'V1_Exc_L4'})).analyse()
    if l23:
        SizeTuningAnalysis(dsv,ParameterSet({'neurons' : center23.tolist(), 'sheet_name' : 'V1_Exc_L2/3'})).analyse()        
    data_store.save()

    dsv = param_filter_query(data_store,st_name='DriftingSinusoidalGratingDisk',analysis_algorithm=['TrialAveragedFiringRate'])    
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'radius', 'neurons': list(center4), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='SizeTuningExcL4.png',fig_param={'dpi' : 100,'figsize': (32,7)}).plot()
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'radius', 'neurons': list(center4), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : True, 'polar' : False, 'pool'  : False}),plot_file_name='SizeTuningExcL4M.png',fig_param={'dpi' : 100,'figsize': (32,7)}).plot()
    
    if l23:
        PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'radius', 'neurons': list(center23), 'sheet_name' : 'V1_Exc_L2/3','centered'  : False,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='SizeTuningExcL23.png',fig_param={'dpi' : 100,'figsize': (32,7)}).plot()        
        PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'radius', 'neurons': list(center23), 'sheet_name' : 'V1_Exc_L2/3','centered'  : False,'mean' : True, 'polar' : False, 'pool'  : False}),plot_file_name='SizeTuningExcL23M.png',fig_param={'dpi' : 100,'figsize': (32,7)}).plot()        
    
    if l23:
        SizeTuningOverview(data_store,ParameterSet({'l4_neurons' : list(center4),'l23_neurons' : list(center23),'l4_neurons_analog' : list(analog_center4),'l23_neurons_analog' : list(analog_center23)}),plot_file_name='SizeTuningOverview.png',fig_param={'dpi' : 300,'figsize': (18,8)}).plot()
    else:
        SizeTuningOverview(data_store,ParameterSet({'l4_neurons' : list(center4),'l23_neurons' : [],'l4_neurons_analog' : list(analog_center4),'l23_neurons_analog' : []}),plot_file_name='SizeTuningOverview.png',fig_param={'dpi' : 300,'figsize': (18,8)}).plot()
    data_store.save()
    
    if True:
        dsv = param_filter_query(data_store,st_name=['DriftingSinusoidalGratingDisk'])    
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : list(analog_center4)[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_ExcL4_1.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : list(analog_center4)[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_ExcL4_2.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : list(analog_center4)[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_ExcL4_3.png').plot()
    
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : list(analog_center23)[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_ExcL23_1.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : list(analog_center23)[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_ExcL23_2.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : list(analog_center23)[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_ExcL23_3.png').plot()
    


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


    l4_exc_or_many = numpy.array(spike_ids)[numpy.nonzero(numpy.array([circular_dist(l4_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids]) < 0.3)[0]]
    l23_exc_or_many = numpy.array(spike_ids23)[numpy.nonzero(numpy.array([circular_dist(l23_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids23]) < 0.3)[0]]
    l4_inh_or_many = numpy.array(spike_ids_inh)[numpy.nonzero(numpy.array([circular_dist(l4_inh_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids_inh]) < 0.3)[0]]
    l23_inh_or_many = numpy.array(spike_ids_inh23)[numpy.nonzero(numpy.array([circular_dist(l23_inh_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids_inh23]) < 0.3)[0]]


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
    logger.info(center23)
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'surround_orientation', 'neurons': list(center4), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : True, 'polar' : False, 'pool'  : False}),plot_file_name='OCTCExcL4M.png',fig_param={'dpi' : 100,'figsize': (32,15)}).plot()
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'surround_orientation', 'neurons': list(center23), 'sheet_name' : 'V1_Exc_L2/3','centered'  : False,'mean' : True, 'polar' : False, 'pool'  : False}),plot_file_name='OCTCExcL23M.png',fig_param={'dpi' : 100,'figsize': (32,15)}).plot()        
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'surround_orientation', 'neurons': list(center4), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='OCTCExcL4.png',fig_param={'dpi' : 100,'figsize': (32,15)}).plot()
    PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'surround_orientation', 'neurons': list(center23), 'sheet_name' : 'V1_Exc_L2/3','centered'  : False,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='OCTCExcL23.png',fig_param={'dpi' : 100,'figsize': (32,15)}).plot()        
   
    x = data_store.get_neuron_postions()['V1_Inh_L4'][0][idx_inh4]
    y = data_store.get_neuron_postions()['V1_Inh_L4'][1][idx_inh4]
    center4_inh = l4_inh_or_many[numpy.nonzero(numpy.sqrt(numpy.multiply(x,x)+numpy.multiply(y,y)) < 0.1)[0]]
    x = data_store.get_neuron_postions()['V1_Inh_L2/3'][0][idx_inh23]
    y = data_store.get_neuron_postions()['V1_Inh_L2/3'][1][idx_inh23]
    center23_inh = l23_inh_or_many[numpy.nonzero(numpy.sqrt(numpy.multiply(x,x)+numpy.multiply(y,y)) < 0.1)[0]]
    
    analog_center4_inh = set(center4_inh).intersection(analog_ids_inh)
    analog_center23_inh = set(center23_inh).intersection(analog_ids_inh23)

    
    if True:
        dsv = param_filter_query(data_store,st_name=['DriftingSinusoidalGratingCenterSurroundStimulus'],st_surround_orientation=[0,numpy.pi/2])    
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : list(analog_center23)[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_ExcL23_1.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : list(analog_center23)[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_ExcL23_2.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : list(analog_center23)[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_ExcL23_3.png').plot()

        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : list(analog_center23_inh)[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_InhL23_1.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : list(analog_center23_inh)[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_InhL23_2.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : list(analog_center23_inh)[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Overview_InhL23_3.png').plot()
        


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

    if l23_flag:
        l23_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L2/3')[0]
        l23_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L2/3')[0]    
        
    
    l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')
    l4_exc_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Exc_L4')
    l4_exc = analog_ids[numpy.argmin([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_exc_or[0].get_value_by_id(analog_ids),l4_exc_phase[0].get_value_by_id(analog_ids))])]
    l4_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L4')
    l4_inh_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Inh_L4')
    l4_inh = analog_ids_inh[numpy.argmin([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_inh_or[0].get_value_by_id(analog_ids_inh),l4_inh_phase[0].get_value_by_id(analog_ids_inh))])]
    l4_exc_or_many = numpy.array(l4_exc_or[0].ids)[numpy.nonzero(numpy.array([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_exc_or[0].values,l4_exc_phase[0].values)]) < 0.1)[0]]
    
    l4_exc_or_many = list(set(l4_exc_or_many) &  set(spike_ids))
    
    if l23_flag:
        l23_exc_or_many = numpy.array(l23_exc_or.ids)[numpy.nonzero(numpy.array([circular_dist(o,0,numpy.pi)  for o in l23_exc_or.values]) < 0.1)[0]]
        l23_exc_or_many = list(set(l23_exc_or_many) &  set(spike_ids23))    
    
    
    
    
    
    orr = list(set([MozaikParametrized.idd(s).orientation for s in queries.param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_contrast=100).get_stimuli()]))                
        
    l4_exc_or_many_analog = numpy.array(analog_ids)[numpy.nonzero(numpy.array([circular_dist(l4_exc_or[0].get_value_by_id(i),0,numpy.pi)  for i in analog_ids]) < 0.1)[0]]
    l4_inh_or_many_analog = numpy.array(analog_ids_inh)[numpy.nonzero(numpy.array([circular_dist(l4_inh_or[0].get_value_by_id(i),0,numpy.pi)  for i in analog_ids_inh]) < 0.15)[0]]
    
    if l23_flag:
        l23_inh_or_many_analog = numpy.array(analog_ids_inh23)[numpy.nonzero(numpy.array([circular_dist(l23_inh_or.get_value_by_id(i),0,numpy.pi)  for i in analog_ids_inh23]) < 0.15)[0]]
        l23_exc_or_many_analog = numpy.array(analog_ids23)[numpy.nonzero(numpy.array([circular_dist(l23_exc_or.get_value_by_id(i),0,numpy.pi)  for i in analog_ids23]) < 0.1)[0]]


    #if l23_flag:
    #    SpontActOverview(data_store,ParameterSet({'l4_exc_neuron' : analog_ids[0], 'l4_inh_neuron' : analog_ids_inh[0],'l23_exc_neuron' : analog_ids23[0],'l23_inh_neuron' : analog_ids_inh23[0]}),plot_file_name='SpontActOverview.png', fig_param={'dpi' : 200,'figsize': (18,14.5)}).plot()
    #else:
    #    SpontActOverview(data_store,ParameterSet({'l4_exc_neuron' : analog_ids[0], 'l4_inh_neuron' : analog_ids_inh[0],'l23_exc_neuron' : -1,'l23_inh_neuron' : -1}),plot_file_name='SpontActOverview.png', fig_param={'dpi' : 200,'figsize': (18,14.5)}).plot()
    #0/0

    #if l23_flag:
    #    MRfig(param_filter_query(data_store,sheet_name=['V1_Exc_L2/3','V1_Exc_L4'],st_contrast=[100],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'V1_Exc_L2/3'}),plot_file_name='MR.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()
    #else:
    #    MRfig(param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'None'}),plot_file_name='MR.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()
    
    #Kremkow_plots.OrientationTuningSummary(data_store,ParameterSet({'exc_sheet_name': 'V1_Exc_L4','inh_sheet_name': 'V1_Inh_L4'}),fig_param={'dpi' : 100,'figsize': (15,9)},plot_file_name='OrientationTuningSummaryL4.png').plot()            
    #0/0
    
    if True:
       if l23_flag:
        analysis(data_store,analog_ids,analog_ids_inh,analog_ids23,analog_ids_inh23)
       else:
        analysis(data_store,analog_ids,analog_ids_inh)

        
    if False:
            pref_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
            ort_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
            pref_fr_50 = 0#param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[5],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
            ort_fr_50 = 0#param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[5],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
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

    
    if False:
        if l23_flag:
            LSV1MReponseOverview(data_store,ParameterSet({'l4_exc_neuron' : l4_exc_or_many_analog[1], 'l4_inh_neuron' : l4_inh_or_many_analog[1],'l23_exc_neuron' : l23_exc_or_many_analog[0], 'l23_inh_neuron' : l23_inh_or_many_analog[0]}),fig_param={'dpi' : 70,'figsize': (18,9)},plot_file_name='SingleCellOverview.png').plot({'*.x_ticks':[0,2.0], "*.title" : None,'*.Vm_plot.y_lim' : (-80,-50),'*.Conductance_plot.y_lim' : (0,70)})
        
        dsv = param_filter_query(data_store,st_name='NaturalImageWithEyeMovement')
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[16], 'sheet_activity' : {}, 'spontaneous' : False}), fig_param={'dpi' : 100,'figsize': (18,10)},plot_file_name="NatExcL4.png").plot({'Vm_plot.y_lim' : (-80,-50),'*.x_ticks':[0,2.0]})
        if l23_flag:
    	    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[15], 'sheet_activity' : {}, 'spontaneous' : False}), fig_param={'dpi' : 100,'figsize': (18,10)},plot_file_name="NatExcL23.png").plot({'Vm_plot.y_lim' : (-80,-50),'*.x_ticks':[0,2.0]})

        #SpontRasterOverview(data_store,ParameterSet({}),fig_param={'dpi' : 70,'figsize': (20,20)},plot_file_name='SpontRasterOverview.png').plot()
        #SpontAnalogOverview(data_store,ParameterSet({}),fig_param={'dpi' : 70,'figsize': (20,20)},plot_file_name='SpontAnalogOverview.png').plot()
      
    
    if True: # PLOTTING
            activity_plot_param =    {
                   'frame_rate' : 5,  
                   'bin_width' : 5.0, 
                   'scatter' :  True,
                   'resolution' : 0
            }   
            
            #Kremkow_plots.OrientationTuningSummary(data_store,ParameterSet({'exc_sheet_name': 'V1_Exc_L4','inh_sheet_name': 'V1_Inh_L4'}),fig_param={'dpi' : 100,'figsize': (15,9)},plot_file_name='OrientationTuningSummaryL4.png').plot()            
            #if l23_flag:
            #    Kremkow_plots.OrientationTuningSummary(data_store,ParameterSet({'exc_sheet_name': 'V1_Exc_L2/3','inh_sheet_name': 'V1_Inh_L2/3'}),fig_param={'dpi' : 100,'figsize': (15,9)},plot_file_name='OrientationTuningSummaryL23.png').plot()            
            
            
            
            # self sustained plotting
            dsv = param_filter_query(data_store,st_name=['InternalStimulus'],st_direct_stimulation_name="None")    
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalog.png').plot()
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhAnalog.png').plot()    
            
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})
            if l23_flag:
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalog23.png').plot()
	        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhAnalog23.png').plot()    

                RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neurons' : spike_ids23,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcRasterL23.png').plot({'SpikeRasterPlot.group_trials':True})
                RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neurons' : spike_ids_inh23,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhRasterL23.png').plot({'SpikeRasterPlot.group_trials':True})
            
            if False:
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

            if l23_flag:
                        TrialToTrialVariabilityComparison(data_store,ParameterSet({'sheet_name1' : 'V1_Exc_L4','sheet_name2' : 'V1_Exc_L2/3','data_dg' : 0.93 , 'data_ni' : 1.19}),fig_param={'dpi' : 200,'figsize': (15,7.5)},plot_file_name='TrialToTrialVariabilityComparison.png').plot()
            else:
                        TrialToTrialVariabilityComparison(data_store,ParameterSet({'sheet_name1' : 'V1_Exc_L4','sheet_name2' : 'V1_Exc_L4','data_dg' : 0.93 , 'data_ni' : 1.19}),fig_param={'dpi' : 200,'figsize': (15,7.5)},plot_file_name='TrialToTrialVariabilityComparison.png').plot()
            
            if False:            
                dsv = param_filter_query(data_store,st_name='NaturalImageWithEyeMovement')            
                KremkowOverviewFigure(dsv,ParameterSet({'neuron' : l4_exc,'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='NMOverview.png').plot()

                SNRAnalysis(data_store,ParameterSet({"neuron" : analog_ids[0]}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='SNR1.png').plot()                        
                SNRAnalysis(data_store,ParameterSet({"neuron" : analog_ids[1]}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='SNR2.png').plot()                        
                SNRAnalysis(data_store,ParameterSet({"neuron" : analog_ids[2]}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name='SNR3.png').plot()                        

                StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids[0],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (18,12)},plot_file_name='StimulusResponseComparison1.png').plot()
                StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids[1],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (18,12)},plot_file_name='StimulusResponseComparison2.png').plot()
                StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids[2],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (18,12)},plot_file_name='StimulusResponseComparison3.png').plot()
                StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids[3],'sheet_name' : 'V1_Exc_L4'}),fig_param={'dpi' : 100,'figsize': (18,12)},plot_file_name='StimulusResponseComparison4.png').plot()
                if l23_flag:
                    StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids23[0],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (18,12)},plot_file_name='StimulusResponseComparison1L23.png').plot()
                    StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids23[1],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (18,12)},plot_file_name='StimulusResponseComparison2L23.png').plot()
                    StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids23[2],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (18,12)},plot_file_name='StimulusResponseComparison3L23.png').plot()
                    StimulusResponseComparison(data_store,ParameterSet({'neuron' : analog_ids23[3],'sheet_name' : 'V1_Exc_L2/3'}),fig_param={'dpi' : 100,'figsize': (18,12)},plot_file_name='StimulusResponseComparison4L23.png').plot()

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
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            
            
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc1.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc2.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc3.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[3], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc4.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[4], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc5.png").plot({'Vm_plot.y_lim' : (-90,-50)})

            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh1.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh2.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh3.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[3], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh4.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[4], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh5.png").plot({'Vm_plot.y_lim' : (-90,-50)})
            
            if l23_flag:
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL231.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL232.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL233.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[3], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL234.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[4], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL235.png").plot({'Vm_plot.y_lim' : (-80,-50)})

                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="InhL231.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="InhL232.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="InhL233.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[3], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="InhL234.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[4], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="InhL235.png").plot({'Vm_plot.y_lim' : (-80,-50)})


            
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0On.png").plot()
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_esyn_ids())[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0Off.png").plot()
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids())[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN1On.png").plot()
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : sorted(param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_esyn_ids())[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN1Off.png").plot()
            
            
            
            

            SpontStatisticsOverview(data_store,ParameterSet({}), fig_param={'dpi' : 200,'figsize': (18,12)},plot_file_name='SpontStatisticsOverview.png').plot()
            if l23_flag:
                MRfig(param_filter_query(data_store,sheet_name=['V1_Exc_L2/3','V1_Exc_L4'],st_contrast=[100],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'V1_Exc_L2/3'}),plot_file_name='MR.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()
                MRfigReal(param_filter_query(data_store,sheet_name=['V1_Exc_L2/3','V1_Exc_L4'],st_contrast=[100],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'V1_Exc_L2/3'}),plot_file_name='MRReal.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()
            else:
                MRfig(param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'None'}),plot_file_name='MR.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()

            
            if l23_flag:
                SpontActOverview(data_store,ParameterSet({'l4_exc_neuron' : analog_ids[0], 'l4_inh_neuron' : analog_ids_inh[0],'l23_exc_neuron' : analog_ids23[0],'l23_inh_neuron' : analog_ids_inh23[0]}),plot_file_name='SpontActOverview.png', fig_param={'dpi' : 200,'figsize': (18,14.5)}).plot()
                OrientationTuningSummaryAnalogSignals(data_store,ParameterSet({'exc_sheet_name1': 'V1_Exc_L4','inh_sheet_name1': 'V1_Inh_L4','exc_sheet_name2': 'V1_Exc_L2/3','inh_sheet_name2': 'V1_Inh_L2/3'}),fig_param={'dpi' : 200,'figsize': (18,12)},plot_file_name='OrientationTuningSummaryAnalogSignals.png').plot({'*.fontsize' : 19,'*.y_lim' : (0,None)})            
                OrientationTuningSummaryFiringRates(data_store,ParameterSet({'exc_sheet_name1': 'V1_Exc_L4','inh_sheet_name1': 'V1_Inh_L4','exc_sheet_name2': 'V1_Exc_L2/3','inh_sheet_name2': 'V1_Inh_L2/3'}),fig_param={'dpi' : 200,'figsize': (18,12)},plot_file_name='OrientationTuningSummary.png').plot({'*.fontsize' : 19})            
            else:
                SpontActOverview(data_store,ParameterSet({'l4_exc_neuron' : analog_ids[0], 'l4_inh_neuron' : analog_ids_inh[0],'l23_exc_neuron' : -1,'l23_inh_neuron' : -1}),plot_file_name='SpontActOverview.png', fig_param={'dpi' : 200,'figsize': (18,14.5)}).plot()
                Kremkow_plots.OrientationTuningSummary(data_store,ParameterSet({'exc_sheet_name': 'V1_Exc_L4','inh_sheet_name': 'V1_Inh_L4'}),fig_param={'dpi' : 100,'figsize': (15,9)},plot_file_name='OrientationTuningSummaryL4.png').plot()            
                OrientationTuningSummaryAnalogSignals(data_store,ParameterSet({'exc_sheet_name1': 'V1_Exc_L4','inh_sheet_name1': 'V1_Inh_L4','exc_sheet_name2': 'None','inh_sheet_name2': 'None'}),fig_param={'dpi' : 200,'figsize': (18,12)},plot_file_name='OrientationTuningSummaryAnalogSignals.png').plot({'*.fontsize' : 19,'*.y_lim' : (0,None)})            
            
            dsv = param_filter_query(data_store,st_name='NaturalImageWithEyeMovement')    
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}, 'spontaneous' : True}),plot_file_name='NMExc.png',fig_param={'dpi' : 100,'figsize': (28,12)}).plot({'Vm_plot.y_lim' : (-70,-50),'Conductance_plot.y_lim' : (0,50.0)})
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}, 'spontaneous' : True}),plot_file_name='NMInh.png',fig_param={'dpi' : 100,'figsize': (28,12)}).plot({'Vm_plot.y_lim' : (-70,-50),'Conductance_plot.y_lim' : (0,50.0)})
                        
            TrialCrossCorrelationAnalysis(data_store,ParameterSet({'neurons1' : list(analog_ids),'sheet_name1' : 'V1_Exc_L4','neurons2' : list(analog_ids23),'sheet_name2' : 'V1_Exc_L2/3', 'window_length' : 250}),fig_param={"dpi" : 100,"figsize": (15,6.5)},plot_file_name="trial-to-trial-cross-correlation.png").plot({'*.Vm.title' : None, '*.fontsize' : 19})
            
            dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L4')    
            MeanVsVarainceOfVM(dsv,ParameterSet({'neurons' : list(analog_ids[:5])}),fig_param={'dpi' : 100,'figsize': (15,7.5)},plot_file_name='TrialToTrialMeanVsVarianceOfVM.png').plot()

            if l23_flag:
                dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L2/3')    
                MeanVsVarainceOfVM(dsv,ParameterSet({'neurons' : list(analog_ids23[:5])}),fig_param={'dpi' : 100,'figsize': (15,7.5)},plot_file_name='TrialToTrialMeanVsVarianceOfVM.png').plot()

            
            dsv = queries.param_filter_query(data_store,value_name=['orientation HWHH of Firing rate','orientation CV(Firing rate)'],sheet_name=["V1_Exc_L2/3"],st_contrast=100)
            PerNeuronValueScatterPlot(dsv, ParameterSet({'only_matching_units' : False,'ignore_nan' : True}),plot_file_name='CVvsHWHH.png').plot({'*.x_lim' : (0,90),'*.y_lim' : (0,1.0)})

            # tuninc curves
            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate','Analog_F0andF1'])    
            PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[:6]), 'sheet_name' : 'V1_Exc_L4','centered'  : True,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='TCsExcL4.png',fig_param={'dpi' : 100,'figsize': (15,7.5)}).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})
            PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh[:6]), 'sheet_name' : 'V1_Inh_L4','centered'  : True,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='TCsInhL4.png',fig_param={'dpi' : 100,'figsize': (15,7.5)}).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})
            if l23_flag:
                PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids23[:6]), 'sheet_name' : 'V1_Exc_L2/3','centered'  : True,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='TCsExcL23.png',fig_param={'dpi' : 100,'figsize': (15,7.5)}).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})
                PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh23[:6]), 'sheet_name' : 'V1_Inh_L2/3','centered'  : True,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='TCsInhL23.png',fig_param={'dpi' : 100,'figsize': (15,7.5)}).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})
            
            
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
    sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3']))
    l23_flag = 'V1_Exc_L2/3' in set(sheets)

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
        analog_ids23 = []
        analog_ids_inh23 = []
        spike_ids23 = []
        spike_ids_inh23 = []


    dsv = param_filter_query(data_store,st_direct_stimulation_name="None",st_name="InternalStimulus")   
    TrialAveragedFiringRate(dsv,ParameterSet({})).analyse()
    Analog_MeanSTDAndFanoFactor(param_filter_query(dsv,st_direct_stimulation_name="None"),ParameterSet({})).analyse()
    TrialAveragedVarianceAndVarianceRatioOfConductances(param_filter_query(dsv,st_direct_stimulation_name="None"),ParameterSet({})).analyse()
    PSTH(param_filter_query(dsv,st_direct_stimulation_name="None"),ParameterSet({'bin_length' : 10.0})).analyse()
    Irregularity(param_filter_query(data_store,st_direct_stimulation_name="None"),ParameterSet({})).analyse()
    CrossCorrelationOfExcitatoryAndInhibitoryConductances(param_filter_query(data_store,st_direct_stimulation_name="None"),ParameterSet({})).analyse()
    NeuronToNeuronAnalogSignalCorrelations(param_filter_query(data_store,analysis_algorithm='PSTH'),ParameterSet({'convert_nan_to_zero' : True})).analyse()
    PopulationMeanAndVar(data_store,ParameterSet({})).analyse()
    data_store.save()
    
    dsv = param_filter_query(data_store,st_name=['InternalStimulus'],st_direct_stimulation_name="None")    
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalogL4.png').plot()
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhAnalogL4.png').plot()    
            
    RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})
    RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})

    if l23_flag:
        dsv = param_filter_query(data_store,st_name=['InternalStimulus'],st_direct_stimulation_name="None")    
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalogL23.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhAnalogL23.png').plot()    
                
        RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neurons' : spike_ids23,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcRasterL23.png').plot({'SpikeRasterPlot.group_trials':True})
        RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neurons' : spike_ids_inh23,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhRasterL23.png').plot({'SpikeRasterPlot.group_trials':True})
    
    #SpontStatisticsOverview(data_store,ParameterSet({}), fig_param={'dpi' : 200,'figsize': (15,10)},plot_file_name='SpontStatisticsOverview.png').plot()
    
    if l23_flag:
        SpontActOverview(data_store,ParameterSet({'l4_exc_neuron' : analog_ids[0], 'l4_inh_neuron' : analog_ids_inh[0],'l23_exc_neuron' : analog_ids23[0],'l23_inh_neuron' : analog_ids_inh23[0]}),plot_file_name='SpontActOverview.png', fig_param={'dpi' : 200,'figsize': (15,12)}).plot()
        SpontStatisticsOverview(data_store,ParameterSet({}), fig_param={'dpi' : 200,'figsize': (15,10)},plot_file_name='SpontStatisticsOverview.png').plot()
    else:
        SpontActOverview(data_store,ParameterSet({'l4_exc_neuron' : analog_ids[0], 'l4_inh_neuron' : analog_ids_inh[0],'l23_exc_neuron' : -1,'l23_inh_neuron' : -1}),plot_file_name='SpontActOverview.png', fig_param={'dpi' : 200,'figsize': (15,12)}).plot()
    


def perform_analysis_and_visualization_or_small(data_store):
    l23 = 'V1_Exc_L2/3' in set(data_store.sheets()) 

    analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()
    spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()

    if l23:
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


    l4_exc_or_many = numpy.array(spike_ids)[numpy.nonzero(numpy.array([circular_dist(l4_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids]) < 0.1)[0]]
    l4_exc_or_many_analog = numpy.array(analog_ids)[numpy.nonzero(numpy.array([circular_dist(l4_exc_or.get_value_by_id(i),0,numpy.pi)  for i in analog_ids]) < 0.1)[0]]
    l4_exc_or_many_analog_inh = numpy.array(analog_ids_inh)[numpy.nonzero(numpy.array([circular_dist(l4_inh_or.get_value_by_id(i),0,numpy.pi)  for i in analog_ids_inh]) < 0.1)[0]]
    if l23:
        l23_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L2/3')[0]
        l23_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L2/3')[0]
        l23_exc_or_many = numpy.array(spike_ids23)[numpy.nonzero(numpy.array([circular_dist(l23_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids23]) < 0.1)[0]]
        l23_exc_or_many_analog = numpy.array(analog_ids23)[numpy.nonzero(numpy.array([circular_dist(l23_exc_or.get_value_by_id(i),0,numpy.pi)  for i in analog_ids23]) < 0.1)[0]]
        l23_exc_or_many_analog_inh = numpy.array(analog_ids_inh23)[numpy.nonzero(numpy.array([circular_dist(l23_inh_or.get_value_by_id(i),0,numpy.pi)  for i in analog_ids_inh23]) < 0.1)[0]]

    if True:
        if l23:
           analysis(data_store,l4_exc_or_many_analog,l4_exc_or_many_analog_inh,l23_exc_or_many_analog,l23_exc_or_many_analog_inh)
        else:
           analysis(data_store,l4_exc_or_many_analog,l4_exc_or_many_analog_inh)


    if True:
            # self sustained plotting
            dsv = param_filter_query(data_store,st_name=['InternalStimulus'])    
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})

            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalog.png').plot()
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhAnalog.png').plot()    
            if l23:
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : analog_ids23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSL23ExcAnalog.png').plot()
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSL23InhAnalog.png').plot()    
                RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neurons' : spike_ids23,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcRasterL23.png').plot({'SpikeRasterPlot.group_trials':True})
                RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neurons' : spike_ids_inh23,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhRasterL23.png').plot({'SpikeRasterPlot.group_trials':True})

            dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating')
            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc_or_many_analog[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Exc1.png").plot({'Vm_plot.y_lim' : (-80,-50)})

    	    if l23:
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc_or_many_analog[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL231.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc_or_many_analog[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL232.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc_or_many_analog[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="ExcL233.png").plot({'Vm_plot.y_lim' : (-80,-50)})


            OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh1.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="Inh2.png").plot({'Vm_plot.y_lim' : (-80,-50)})
            if l23:
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="InhL231.png").plot({'Vm_plot.y_lim' : (-80,-50)})
                OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : analog_ids_inh23[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (25,12)},plot_file_name="InhL232.png").plot({'Vm_plot.y_lim' : (-80,-50)})


            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedExcRaster.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedInhRaster.png').plot({'SpikeRasterPlot.group_trials':True})

            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : list(analog_ids),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedExcRasterWithAnalog.png').plot({'SpikeRasterPlot.group_trials':True})
            RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : list(analog_ids_inh),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedInhRasterWithAnalog.png').plot({'SpikeRasterPlot.group_trials':True})
            if l23:
                RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neurons' : list(analog_ids23),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedExcL23RasterWithAnalog.png').plot({'SpikeRasterPlot.group_trials':True})
                RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neurons' : list(analog_ids_inh23),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedInhL23RasterWithAnalog.png').plot({'SpikeRasterPlot.group_trials':True})



            #dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L4',st_contrast=100)    
            #FanoFactor_Baudot_et_al(dsv,ParameterSet({}),plot_file_name='FanoFactorL4.png').plot()
            
            #dsv = param_filter_query(data_store,sheet_name = 'V1_Exc_L2/3',st_contrast=100)    
            #FanoFactor_Baudot_et_al(dsv,ParameterSet({}),plot_file_name='FanoFactorL23.png').plot()
            if True:
                pref_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
                ort_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[numpy.pi/2],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)
                spont = param_filter_query(data_store,st_name='InternalStimulus',sheet_name=['V1_Exc_L4'],analysis_algorithm=['TrialAveragedFiringRate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc_or_many)

                pylab.figure()
                pylab.bar([1,2,3,4,5],[numpy.mean(pref_fr_100),numpy.mean(ort_fr_100),0,0,numpy.mean(spont)])
                pylab.xticks([1.4,2.4,3.4,4.4,5.4],['PREF100','ORT100','PREF50','ORT50','SPONT'])
                pylab.xlim([0.8,6.0])
                pylab.ylabel("Firing rate")

                pylab.savefig(Global.root_directory+"Orientation_responseL4.png")
                
                if l23:
                    pref_fr_100 = param_filter_query(data_store,sheet_name=['V1_Exc_L2/3'],st_contrast=[100],analysis_algorithm=['TrialAveragedFiringRate'],st_orientation=[0],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)
                    spont = param_filter_query(data_store,st_name='InternalStimulus',sheet_name=['V1_Exc_L2/3'],analysis_algorithm=['TrialAveragedFiringRate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(l23_exc_or_many)

                    pylab.figure()
                    pylab.bar([1,2,3,4,5],[numpy.mean(pref_fr_100),0,0,0,numpy.mean(spont)])
                    pylab.xticks([1.4,2.4,3.4,4.4,5.4],['PREF100','ORT100','PREF50','ORT50','SPONT'])
                    pylab.xlim([0.8,6.0])
                    pylab.ylabel("Firing rate")
                
                    pylab.savefig(Global.root_directory+"Orientation_responseL23.png")


            if l23:
                SpontActOverview(data_store,ParameterSet({'l4_exc_neuron' : analog_ids[0], 'l4_inh_neuron' : analog_ids_inh[0],'l23_exc_neuron' : analog_ids23[0],'l23_inh_neuron' : analog_ids_inh23[0]}),plot_file_name='SpontActOverview.png', fig_param={'dpi' : 200,'figsize': (15,12)}).plot()
            else:
                SpontActOverview(data_store,ParameterSet({'l4_exc_neuron' : analog_ids[0], 'l4_inh_neuron' : analog_ids_inh[0],'l23_exc_neuron' : -1,'l23_inh_neuron' : -1}),plot_file_name='SpontActOverview.png', fig_param={'dpi' : 200,'figsize': (15,12)}).plot()
            
            SpontStatisticsOverview(data_store,ParameterSet({}), fig_param={'dpi' : 200,'figsize': (15,10)},plot_file_name='SpontStatisticsOverview.png').plot()
            
            if l23:
                MRfig(param_filter_query(data_store,sheet_name=['V1_Exc_L2/3','V1_Exc_L4'],st_contrast=[100],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'V1_Exc_L2/3'}),plot_file_name='MR.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()
            else:
                MRfig(param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'None'}),plot_file_name='MR.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()

            if l23:      
                dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'],value_name='LGNAfferentOrientation')   
                PerNeuronValuePlot(dsv,ParameterSet({'cortical_view' : True}),plot_file_name='Map.png').plot({'*.colorbar' : False,'*.x_axis' : None, '*.y_axis' : None})
            if l23:      
                TrialToTrialVariabilityComparison(data_store,ParameterSet({'sheet_name1' : 'V1_Exc_L4','sheet_name2' : 'V1_Exc_L2/3','data_dg' : 0.93 , 'data_ni' : 1.19}),fig_param={'dpi' : 200,'figsize': (10,5)},plot_file_name='TrialToTrialVariabilityComparison.png').plot()
            else:
                TrialToTrialVariabilityComparison(data_store,ParameterSet({'sheet_name1' : 'V1_Exc_L4','sheet_name2' : 'None','data_dg' : 0.93 , 'data_ni' : 1.19}),fig_param={'dpi' : 200,'figsize': (10,5)},plot_file_name='TrialToTrialVariabilityComparison.png').plot()
