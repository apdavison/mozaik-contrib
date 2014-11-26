import numpy
import mozaik
import pylab
from mozaik.visualization.plotting import *
from mozaik.analysis.technical import NeuronAnnotationsToPerNeuronValues
from mozaik.analysis.analysis import *
from mozaik.analysis.vision import *
from mozaik.storage.queries import *
from mozaik.storage.datastore import PickledDataStore
from mozaik.tools.circ_stat import circular_dist



class STC(Plotting):
    required_parameters = ParameterSet({
        'l4_exc_neurons': list,
        'l4_inh_neurons': list,
        'l23_exc_neurons': list,
        'l23_inh_neurons': list,
    })

    def subplot(self, subplotspec):
        gs = gridspec.GridSpecFromSubplotSpec(12, 18, subplot_spec=subplotspec,
                                              hspace=1.0, wspace=1.0)
        
        return {
                'Layer4Exc' : (PlotTuningCurve(self.datastore,ParameterSet({'parameter_name' : 'radius', 'neurons': self.parameters.l4_exc_neurons, 'sheet_name' : 'V1_Exc_L4'})),gs[0:3,:],{}),
                'Layer4Inh' : (PlotTuningCurve(self.datastore,ParameterSet({'parameter_name' : 'radius', 'neurons': self.parameters.l4_inh_neurons, 'sheet_name' : 'V1_Inh_L4'})),gs[3:6,:],{}),
                'Layer23Exc' : (PlotTuningCurve(self.datastore,ParameterSet({'parameter_name' : 'radius', 'neurons': self.parameters.l23_exc_neurons, 'sheet_name' : 'V1_Exc_L2/3'})),gs[6:9,:],{}),
                'Layer23Inh' : (PlotTuningCurve(self.datastore,ParameterSet({'parameter_name' : 'radius', 'neurons': self.parameters.l23_inh_neurons, 'sheet_name' : 'V1_Inh_L2/3'})),gs[9:12,:],{}),
        }

class CTC(Plotting):
    required_parameters = ParameterSet({
        'l4_exc_neurons': list,
        'l4_inh_neurons': list,
        'l23_exc_neurons': list,
        'l23_inh_neurons': list,
    })

    def subplot(self, subplotspec):
        gs = gridspec.GridSpecFromSubplotSpec(12, 18, subplot_spec=subplotspec,
                                              hspace=1.0, wspace=1.0)
        
        return {
                'Layer4Exc' : (PlotTuningCurve(self.datastore,ParameterSet({'parameter_name' : 'contrast', 'neurons': self.parameters.l4_exc_neurons, 'sheet_name' : 'V1_Exc_L4'})),gs[0:3,:],{}),
                'Layer4Inh' : (PlotTuningCurve(self.datastore,ParameterSet({'parameter_name' : 'contrast', 'neurons': self.parameters.l4_inh_neurons, 'sheet_name' : 'V1_Inh_L4'})),gs[3:6,:],{}),
                'Layer23Exc' : (PlotTuningCurve(self.datastore,ParameterSet({'parameter_name' : 'contrast', 'neurons': self.parameters.l23_exc_neurons, 'sheet_name' : 'V1_Exc_L2/3'})),gs[6:9,:],{}),
                'Layer23Inh' : (PlotTuningCurve(self.datastore,ParameterSet({'parameter_name' : 'contrast', 'neurons': self.parameters.l23_inh_neurons, 'sheet_name' : 'V1_Inh_L2/3'})),gs[9:12,:],{}),
        }


class OCTC(Plotting):
    required_parameters = ParameterSet({
        'l4_exc_neurons': list,
        'l4_inh_neurons': list,
        'l23_exc_neurons': list,
        'l23_inh_neurons': list,
    })

    def subplot(self, subplotspec):
        gs = gridspec.GridSpecFromSubplotSpec(12, 18, subplot_spec=subplotspec,
                                              hspace=1.0, wspace=1.0)
        
        return {
                'Layer4Exc' : (PlotTuningCurve(self.datastore,ParameterSet({'parameter_name' : 'surround_orientation', 'neurons': self.parameters.l4_exc_neurons, 'sheet_name' : 'V1_Exc_L4'})),gs[0:3,:],{}),
                'Layer4Inh' : (PlotTuningCurve(self.datastore,ParameterSet({'parameter_name' : 'surround_orientation', 'neurons': self.parameters.l4_inh_neurons, 'sheet_name' : 'V1_Inh_L4'})),gs[3:6,:],{}),
                'Layer23Exc' : (PlotTuningCurve(self.datastore,ParameterSet({'parameter_name' : 'surround_orientation', 'neurons': self.parameters.l23_exc_neurons, 'sheet_name' : 'V1_Exc_L2/3'})),gs[6:9,:],{}),
                'Layer23Inh' : (PlotTuningCurve(self.datastore,ParameterSet({'parameter_name' : 'surround_orientation', 'neurons': self.parameters.l23_inh_neurons, 'sheet_name' : 'V1_Inh_L2/3'})),gs[9:12,:],{}),
        }


class OR(Plotting):
    required_parameters = ParameterSet({
        'l4_exc_neurons': list,
        'l4_inh_neurons': list,
        'l23_exc_neurons': list,
        'l23_inh_neurons': list,
    })

    def subplot(self, subplotspec):
        gs = gridspec.GridSpecFromSubplotSpec(12, 18, subplot_spec=subplotspec,
                                              hspace=1.0, wspace=1.0)
        
        return {
                'Layer4Exc' : (PlotTuningCurve(self.datastore,ParameterSet({'parameter_name' : 'orientation', 'neurons': self.parameters.l4_exc_neurons, 'sheet_name' : 'V1_Exc_L4'})),gs[0:3,:],{}),
                'Layer4Inh' : (PlotTuningCurve(self.datastore,ParameterSet({'parameter_name' : 'orientation', 'neurons': self.parameters.l4_inh_neurons, 'sheet_name' : 'V1_Inh_L4'})),gs[3:6,:],{}),
                'Layer23Exc' : (PlotTuningCurve(self.datastore,ParameterSet({'parameter_name' : 'orientation', 'neurons': self.parameters.l23_exc_neurons, 'sheet_name' : 'V1_Exc_L2/3'})),gs[6:9,:],{}),
                'Layer23Inh' : (PlotTuningCurve(self.datastore,ParameterSet({'parameter_name' : 'orientation', 'neurons': self.parameters.l23_inh_neurons, 'sheet_name' : 'V1_Inh_L2/3'})),gs[9:12,:],{}),
        }


class TuningComparison(Plotting):

    def subplot(self, subplotspec):
        gs = gridspec.GridSpecFromSubplotSpec(2,4, subplot_spec=subplotspec)
        plots = {}    
        dsv = param_filter_query(self.datastore,sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'])   
        for i,sheet in enumerate(dsv.sheets()):
            dsv1 = param_filter_query(dsv,sheet_name=sheet,value_name='LGNAfferentOrientation')   
            plots[sheet + 'A'] = (PerNeuronValuePlot(dsv1,ParameterSet({})),gs[0,i],{})

            dsv1 = param_filter_query(dsv,sheet_name=sheet,value_name='orientation preference',analysis_algorithm='PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage',st_contrast=100)    
            plots[sheet + 'B'] = (PerNeuronValuePlot(dsv1,ParameterSet({})),gs[1,i],{'title' : None})
        return plots    

class HWHH(Plotting):

    def subplot(self, subplotspec):
        gs = gridspec.GridSpecFromSubplotSpec(16,1, subplot_spec=subplotspec)
        plots = {}    
        dsv = param_filter_query(self.datastore,sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'])   
        print dsv.sheets()
        print ['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3']
        for i,sheet in enumerate(['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3']):
            dsv1 = param_filter_query(data_store,value_name=['orientation HWHH'],sheet_name=sheet)    
            plots[sheet] = (PerNeuronValueScatterPlot(dsv1,ParameterSet({})),gs[i*4:i*4+3,0],{})
           
        return plots    





def perform_analysis_and_visualization(data_store):
    pref_or = numpy.pi/2
    #find neuron with preference closet to pref_or
    l4_analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
    l4_analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()
    l23_analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_esyn_ids()
    l23_analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_esyn_ids()
    NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
    l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')
    l4_exc = l4_analog_ids[numpy.argmin([circular_dist(o,pref_or,numpy.pi)  for o in l4_exc_or[0].get_value_by_id(l4_analog_ids)])]
    l4_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L4')
    l4_inh = l4_analog_ids_inh[numpy.argmin([circular_dist(o,pref_or,numpy.pi)  for o in l4_inh_or[0].get_value_by_id(l4_analog_ids_inh)])]
    l23_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L2/3')
    l23_exc = l23_analog_ids[numpy.argmin([circular_dist(o,pref_or,numpy.pi)  for o in l23_exc_or[0].get_value_by_id(l23_analog_ids)])]
    l23_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L2/3')
    l23_inh = l23_analog_ids_inh[numpy.argmin([circular_dist(o,pref_or,numpy.pi)  for o in l23_inh_or[0].get_value_by_id(l23_analog_ids_inh)])]
    #l4_exc_or_many = numpy.array(l4_exc_or[0].ids)[numpy.nonzero(numpy.array([circular_dist(o,pref_or,numpy.pi)  for (o,p) in zip(l4_exc_or[0].values,l4_exc_phase[0].values)]) < 0.1)[0]]

    # Find neurons for which spikes were recorded and whose orientation is close to pref_or
    l4_spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
    l4_spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()
    l23_spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_spike_train_ids()
    l23_spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L2/3").get_segments()[0].get_stored_spike_train_ids()

    l4_spike_exc_or_close_to_pi_half = numpy.array(l4_spike_ids)[numpy.nonzero(numpy.array([circular_dist(o,pref_or,numpy.pi)  for o in l4_exc_or[0].get_value_by_id(l4_spike_ids)]) < 0.1)[0]].tolist()
    l4_spike_inh_or_close_to_pi_half = numpy.array(l4_spike_ids_inh)[numpy.nonzero(numpy.array([circular_dist(o,pref_or,numpy.pi)  for o in l4_inh_or[0].get_value_by_id(l4_spike_ids_inh)]) < 0.1)[0]].tolist()
    l23_spike_exc_or_close_to_pi_half = numpy.array(l23_spike_ids)[numpy.nonzero(numpy.array([circular_dist(o,pref_or,numpy.pi)  for o in l23_exc_or[0].get_value_by_id(l23_spike_ids)]) < 0.1)[0]].tolist()
    l23_spike_inh_or_close_to_pi_half = numpy.array(l23_spike_ids_inh)[numpy.nonzero(numpy.array([circular_dist(o,pref_or,numpy.pi)  for o in l23_inh_or[0].get_value_by_id(l23_spike_ids_inh)]) < 0.1)[0]].tolist()


    print "Prefered orientation of plotted exc neurons:"
    print 'id: ' + str(l4_exc)
    print "Prefered orientation of plotted inh neurons:"
    print 'id: ' + str(l4_inh)


    if True:  #ANALYSIS
        #TrialAveragedFiringRate(data_store,ParameterSet({'stimulus_type':"DriftingSinusoidalGratingCenterSurroundStimulus"})).analyse()
        TrialAveragedFiringRate(data_store,ParameterSet({})).analyse()
        #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate')    
        #GaussianTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
        #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'])   
        #Analog_F0andF1(dsv,ParameterSet({})).analyse()
        #TrialVariability(dsv,ParameterSet({'vm': True,  'cond_exc': False, 'cond_inh': False})).analyse()
        #TrialMean(dsv,ParameterSet({'vm': True,  'cond_exc': False, 'cond_inh': False})).analyse()
        
        #GSTA(param_filter_query(data_store,sheet_name='V1_Exc_L4'),ParameterSet({'neurons' : [l4_exc], 'length' : 250.0 }),tags=['GSTA']).analyse()
        #Precision(param_filter_query(data_store,sheet_name='V1_Exc_L4'),ParameterSet({'neurons' : [l4_exc], 'bin_length' : 10.0 })).analyse()
        #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'])    
        #PSTH(dsv,ParameterSet({'bin_length' : 10.0})).analyse()
       
        #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'],y_axis_name='psth (bin=10.0)',st_contrast=100,analysis_algorithm='PSTH')
        #dsv1 = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'],value_name='orientation preference',st_contrast=100,analysis_algorithm='PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage')
        #ModulationRatio(dsv+dsv1,ParameterSet({})).analyse()
        
        #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'])  
        #PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
        #PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage(dsv,ParameterSet({'parameter_name' : 'phase'})).analyse()
        
        #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name=['V1_Exc_L4'],value_name='orientation preference',st_contrast=100,analysis_algorithm='PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage')
        #LocalHomogeneityIndex(dsv,ParameterSet({'sigma' : 0.1})).analyse()
        #dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],value_name='LGNAfferentOrientation')
        #LocalHomogeneityIndex(dsv,ParameterSet({'sigma' : 100.0})).analyse()
        
        data_store.save()
        



    if True: # PLOTTING
        #verify_connectivity(data_store)
        
        activity_plot_param =    {
               'frame_rate' : 5,  
               'bin_width' : 20.0, 
               'scatter' :  True,
               'resolution' : 0
        }       
        
        #dsv = param_filter_query(data_store,st_name='NaturalImageWithEyeMovement')    
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : l23_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        
        #dsv = param_filter_query(data_store,st_name='DriftingGratingWithEyeMovement')    
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : l23_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        
        
        #dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],value_name='LocalHomogeneityIndex' + '(0.1:orientation preference)')   
        #PerNeuronValuePlot(dsv,ParameterSet({})).plot({'*.colorbar' : False,'*.x_axis' : None, '*.y_axis' : None})

        #dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4'],value_name='LocalHomogeneityIndex' + '(100.0:LGNAfferentOrientation)')   
        #PerNeuronValuePlot(dsv,ParameterSet({})).plot({'*.colorbar' : False,'*.x_axis' : None, '*.y_axis' : None})
        
        
        #dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'],value_name='LGNAfferentOrientation')   
        #PerNeuronValuePlot(dsv,ParameterSet({})).plot({'*.colorbar' : False,'*.x_axis' : None, '*.y_axis' : None})
        #dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'],value_name='orientation preference',analysis_algorithm='PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage',st_contrast=100)    
        #PerNeuronValuePlot(dsv,ParameterSet({})).plot({'*.colorbar' : False, '*.x_axis' : None, '*.y_axis' : None})
        
        #TuningComparison(data_store,ParameterSet({}),fig_param={'dpi' : 100,'figsize': (17,7.5)},plot_file_name="/home/antolikjan/Doc/Talks/JCUnicMar2013/TuningComparison.png").plot({'*.colorbar' : False,'*.x_axis' : None, '*.y_axis' : None})
            
       
        #dsv = param_filter_query(data_store,sheet_name=['V1_Exc_L4','V1_Inh_L4'],analysis_algorithm='TrialVariability',st_orientation = 0)    
        #PerNeuronValueScatterPlot(dsv,ParameterSet({})).plot({'ScatterPlot.identity_line' : True, 'ScatterPlot.mark_means' : True})
        
        #dsv = param_filter_query(data_store,value_name=['orientation HWHH'],sheet_name=['V1_Exc_L4','V1_Inh_L4','V1_Exc_L2/3','V1_Inh_L2/3'])    
        #PerNeuronValueScatterPlot(dsv,ParameterSet({}),plot_file_name="/home/antolikjan/Doc/Talks/JCUnicMar2013/MR.png").plot({ 'ScatterPlot.x_lim' : (0,90), 'ScatterPlot.y_lim' : (0,90), 'ScatterPlot.identity_line' : True})
        
        #HWHH(data_store,ParameterSet({}),fig_param={'dpi' : 100,'figsize': (8,12)},plot_file_name="/home/antolikjan/Doc/Talks/JCUnicMar2013/HWHH.png").plot({ '*.x_lim' : (0,90), '*.y_lim' : (0,90), '*.identity_line' : True, '*.title' : None})
        
        # LHI orientation selectivity correlation
        #dsv1 = param_filter_query(data_store,value_name='orientation selectivity',sheet_name=['V1_Exc_L4'],analysis_algorithm='PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage',st_contrast=100)    
        #dsv2 = param_filter_query(data_store,value_name='LocalHomogeneityIndex' + '(0.1:orientation preference)',sheet_name=['V1_Exc_L4'])
        #PerNeuronValueScatterPlot(dsv1+dsv2,ParameterSet({}),plot_file_name="/home/antolikjan/Doc/Talks/JCUnicMar2013/LHIvsSelectivity.png").plot()
        
        
        #dsv = param_filter_query(data_store,st_orientation=[0,numpy.pi/2],st_name='DriftingSinusoidalGratingDisk',st_contrast=100)    
        #dsv = param_filter_query(data_store,st_orientation=[0,numpy.pi/2],st_name=['FullfieldDriftingSinusoidalGrating'],st_contrast=100)    
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (16,12)},plot_file_name="OverviewEXCL4.png").plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (16,12)},plot_file_name="OverviewEXCL23.png").plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (16,12)},plot_file_name="OverviewINHL4.png").plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : l23_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (16,12)},plot_file_name="OverviewINHL23.png").plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        
        #dsv = param_filter_query(data_store,st_name='Null')    
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : l23_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : l23_exc, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})
        #OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : l23_inh, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot({'Vm_plot.y_lim' : (-67,-56),'Conductance_plot.y_lim' : (0,35.0)})


        #  PLOT ORIENTATION TUNING
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate'])    
        PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'contrast', 'neurons': l4_spike_exc_or_close_to_pi_half, 'sheet_name' : 'V1_Exc_L4'}),plot_file_name="CTC.png",fig_param={'dpi' : 200,'figsize': (14,3)}).plot({'TuningCurve_firing_rate.y_lim' : (0,50)})
        CTC(dsv,ParameterSet({'l4_exc_neurons' : numpy.array(l4_spike_exc_or_close_to_pi_half)[:4].tolist(),'l4_inh_neurons' : numpy.array(l4_spike_inh_or_close_to_pi_half)[:4].tolist(),'l23_exc_neurons' : numpy.array(l23_spike_exc_or_close_to_pi_half)[:4].tolist(),'l23_inh_neurons' : numpy.array(l23_spike_inh_or_close_to_pi_half)[:4].tolist()}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="CTC.svg").plot({'*.title' : None})
        #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'contrast', 'neurons': l4_spike_inh_or_close_to_pi_half, 'sheet_name' : 'V1_Inh_L4'})).plot({'TuningCurve_firing_rate.y_lim' : (0,50)})
        #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'contrast', 'neurons': l23_spike_exc_or_close_to_pi_half, 'sheet_name' : 'V1_Exc_L2/3'})).plot({'TuningCurve_firing_rate.y_lim' : (0,50)})
        #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'contrast', 'neurons': l23_spike_inh_or_close_to_pi_half, 'sheet_name' : 'V1_Inh_L2/3'})).plot({'TuningCurve_firing_rate.y_lim' : (0,50)})

        
        
        #  PLOT SIZETUNING
        #dsv = param_filter_query(data_store,st_name='DriftingSinusoidalGratingDisk',analysis_algorithm=['TrialAveragedFiringRate'])    
        #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'radius', 'neurons': l4_spike_exc_or_close_to_pi_half, 'sheet_name' : 'V1_Exc_L4'})).plot({'TuningCurve_firing_rate.y_lim' : (0,50)})
        #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'radius', 'neurons': l4_spike_inh_or_close_to_pi_half, 'sheet_name' : 'V1_Inh_L4'})).plot({'TuningCurve_firing_rate.y_lim' : (0,50)})
        #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'radius', 'neurons': l23_spike_exc_or_close_to_pi_half, 'sheet_name' : 'V1_Exc_L2/3'})).plot({'TuningCurve_firing_rate.y_lim' : (0,50)})
        #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'radius', 'neurons': l23_spike_inh_or_close_to_pi_half, 'sheet_name' : 'V1_Inh_L2/3'})).plot({'TuningCurve_firing_rate.y_lim' : (0,50)})
        #STC(dsv,ParameterSet({'l4_exc_neurons' : numpy.array(l4_spike_exc_or_close_to_pi_half)[:6].tolist(),'l4_inh_neurons' : numpy.array(l4_spike_inh_or_close_to_pi_half)[:6].tolist(),'l23_exc_neurons' : numpy.array(l23_spike_exc_or_close_to_pi_half)[:6].tolist(),'l23_inh_neurons' : numpy.array(l23_spike_inh_or_close_to_pi_half)[:6].tolist()}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="/home/antolikjan/Doc/Talks/JCUnicMar2013/STC.png").plot({'*.title' : None})
        
        # PLOT OCTC
        #dsv = param_filter_query(data_store,st_name='DriftingSinusoidalGratingCenterSurroundStimulus',analysis_algorithm=['TrialAveragedFiringRate'])    
        #OCTC(dsv,ParameterSet({'l4_exc_neurons' : numpy.array(l4_spike_exc_or_close_to_pi_half)[:5].tolist(),'l4_inh_neurons' : numpy.array(l4_spike_inh_or_close_to_pi_half)[:5].tolist(),'l23_exc_neurons' : numpy.array(l23_spike_exc_or_close_to_pi_half).tolist()[:5],'l23_inh_neurons' : numpy.array(l23_spike_inh_or_close_to_pi_half)[:5].tolist()}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="/home/antolikjan/Doc/Talks/JCUnicMar2013/OCTC.png").plot({'*.title' : None})
        
        # tuninc curves
        #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate','Analog_F0andF1'])    
        #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(l4_analog_ids), 'sheet_name' : 'V1_Exc_L4'})).plot()
        #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate','Analog_F0andF1'])    
        #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(l23_analog_ids), 'sheet_name' : 'V1_Exc_L2/3'})).plot()
        #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate','Analog_F0andF1'])    
        #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(l4_analog_ids_inh), 'sheet_name' : 'V1_Inh_L4'})).plot()
        #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate','Analog_F0andF1'])    
        #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(l23_analog_ids_inh), 'sheet_name' : 'V1_Inh_L2/3'})).plot()
        #OR(dsv,ParameterSet({'l4_exc_neurons' : numpy.array(l4_spike_exc_or_close_to_pi_half)[:4].tolist(),'l4_inh_neurons' : numpy.array(l4_spike_inh_or_close_to_pi_half)[:4].tolist(),'l23_exc_neurons' : numpy.array(l23_spike_exc_or_close_to_pi_half)[:4].tolist(),'l23_inh_neurons' : numpy.array(l23_spike_inh_or_close_to_pi_half)[:4].tolist()}),fig_param={'dpi' : 100,'figsize': (16,12)},plot_file_name="/home/antolikjan/Doc/Talks/JCUnicMar2013/OR.png").plot({'*.title' : None})
        
        #MRfig(data_store,ParameterSet({'SimpleSheetName' : 'V1_Exc_L4', 'ComplexSheetName' : 'V1_Exc_L2/3',}),fig_param={'dpi' : 100,'figsize': (8,6)},plot_file_name="/home/antolikjan/Doc/Talks/JCUnicMar2013/MR.png").plot()
        
        #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialVariability'])    
        #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(l4_analog_ids), 'sheet_name' : 'V1_Exc_L4'})).plot()

        #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialVariability'])    
        #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(l4_analog_ids_inh), 'sheet_name' : 'V1_Inh_L4'})).plot()


        import pylab
        pylab.show()







