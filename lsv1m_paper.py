import pylab
import numpy

import mozaik.storage.queries as queries
import matplotlib.gridspec as gridspec
from mozaik.visualization.plotting import Plotting
from mozaik.visualization.helper_functions import *
from parameters import ParameterSet
from mozaik.storage.queries import *
from mozaik.controller import Global
from mozaik.visualization.plotting import (Plotting, GSynPlot,RasterPlot,PerNeuronAnalogSignalScatterPlot,
                                           VmPlot, ConductanceSignalListPlot,ScatterPlot,
                                           AnalogSignalListPlot,OverviewPlot,PerNeuronValueScatterPlot,PlotTuningCurve,PerNeuronValuePlot,CorticalColumnRasterPlot)


class MRfig(Plotting):
      required_parameters = ParameterSet({
            'SimpleSheetName' : str,  #the name of the sheet for which to plot
            'ComplexSheetName' : str, # which neuron to show
      })

      def plot(self):
          self.fig = pylab.figure(facecolor='w', **self.fig_param)
          gs = gridspec.GridSpec(1, 1)
          gs.update(left=0.07, right=0.97, top=0.9, bottom=0.1)
          gs = gs[0,0]
        
          dsv_simple = self.datastore.get_analysis_result(identifier='PerNeuronValue',sheet_name=self.parameters.SimpleSheetName,analysis_algorithm='ModulationRatio')
          dsv_complex = self.datastore.get_analysis_result(identifier='PerNeuronValue',sheet_name=self.parameters.ComplexSheetName,analysis_algorithm='ModulationRatio')
          
          
          dsv = queries.param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0)
          dsv_simple_v_F0 = dsv.get_analysis_result(identifier='PerNeuronValue',sheet_name=self.parameters.SimpleSheetName,value_name='F0_Vm-Mean(VM)')
          dsv_complex_v_F0 = dsv.get_analysis_result(identifier='PerNeuronValue',sheet_name=self.parameters.ComplexSheetName,value_name='F0_Vm-Mean(VM)')
          dsv_simple_v_F1 = dsv.get_analysis_result(identifier='PerNeuronValue',sheet_name=self.parameters.SimpleSheetName,value_name='F1_Vm')
          dsv_complex_v_F1 = dsv.get_analysis_result(identifier='PerNeuronValue',sheet_name=self.parameters.ComplexSheetName,value_name='F1_Vm')
          
          assert len(dsv_simple) == 1
          assert len(dsv_complex) == 1
          assert len(dsv_simple_v_F0) == 1
          assert len(dsv_complex_v_F0) == 1
          assert len(dsv_simple_v_F1) == 1
          assert len(dsv_complex_v_F1) == 1
            
          s_ids = dsv_simple_v_F0[0].ids
          c_ids = dsv_complex_v_F0[0].ids

          simple_v_mr = 2*numpy.array(dsv_simple_v_F1[0].get_value_by_id(s_ids))/abs(numpy.array(dsv_simple_v_F0[0].get_value_by_id(s_ids)))
          complex_v_mr = 2*numpy.array(dsv_complex_v_F1[0].get_value_by_id(c_ids))/abs(numpy.array(dsv_complex_v_F0[0].get_value_by_id(c_ids)))

          dsv_simple = dsv_simple[0]
          dsv_complex = dsv_complex[0]

          gs = gridspec.GridSpecFromSubplotSpec(3, 7,subplot_spec=gs,wspace=0.3)
          ax = pylab.subplot(gs[0,0])
          ax.hist(dsv_simple.values,bins=numpy.arange(0,2.2,0.2),color='w',rwidth=0.8)
          disable_top_right_axis(ax)
          disable_left_axis(ax)
          pylab.ylim(0,450)
          disable_xticks(ax)
          remove_x_tick_labels()
          remove_y_tick_labels()
          pylab.ylabel('Layer 4',fontsize=19)
          ax = pylab.subplot(gs[1,0])
          ax.hist(dsv_complex.values,bins=numpy.arange(0,2.2,0.2),color='k',rwidth=0.8)
          disable_top_right_axis(ax)
          disable_left_axis(ax)
          pylab.ylim(0,450)
          disable_xticks(ax)
          remove_x_tick_labels()
          remove_y_tick_labels()
          pylab.ylabel('Layer 2/3',fontsize=19)
          ax = pylab.subplot(gs[2,0])
          ax.hist([dsv_complex.values,dsv_simple.values],bins=numpy.arange(0,2.2,0.2),histtype='barstacked',color=['k','w'],rwidth=0.8)
          disable_top_right_axis(ax) 
          disable_left_axis(ax)  
          pylab.ylim(0,450)
          pylab.ylabel('Pooled',fontsize=19)
          three_tick_axis(ax.xaxis)
          remove_y_tick_labels()
          pylab.xlabel('F0/F1 spikes',fontsize=19)
          for label in ax.get_xticklabels() + ax.get_yticklabels(): 
              label.set_fontsize(19) 
          disable_top_right_axis(ax)
          disable_left_axis(ax)      

          ax = pylab.subplot(gs[0,1])
          ax.hist(simple_v_mr,bins=numpy.arange(0,3.3,0.3),color='w',rwidth=0.8)
          disable_top_right_axis(ax)    
          disable_left_axis(ax)      
          disable_xticks(ax)
          remove_x_tick_labels()
          remove_y_tick_labels()
          ax = pylab.subplot(gs[1,1])
          ax.hist(complex_v_mr,bins=numpy.arange(0,3.3,0.3),color='k',rwidth=0.8)
          disable_top_right_axis(ax)
          disable_left_axis(ax)
          disable_xticks(ax)
          remove_x_tick_labels()
          remove_y_tick_labels()
          ax = pylab.subplot(gs[2,1])
          ax.hist([complex_v_mr,simple_v_mr],bins=numpy.arange(0,3.3,0.3),histtype='barstacked',color=['k','w'],rwidth=0.8)
          three_tick_axis(ax.xaxis)
          remove_y_tick_labels()
          pylab.xlabel('F0/F1 Vm',fontsize=19)
          for label in ax.get_xticklabels() + ax.get_yticklabels(): 
              label.set_fontsize(19) 
          disable_top_right_axis(ax) 
          disable_left_axis(ax)                    
                
          ax = pylab.subplot(gs[0,2])
          ax.hist(abs(dsv_simple_v_F0[0].values),bins=numpy.arange(0,16,2.0),color='w',rwidth=0.8)
          disable_top_right_axis(ax)
          disable_left_axis(ax)                
          disable_left_axis(ax)      
          disable_xticks(ax)
          remove_x_tick_labels()
          remove_y_tick_labels()
          ax = pylab.subplot(gs[1,2])
          ax.hist(abs(dsv_complex_v_F0[0].values),bins=numpy.arange(0,16,2.0),color='k',rwidth=0.8)
          disable_top_right_axis(ax) 
          disable_left_axis(ax)                
          disable_xticks(ax)
          remove_x_tick_labels()
          remove_y_tick_labels()
          ax = pylab.subplot(gs[2,2])
          ax.hist([abs(dsv_complex_v_F0[0].values),abs(dsv_simple_v_F0[0].values)],bins=numpy.arange(0,16,2.0),histtype='barstacked',color=['k','w'],rwidth=0.8)
          three_tick_axis(ax.xaxis)
          remove_y_tick_labels()
          pylab.xlabel('F0 Vm (mV)',fontsize=19)
          for label in ax.get_xticklabels() + ax.get_yticklabels(): 
              label.set_fontsize(19) 
          disable_top_right_axis(ax)
          disable_left_axis(ax)

                    
          ax = pylab.subplot(gs[0,3])
          ax.hist(abs(dsv_simple_v_F1[0].values),bins=numpy.arange(0,16,2.0),color='w',rwidth=0.8)
          disable_top_right_axis(ax)      
          disable_left_axis(ax)
          disable_xticks(ax)
          remove_x_tick_labels()
          remove_y_tick_labels()
          ax = pylab.subplot(gs[1,3])
          ax.hist(abs(dsv_complex_v_F1[0].values),bins=numpy.arange(0,16,2.0),color='k',rwidth=0.8)
          disable_top_right_axis(ax)
          disable_left_axis(ax)
          disable_xticks(ax)
          remove_x_tick_labels()
          remove_y_tick_labels()
          ax = pylab.subplot(gs[2,3])
          ax.hist([abs(dsv_complex_v_F1[0].values),abs(dsv_simple_v_F0[0].values)],bins=numpy.arange(0,16,2.0),histtype='barstacked',color=['k','w'],rwidth=0.8)
          three_tick_axis(ax.xaxis)
          remove_y_tick_labels()
          pylab.xlabel('F1 Vm (mV)',fontsize=19)
          for label in ax.get_xticklabels() + ax.get_yticklabels(): 
              label.set_fontsize(19) 
          disable_top_right_axis(ax) 
          disable_left_axis(ax)
        
          ggs = gridspec.GridSpecFromSubplotSpec(20, 20, gs[:,4:7])
          ax = pylab.subplot(ggs[3:18,3:18])
          ax.plot(simple_v_mr,dsv_simple.get_value_by_id(s_ids),'ow',label='layer 4')
          ax.plot(complex_v_mr,dsv_complex.get_value_by_id(c_ids),'ok',label='layer 2/3')
          pylab.xlabel('F0/F1 Vm',fontsize=19)
          pylab.ylabel('F0/F1 Spikes',fontsize=19)
          pylab.xlim(0,5.0)  
          pylab.ylim(0,2.0)  
          for label in ax.get_xticklabels() + ax.get_yticklabels(): 
              label.set_fontsize(19) 
            
            
          if self.plot_file_name:
                        pylab.savefig(Global.root_directory+self.plot_file_name)


class LSV1MReponseOverview(Plotting):
    required_parameters = ParameterSet({
        'l4_exc_neuron' : int,
        'l4_inh_neuron' : int,
        'l23_exc_neuron' : int,
        'l23_inh_neuron' : int,
    })

    def subplot(self, subplotspec):
        plots = {}
        gs = gridspec.GridSpecFromSubplotSpec(17, 68, subplot_spec=subplotspec,hspace=1.0, wspace=100.0)
        
        dsv = param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=[0],st_contrast=[100])
        plots['ExcOr0L4'] = (OverviewPlot(dsv, ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : self.parameters.l4_exc_neuron, 'sheet_activity' : {}, 'spontaneous' : False})),gs[0:8,0:16],{'x_label': None,'x_axis' : False, 'x_ticks' : False })

        dsv = param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=[numpy.pi/2],st_contrast=[100])
        plots['ExcOrPiHL4'] = (OverviewPlot(dsv, ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : self.parameters.l4_exc_neuron, 'sheet_activity' : {}, 'spontaneous' : False})),gs[0:8,17:33],{'x_label': None,'y_label': None,'x_axis' : False, 'x_ticks' : False,'y_axis' : False, 'y_ticks' : False})

        dsv = param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=[0],st_contrast=[100])
        plots['InhOr0L4'] = (OverviewPlot(dsv, ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : self.parameters.l4_inh_neuron, 'sheet_activity' : {}, 'spontaneous' : False})),gs[0:8,35:51],{'x_label': None,'y_label': None, 'title' : None,'y_axis' : False, 'y_ticks' : False,'x_axis' : False, 'x_ticks' : False})

        dsv = param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=[numpy.pi/2],st_contrast=[100])
        plots['InhOrPiHL4'] = (OverviewPlot(dsv, ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : self.parameters.l4_inh_neuron, 'sheet_activity' : {}, 'spontaneous' : False})),gs[0:8,52:68],{'x_label': None,'y_label': None, 'title' : None,'y_axis' : False, 'y_ticks' : False,'x_axis' : False, 'x_ticks' : False})


        dsv = param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=[0],st_contrast=[100])
        plots['ExcOr0L23'] = (OverviewPlot(dsv, ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : self.parameters.l23_exc_neuron, 'sheet_activity' : {}, 'spontaneous' : False})),gs[9:17,0:16],{'title' : None})

        dsv = param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=[numpy.pi/2],st_contrast=[100])
        plots['ExcOrPiHL23'] = (OverviewPlot(dsv, ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : self.parameters.l23_exc_neuron, 'sheet_activity' : {}, 'spontaneous' : False})),gs[9:17,17:33],{'y_label': None, 'title' : None,'y_axis' : False, 'y_ticks' : False})

        dsv = param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=[0],st_contrast=[100])
        plots['InhOr0L23'] = (OverviewPlot(dsv, ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : self.parameters.l23_inh_neuron, 'sheet_activity' : {}, 'spontaneous' : False})),gs[9:17,35:51],{'y_label': None, 'title' : None,'y_axis' : False, 'y_ticks' : False})

        dsv = param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=[numpy.pi/2],st_contrast=[100])
        plots['InhOrPiHL23'] = (OverviewPlot(dsv, ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : self.parameters.l23_inh_neuron, 'sheet_activity' : {}, 'spontaneous' : False})),gs[9:17,52:68],{'y_label': None, 'title' : None,'y_axis' : False, 'y_ticks' : False})
        
        return plots


class SpontActOverview(Plotting):
    required_parameters = ParameterSet({
        'l4_exc_neuron' : int,
        'l4_inh_neuron' : int,
        'l23_exc_neuron' : int,
        'l23_inh_neuron' : int,
    })

    def subplot(self, subplotspec):
        plots = {}
        gs = gridspec.GridSpecFromSubplotSpec(8,3, subplot_spec=subplotspec,hspace=0.3, wspace=0.45)
        dsv = param_filter_query(data_store,st_direct_stimulation_name="None",st_name=['InternalStimulus'])    

        fontsize=17
        
        plots['SpikingOverview'] = (CorticalColumnRasterPlot(dsv,ParameterSet({'spontaneous' : False, 'sheet_names' : ['V1_Inh_L4','V1_Exc_L4','V1_Inh_L2/3','V1_Exc_L2/3'], 'colors' : ['#666666', '#000000' , '#666666', '#000000'], 'labels' : ["L4i","L4e" , "L2/3i", "L2/3e"]})),gs[:,0],{'fontsize' : fontsize})
        plots['ExcL2/3Cond'] = (GSynPlot(dsv, ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : self.parameters.l23_exc_neuron, 'spontaneous' : False})),gs[0,1:],{'x_label': None,'fontsize' : fontsize, 'x_ticks' : [],'title' : None})
        plots['ExcL2/3Vm'] = (VmPlot(dsv, ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neuron' : self.parameters.l23_exc_neuron, 'spontaneous' : False})),gs[1,1:],{'x_label': None,'fontsize' : fontsize, 'x_ticks' : [],'title' : None})
        plots['InhL2/3Cond'] = (GSynPlot(dsv, ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : self.parameters.l23_inh_neuron, 'spontaneous' : False})),gs[2,1:],{'x_label': None,'fontsize' : fontsize, 'x_ticks' : [],'title' : None})
        plots['InhL2/3Vm'] = (VmPlot(dsv, ParameterSet({'sheet_name' : 'V1_Inh_L2/3', 'neuron' : self.parameters.l23_inh_neuron, 'spontaneous' : False})),gs[3,1:],{'x_label': None,'fontsize' : fontsize, 'x_ticks' : [],'title' : None})
        plots['ExcL4Cond'] = (GSynPlot(dsv, ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : self.parameters.l4_exc_neuron, 'spontaneous' : False})),gs[4,1:],{'x_label': None,'fontsize' : fontsize, 'x_ticks' : [],'title' : None})
        plots['ExcL4Vm'] = (VmPlot(dsv, ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : self.parameters.l4_exc_neuron, 'spontaneous' : False})),gs[5,1:],{'x_label': None,'fontsize' : fontsize, 'x_ticks' : [],'title' : None})
        plots['InhL4Cond'] = (GSynPlot(dsv, ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : self.parameters.l4_inh_neuron, 'spontaneous' : False})),gs[6,1:],{'x_label': None,'fontsize' : fontsize, 'x_ticks' : [],'title' : None})
        plots['InhL4Vm'] = (VmPlot(dsv, ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : self.parameters.l4_inh_neuron, 'spontaneous' : False})),gs[7,1:],{'fontsize' : fontsize,'title' : None})
                
        return plots

class SpontStatisticsOverview(Plotting):
    required_parameters = ParameterSet({

    })

    def subplot(self, subplotspec):
        plots = {}
        gs = gridspec.GridSpecFromSubplotSpec(12,4, subplot_spec=subplotspec,hspace=10.0, wspace=0.5)
        dsv = param_filter_query(data_store,st_direct_stimulation_name="None",st_name=['InternalStimulus'])    

        fontsize=17
        
        mean_firing_rate_L4E = param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus',analysis_algorithm='PopulationMean',sheet_name='V1_Exc_L4',identifier='SingleValue',value_name='Mean(Firing rate)',ads_unique=True).get_analysis_result()[0].value
        mean_firing_rate_L4I = param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus',analysis_algorithm='PopulationMean',sheet_name='V1_Inh_L4',identifier='SingleValue',value_name='Mean(Firing rate)',ads_unique=True).get_analysis_result()[0].value
        mean_firing_rate_L23E = param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus',analysis_algorithm='PopulationMean',sheet_name='V1_Exc_L2/3',identifier='SingleValue',value_name='Mean(Firing rate)',ads_unique=True).get_analysis_result()[0].value
        mean_firing_rate_L23I = param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus',analysis_algorithm='PopulationMean',sheet_name='V1_Inh_L2/3',identifier='SingleValue',value_name='Mean(Firing rate)',ads_unique=True).get_analysis_result()[0].value
                
        mean_CV_L4E = param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus',analysis_algorithm='PopulationMean',sheet_name='V1_Exc_L4',identifier='SingleValue',value_name='Mean(CV of ISI squared)',ads_unique=True).get_analysis_result()[0].value
        mean_CV_L4I = param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus',analysis_algorithm='PopulationMean',sheet_name='V1_Inh_L4',identifier='SingleValue',value_name='Mean(CV of ISI squared)',ads_unique=True).get_analysis_result()[0].value
        mean_CV_L23E = param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus',analysis_algorithm='PopulationMean',sheet_name='V1_Exc_L2/3',identifier='SingleValue',value_name='Mean(CV of ISI squared)',ads_unique=True).get_analysis_result()[0].value
        mean_CV_L23I = param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus',analysis_algorithm='PopulationMean',sheet_name='V1_Inh_L2/3',identifier='SingleValue',value_name='Mean(CV of ISI squared)',ads_unique=True).get_analysis_result()[0].value
        
        mean_CC_L4E = param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus',analysis_algorithm='PopulationMean',sheet_name='V1_Exc_L4',identifier='SingleValue',value_name='Mean(Correlation coefficient(psth (bin=2.0)))',ads_unique=True).get_analysis_result()[0].value
        mean_CC_L4I = param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus',analysis_algorithm='PopulationMean',sheet_name='V1_Inh_L4',identifier='SingleValue',value_name='Mean(Correlation coefficient(psth (bin=2.0)))',ads_unique=True).get_analysis_result()[0].value
        mean_CC_L23E = param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus',analysis_algorithm='PopulationMean',sheet_name='V1_Exc_L2/3',identifier='SingleValue',value_name='Mean(Correlation coefficient(psth (bin=2.0)))',ads_unique=True).get_analysis_result()[0].value
        mean_CC_L23I = param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus',analysis_algorithm='PopulationMean',sheet_name='V1_Inh_L2/3',identifier='SingleValue',value_name='Mean(Correlation coefficient(psth (bin=2.0)))',ads_unique=True).get_analysis_result()[0].value
        
        mean_VM_L4E = numpy.mean(param_filter_query(data_store,sheet_name='V1_Exc_L4',st_direct_stimulation_name="None",st_name=['InternalStimulus'],analysis_algorithm='Analog_MeanSTDAndFanoFactor',value_name='Mean(VM)',ads_unique=True).get_analysis_result()[0].values)
        mean_VM_L4I= numpy.mean(param_filter_query(data_store,sheet_name='V1_Inh_L4',st_direct_stimulation_name="None",st_name=['InternalStimulus'],analysis_algorithm='Analog_MeanSTDAndFanoFactor',value_name='Mean(VM)',ads_unique=True).get_analysis_result()[0].values)
        mean_VM_L23E= numpy.mean(param_filter_query(data_store,sheet_name='V1_Exc_L2/3',st_direct_stimulation_name="None",st_name=['InternalStimulus'],analysis_algorithm='Analog_MeanSTDAndFanoFactor',value_name='Mean(VM)',ads_unique=True).get_analysis_result()[0].values)
        mean_VM_L23I = numpy.mean(param_filter_query(data_store,sheet_name='V1_Inh_L2/3',st_direct_stimulation_name="None",st_name=['InternalStimulus'],analysis_algorithm='Analog_MeanSTDAndFanoFactor',value_name='Mean(VM)',ads_unique=True).get_analysis_result()[0].values)
        
        mean_CondE_L4E = numpy.mean(param_filter_query(data_store,sheet_name='V1_Exc_L4',st_direct_stimulation_name="None",st_name=['InternalStimulus'],analysis_algorithm='Analog_MeanSTDAndFanoFactor',value_name='Mean(ECond)',ads_unique=True).get_analysis_result()[0].values)
        mean_CondE_L4I= numpy.mean(param_filter_query(data_store,sheet_name='V1_Inh_L4',st_direct_stimulation_name="None",st_name=['InternalStimulus'],analysis_algorithm='Analog_MeanSTDAndFanoFactor',value_name='Mean(ECond)',ads_unique=True).get_analysis_result()[0].values)
        mean_CondE_L23E= numpy.mean(param_filter_query(data_store,sheet_name='V1_Exc_L2/3',st_direct_stimulation_name="None",st_name=['InternalStimulus'],analysis_algorithm='Analog_MeanSTDAndFanoFactor',value_name='Mean(ECond)',ads_unique=True).get_analysis_result()[0].values)
        mean_CondE_L23I = numpy.mean(param_filter_query(data_store,sheet_name='V1_Inh_L2/3',st_direct_stimulation_name="None",st_name=['InternalStimulus'],analysis_algorithm='Analog_MeanSTDAndFanoFactor',value_name='Mean(ECond)',ads_unique=True).get_analysis_result()[0].values)

        mean_CondI_L4E = numpy.mean(param_filter_query(data_store,sheet_name='V1_Exc_L4',st_direct_stimulation_name="None",st_name=['InternalStimulus'],analysis_algorithm='Analog_MeanSTDAndFanoFactor',value_name='Mean(ICond)',ads_unique=True).get_analysis_result()[0].values)
        mean_CondI_L4I= numpy.mean(param_filter_query(data_store,sheet_name='V1_Inh_L4',st_direct_stimulation_name="None",st_name=['InternalStimulus'],analysis_algorithm='Analog_MeanSTDAndFanoFactor',value_name='Mean(ICond)',ads_unique=True).get_analysis_result()[0].values)
        mean_CondI_L23E= numpy.mean(param_filter_query(data_store,sheet_name='V1_Exc_L2/3',st_direct_stimulation_name="None",st_name=['InternalStimulus'],analysis_algorithm='Analog_MeanSTDAndFanoFactor',value_name='Mean(ICond)',ads_unique=True).get_analysis_result()[0].values)
        mean_CondI_L23I = numpy.mean(param_filter_query(data_store,sheet_name='V1_Inh_L2/3',st_direct_stimulation_name="None",st_name=['InternalStimulus'],analysis_algorithm='Analog_MeanSTDAndFanoFactor',value_name='Mean(ICond)',ads_unique=True).get_analysis_result()[0].values)
        
        
        pylab.rc('axes', linewidth=3)
        
        def plot_with_log_normal_fit(values,gs1,gs2,x_label=False,y_label=""):
            valuesnz = values[numpy.nonzero(values)[0]]
            h,bin_edges = numpy.histogram(numpy.log10(valuesnz),range=(-2,2),bins=20,normed=True)
            bin_centers = bin_edges[:-1] + (bin_edges[1:] - bin_edges[:-1])/2.0
            
            m = numpy.mean(numpy.log10(valuesnz))
            nm = numpy.mean(valuesnz)
            s = numpy.std(numpy.log10(valuesnz))

            pylab.subplot(gs1)
            pylab.plot(numpy.logspace(-2,2,100),numpy.exp(-((numpy.log10(numpy.logspace(-2,2,100))-m)**2)/(2*s*s))/(s*numpy.sqrt(2*numpy.pi)),linewidth=4,color="#666666")
            pylab.plot(numpy.power(10,bin_centers),h,'ko',mec=None,mew=3)
            pylab.xlim(10**-2,10**2)
            pylab.gca().set_xscale("log")
            if x_label:
                pylab.xlabel('firing rate [Hz]',fontsize=fontsize)
                pylab.xticks([0.01,0.1,1.0,10,100])
            else:
                pylab.xticks([])
            pylab.ylabel(y_label,fontsize=fontsize)                
            pylab.yticks([0.0,0.5,1.0])
            for label in pylab.gca().get_xticklabels() + pylab.gca().get_yticklabels():
                label.set_fontsize(fontsize)
            phf.disable_top_right_axis(pylab.gca())
            
            pylab.subplot(gs2)
            pylab.plot(numpy.logspace(-1,2,100),numpy.exp(-((numpy.log10(numpy.logspace(-1,2,100))-m)**2)/(2*s*s))/(s*numpy.sqrt(2*numpy.pi)),linewidth=4,color="#666666")
            pylab.plot(numpy.logspace(-1,2,100),numpy.exp(-numpy.logspace(-1,2,100)/nm)/nm,'k--',linewidth=4)
            pylab.plot(numpy.power(10,bin_centers),h,'ko',mec=None,mew=3)
            pylab.xlim(10**-1,10**2)
            pylab.ylim(0.00001,5.0)
            pylab.gca().set_xscale("log")
            pylab.gca().set_yscale("log")
            if x_label:
                pylab.xlabel('firing rate [Hz]',fontsize=fontsize)
                pylab.xticks([0.1,1.0,10,100])
            else:
                pylab.xticks([])
            pylab.yticks([0.0001,0.01,1.0])
            for label in pylab.gca().get_xticklabels() + pylab.gca().get_yticklabels():
                label.set_fontsize(fontsize)
            phf.disable_top_right_axis(pylab.gca())
        
        
        plot_with_log_normal_fit(param_filter_query(data_store,value_name=['Firing rate'],sheet_name=["V1_Exc_L4"],st_direct_stimulation_name="None",st_name=['InternalStimulus'],ads_unique=True).get_analysis_result()[0].values,gs[0:3,2],gs[0:3,3],y_label='L4e')
        plot_with_log_normal_fit(param_filter_query(data_store,value_name=['Firing rate'],sheet_name=["V1_Inh_L4"],st_direct_stimulation_name="None",st_name=['InternalStimulus'],ads_unique=True).get_analysis_result()[0].values,gs[3:6,2],gs[3:6,3],y_label='L4i')
        plot_with_log_normal_fit(param_filter_query(data_store,value_name=['Firing rate'],sheet_name=["V1_Exc_L2/3"],st_direct_stimulation_name="None",st_name=['InternalStimulus'],ads_unique=True).get_analysis_result()[0].values,gs[6:9,2],gs[6:9,3],y_label='L2/3e')
        plot_with_log_normal_fit(param_filter_query(data_store,value_name=['Firing rate'],sheet_name=["V1_Inh_L2/3"],st_direct_stimulation_name="None",st_name=['InternalStimulus'],ads_unique=True).get_analysis_result()[0].values,gs[9:12,2],gs[9:12,3],x_label=True,y_label='L2/3i')
        
        if True:
            pylab.subplot(gs[0:4,0])
            pylab.barh(numpy.array([0.17,0.67])-0.06,[mean_firing_rate_L4E,mean_firing_rate_L23E],height = 0.12,color='#000000')
            pylab.barh(numpy.array([0.33,0.83])-0.06,[mean_firing_rate_L4I,mean_firing_rate_L23I],height = 0.12,color='#666666')
            pylab.ylim(0,1.0)
            pylab.yticks([0.25,0.75],['L4','L2/3'])
            pylab.xlabel('firing rate (Hz)',fontsize=fontsize)
            phf.three_tick_axis(pylab.gca().xaxis)
            for label in pylab.gca().get_xticklabels() + pylab.gca().get_yticklabels():
                label.set_fontsize(fontsize)
            phf.disable_top_right_axis(pylab.gca())
            
            pylab.subplot(gs[4:8,0])
            pylab.barh(numpy.array([0.17,0.67])-0.06,[mean_CV_L4E,mean_CV_L23E],height = 0.12,color='#000000')
            pylab.barh(numpy.array([0.33,0.83])-0.06,[mean_CV_L4I,mean_CV_L23I],height = 0.12,color='#666666')
            pylab.ylim(0,1.0)
            pylab.yticks([0.25,0.75],['L4','L2/3'])
            pylab.xlabel('irregularity',fontsize=fontsize)
            phf.three_tick_axis(pylab.gca().xaxis)
            for label in pylab.gca().get_xticklabels() + pylab.gca().get_yticklabels():
                label.set_fontsize(fontsize)
            phf.disable_top_right_axis(pylab.gca())            

            pylab.subplot(gs[8:12,0])
            pylab.barh(numpy.array([0.17,0.67])-0.06,[mean_CC_L4E,mean_CC_L23E],height = 0.12,color='#000000')
            pylab.barh(numpy.array([0.33,0.83])-0.06,[mean_CC_L4I,mean_CC_L23I],height = 0.12,color='#666666')
            pylab.ylim(0,1.0)
            pylab.yticks([0.25,0.75],['L4','L2/3'])
            pylab.xlabel('synchrony',fontsize=fontsize)
            phf.three_tick_axis(pylab.gca().xaxis)
            for label in pylab.gca().get_xticklabels() + pylab.gca().get_yticklabels():
                label.set_fontsize(fontsize)
            phf.disable_top_right_axis(pylab.gca())


            
            pylab.subplot(gs[0:4,1])
            pylab.barh(numpy.array([0.17,0.67])-0.06,[mean_VM_L4E,mean_VM_L23E],height = 0.12,color='#000000')
            pylab.barh(numpy.array([0.33,0.83])-0.06,[mean_VM_L4I,mean_VM_L23I],height = 0.12,color='#666666')
            pylab.ylim(0,1.0)
            pylab.xlim(40,80)
            pylab.yticks([0.25,0.75],['L4','L2/3'])
            pylab.xlabel('membrane potential (mV)',fontsize=fontsize)
            phf.three_tick_axis(pylab.gca().xaxis)
            for label in pylab.gca().get_xticklabels() + pylab.gca().get_yticklabels():
                label.set_fontsize(fontsize)
            phf.disable_top_right_axis(pylab.gca())

            pylab.subplot(gs[4:8,1])
            pylab.barh(numpy.array([0.17,0.67])-0.06,[mean_CondE_L4E*1000,mean_CondE_L23E*1000],height = 0.12,color='#000000')
            pylab.barh(numpy.array([0.33,0.83])-0.06,[mean_CondE_L4I*1000,mean_CondE_L23I*1000],height = 0.12,color='#666666')
            pylab.ylim(0,1.0)
            pylab.yticks([0.25,0.75],['L4','L2/3'])
            pylab.xlabel('excitatory conductance (nS)',fontsize=fontsize)
            phf.three_tick_axis(pylab.gca().xaxis)
            for label in pylab.gca().get_xticklabels() + pylab.gca().get_yticklabels():
                label.set_fontsize(fontsize)
            phf.disable_top_right_axis(pylab.gca())            

            pylab.subplot(gs[8:12,1])
            pylab.barh(numpy.array([0.17,0.67])-0.06,[mean_CondI_L4E*1000,mean_CondI_L23E*1000],height = 0.12,color='#000000')
            pylab.barh(numpy.array([0.33,0.83])-0.06,[mean_CondI_L4I*1000,mean_CondI_L23I*1000],height = 0.12,color='#666666')
            pylab.ylim(0,1.0)
            pylab.yticks([0.25,0.75],['L4','L2/3'])
            pylab.xlabel('inhibitory conductance (nS)',fontsize=fontsize)
            phf.three_tick_axis(pylab.gca().xaxis)
            for label in pylab.gca().get_xticklabels() + pylab.gca().get_yticklabels():
                label.set_fontsize(fontsize)
            phf.disable_top_right_axis(pylab.gca())
            
            pylab.rc('axes', linewidth=1)
        
        return plots

        
def VMVarianceSummary():
       
        pylab.rc('axes', linewidth=3)
        
        
        
        def pl(dg,ni):
            values = {'DG' : dg, 'NI' : ni}
            
            pylab.figure(figsize=(3.7,5.3),dpi=400)
            pylab.subplots_adjust(left=.15)

            width = 1.0/len(values.keys())/2.0
            x = numpy.linspace(0+width,1-width,len(values.keys()))
            rects = pylab.bar(x-width/2.0,values.values(),width = width,color='k')
            pylab.xlim(0,1.0)
            pylab.ylim(-30,70)
            pylab.xticks(x,["DG","NI"])
            #pylab.yticks([-30,0,70],["70%","100%","170%"])
            pylab.yticks([-30,0,70],["","",""])
            pylab.axhline(0.0,color='k',linewidth=3)
            pylab.savefig("VMVarLayer23.png")


class OrientationTuningSummaryFiringRates(Plotting):
    required_parameters = ParameterSet({})

    required_parameters = ParameterSet({
        'exc_sheet_name1': str,  # the name of the sheet for which to plot
        'inh_sheet_name1': str,  # the name of the sheet for which to plot
        'exc_sheet_name2': str,  # the name of the sheet for which to plot
        'inh_sheet_name2': str,  # the name of the sheet for which to plot
    })


    def subplot(self, subplotspec):
        plots = {}
        gs = gridspec.GridSpecFromSubplotSpec(27, 38, subplot_spec=subplotspec,
                                              hspace=1.0, wspace=5.0)
        
        
        analog_ids1 = sorted(numpy.random.permutation(queries.param_filter_query(self.datastore,sheet_name=self.parameters.exc_sheet_name1).get_segments()[0].get_stored_esyn_ids()))
        analog_ids_inh1 = sorted(numpy.random.permutation(queries.param_filter_query(self.datastore,sheet_name=self.parameters.inh_sheet_name1).get_segments()[0].get_stored_esyn_ids()))
        analog_ids2 = sorted(numpy.random.permutation(queries.param_filter_query(self.datastore,sheet_name=self.parameters.exc_sheet_name2).get_segments()[0].get_stored_esyn_ids()))
        analog_ids_inh2 = sorted(numpy.random.permutation(queries.param_filter_query(self.datastore,sheet_name=self.parameters.inh_sheet_name2).get_segments()[0].get_stored_esyn_ids()))
        
        dsv = queries.param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate'])
        plots['ExcORTCMeanL4'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids1), 'sheet_name' : self.parameters.exc_sheet_name1,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[0:6,:6],{'title' : None,'x_label' : None , 'y_label' : 'Layer 4 (EXC)\n\nfiring rate (sp/s)', 'x_ticks' : None})
        plots['ExcORTC1L4'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids1[0:3]), 'sheet_name' : self.parameters.exc_sheet_name1,'centered'  : True,'mean' : False,'pool' : False,'polar' : False})),gs[0:3,6:15],{'title' : None,'left_border' : None, 'x_label' : None,'y_axis' : False,'x_axis' : False, 'x_ticks' : False})
        plots['ExcORTC2L4'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids1[3:6]), 'sheet_name' : self.parameters.exc_sheet_name1,'centered'  : True,'mean' : False,'pool' : False,'polar': False})),gs[3:6,6:15],{'title' : None,'left_border' : None, 'x_label' : None,'y_axis' : False,'x_axis' : False})

        plots['InhORTCMeanL4'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh1), 'sheet_name' : self.parameters.inh_sheet_name1,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[7:13,:6],{'title' : None, 'x_label' : None ,'y_label' : 'Layer 4 (INH)\n\nfiring rate (sp/s)', 'x_ticks' : None})
        plots['InhORTC1L4'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh1[0:3]), 'sheet_name' : self.parameters.inh_sheet_name1,'centered'  : True,'mean' : False,'pool' : False,'polar' : False})),gs[7:10,6:15],{'title' : None,'left_border' : None, 'x_label' : None,'y_axis' : False,'x_axis' : False})
        plots['InhORTC2L4'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh1[3:6]), 'sheet_name' : self.parameters.inh_sheet_name1,'centered'  : True,'mean' : False,'pool' : False,'polar' : False})),gs[10:13,6:15],{'title' : None,'left_border' : None, 'x_label' : None ,'y_axis' : None,'x_axis' : False, 'x_ticks' : False})

        plots['ExcORTCMeanL23'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids2), 'sheet_name' : self.parameters.exc_sheet_name2,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[14:20,:6],{'title' : None,'x_label' : None , 'y_label' : 'Layer 2/3 (EXC)\n\nfiring rate (sp/s)', 'x_ticks' : None})
        plots['ExcORTC1L23'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids2[0:3]), 'sheet_name' : self.parameters.exc_sheet_name2,'centered'  : True,'mean' : False,'pool' : False,'polar' : False})),gs[14:17,6:15],{'title' : None,'left_border' : None, 'x_label' : None,'y_axis' : False,'x_axis' : False, 'x_ticks' : False})
        plots['ExcORTC2L23'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids2[3:6]), 'sheet_name' : self.parameters.exc_sheet_name2,'centered'  : True,'mean' : False,'pool' : False,'polar': False})),gs[17:20,6:15],{'title' : None,'left_border' : None, 'x_label' : None,'y_axis' : False,'x_axis' : False})

        plots['InhORTCMeanL23'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh2), 'sheet_name' : self.parameters.inh_sheet_name2,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[21:27,:6],{'title' : None, 'y_label' : 'Layer 2/3 (INH)\n\nfiring rate (sp/s)'})
        plots['InhORTC1L23'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh2[0:3]), 'sheet_name' : self.parameters.inh_sheet_name2,'centered'  : True,'mean' : False,'pool' : False,'polar' : False})),gs[21:24,6:15],{'title' : None,'left_border' : None, 'x_label' : None,'y_axis' : False,'x_axis' : False})
        plots['InhORTC2L23'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh2[3:6]), 'sheet_name' : self.parameters.inh_sheet_name2,'centered'  : True,'mean' : False,'pool' : False,'polar' : False})),gs[24:27,6:15],{'title' : None,'left_border' : None, 'y_axis' : None,'x_axis' : False})

        
        dsv = queries.param_filter_query(self.datastore,value_name=['orientation HWHH of Firing rate'],sheet_name=[self.parameters.exc_sheet_name1])    
        plots['HWHHExcL4'] = (PerNeuronValueScatterPlot(dsv, ParameterSet({'only_matching_units' : True, 'ignore_nan' : True})),gs[0:6,17:23],{'x_lim': (0,50),'y_lim' : (0,50),'identity_line' : True, 'x_label' : None,'y_label' : 'HWHH Cont. 5%', 'cmp' : None,'title' : None})
        dsv = queries.param_filter_query(self.datastore,value_name=['orientation HWHH of Firing rate'],sheet_name=[self.parameters.inh_sheet_name1])    
        plots['HWHHInhL4'] = (PerNeuronValueScatterPlot(dsv, ParameterSet({'only_matching_units' : True, 'ignore_nan' : True})),gs[7:13,17:23],{'x_lim': (0,50),'y_lim' : (0,50),'identity_line' : True, 'x_label' : None,'y_label' : 'HWHH Cont. 5%', 'cmp' : None,'title' : None})
        dsv = queries.param_filter_query(self.datastore,value_name=['orientation HWHH of Firing rate'],sheet_name=[self.parameters.exc_sheet_name2])    
        plots['HWHHExcL23'] = (PerNeuronValueScatterPlot(dsv, ParameterSet({'only_matching_units' : True, 'ignore_nan' : True})),gs[14:20,17:23],{'x_lim': (0,50),'y_lim' : (0,50),'identity_line' : True, 'x_label' : None,'y_label' : 'HWHH Cont. 5%', 'cmp' : None,'title' : None})
        dsv = queries.param_filter_query(self.datastore,value_name=['orientation HWHH of Firing rate'],sheet_name=[self.parameters.inh_sheet_name2])    
        plots['HWHHInhL23'] = (PerNeuronValueScatterPlot(dsv, ParameterSet({'only_matching_units' : True, 'ignore_nan' : True})),gs[21:27,17:23],{'x_lim': (0,50),'y_lim' : (0,50),'identity_line' : True, 'x_label' : 'HWHH Cont. 100%','y_label' : 'HWHH Cont. 5%', 'cmp' : None,'title' : None})

        dsv = queries.param_filter_query(self.datastore,value_name=['orientation HWHH of Firing rate'],sheet_name=[self.parameters.exc_sheet_name1],st_contrast=[100])    
        plots['HWHHHistogramExcL4'] = (PerNeuronValuePlot(dsv, ParameterSet({'cortical_view' : False})),gs[0:6,25:31],{ 'x_lim' : (0.0,50.0), 'x_label' : None,'title' : None,'y_label' : '# neurons'})
        dsv = queries.param_filter_query(self.datastore,value_name=['orientation HWHH of Firing rate'],sheet_name=[self.parameters.inh_sheet_name1],st_contrast=[100])    
        plots['HWHHHistogramInhL4'] = (PerNeuronValuePlot(dsv, ParameterSet({'cortical_view' : False})),gs[7:13,25:31],{ 'x_lim' : (0.0,50.0), 'x_label' : None,'title' : None,'y_label' : '# neurons'})
        dsv = queries.param_filter_query(self.datastore,value_name=['orientation HWHH of Firing rate'],sheet_name=[self.parameters.exc_sheet_name2],st_contrast=[100])    
        plots['HWHHHistogramExcL23'] = (PerNeuronValuePlot(dsv, ParameterSet({'cortical_view' : False})),gs[14:20,25:31],{ 'x_lim' : (0.0,50.0), 'x_label' : None,'title' : None,'y_label' : '# neurons'})
        dsv = queries.param_filter_query(self.datastore,value_name=['orientation HWHH of Firing rate'],sheet_name=[self.parameters.inh_sheet_name2],st_contrast=[100])    
        plots['HWHHHistogramInhL23'] = (PerNeuronValuePlot(dsv, ParameterSet({'cortical_view' : False})),gs[21:27,25:31],{ 'x_lim' : (0.0,50.0), 'x_label' : 'HWHH (100% contrast)','title' : None,'y_label' : '# neurons'})

        dsv = queries.param_filter_query(self.datastore,value_name=['orientation CV(Firing rate)'],sheet_name=[self.parameters.exc_sheet_name1],st_contrast=[100])    
        plots['CVHistogramExcL4'] = (PerNeuronValuePlot(dsv, ParameterSet({'cortical_view' : False})),gs[0:6,32:38],{ 'x_lim' : (0.0,1.0), 'x_label' : None,'title' : None,'y_label' : None})
        dsv = queries.param_filter_query(self.datastore,value_name=['orientation CV(Firing rate)'],sheet_name=[self.parameters.inh_sheet_name1],st_contrast=[100])    
        plots['CVHistogramInhL4'] = (PerNeuronValuePlot(dsv, ParameterSet({'cortical_view' : False})),gs[7:13,32:38],{ 'x_lim' : (0.0,1.0), 'x_label' : None,'title' : None,'y_label' : None})
        dsv = queries.param_filter_query(self.datastore,value_name=['orientation CV(Firing rate)'],sheet_name=[self.parameters.exc_sheet_name2],st_contrast=[100])    
        plots['CVHistogramExcL23'] = (PerNeuronValuePlot(dsv, ParameterSet({'cortical_view' : False})),gs[14:20,32:38],{ 'x_lim' : (0.0,1.0), 'x_label' : None,'title' : None,'y_label' : None})
        dsv = queries.param_filter_query(self.datastore,value_name=['orientation CV(Firing rate)'],sheet_name=[self.parameters.inh_sheet_name2],st_contrast=[100])    
        plots['CVHistogramInhL23'] = (PerNeuronValuePlot(dsv, ParameterSet({'cortical_view' : False})),gs[21:27,32:38],{ 'x_lim' : (0.0,1.0), 'x_label' : 'CV (100% contrast)','title' : None,'y_label' : None})
        
        return plots

class OrientationTuningSummaryAnalogSignals(Plotting):

    required_parameters = ParameterSet({
        'exc_sheet_name1': str,  # the name of the sheet for which to plot
        'inh_sheet_name1': str,  # the name of the sheet for which to plot
        'exc_sheet_name2': str,  # the name of the sheet for which to plot
        'inh_sheet_name2': str,  # the name of the sheet for which to plot
    })


    def subplot(self, subplotspec):
        plots = {}
        gs = gridspec.GridSpecFromSubplotSpec(24,35, subplot_spec=subplotspec,
                                              hspace=10.0, wspace=0.5)
        
        analog_ids1 = sorted(numpy.random.permutation(queries.param_filter_query(self.datastore,sheet_name=self.parameters.exc_sheet_name1).get_segments()[0].get_stored_esyn_ids()))
        analog_ids_inh1 = sorted(numpy.random.permutation(queries.param_filter_query(self.datastore,sheet_name=self.parameters.inh_sheet_name1).get_segments()[0].get_stored_esyn_ids()))
        analog_ids2 = sorted(numpy.random.permutation(queries.param_filter_query(self.datastore,sheet_name=self.parameters.exc_sheet_name2).get_segments()[0].get_stored_esyn_ids()))
        analog_ids_inh2 = sorted(numpy.random.permutation(queries.param_filter_query(self.datastore,sheet_name=self.parameters.inh_sheet_name2).get_segments()[0].get_stored_esyn_ids()))
        
        # L4 EXC
        dsv = queries.param_filter_query(self.datastore,value_name=['F0_Vm-Mean(VM)'])    
        plots['L4E_F0_Vm'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids1), 'sheet_name' : self.parameters.exc_sheet_name1,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[0:6,:5],{'title' : None,'x_label' : None , 'y_label' : 'Layer 4 (EXC)','x_axis' : False, 'x_ticks' : False, 'title' : 'F0 of Vm (mV)'})
        dsv = queries.param_filter_query(self.datastore,value_name=['F1_Vm'])    
        plots['L4E_F1_Vm'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids1), 'sheet_name' : self.parameters.exc_sheet_name1,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[0:6,6:11],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False, 'title' : 'F1 of Vm (mV)'})
        dsv = queries.param_filter_query(self.datastore,value_name=['F0_Exc_Cond-Mean(ECond)'])    
        plots['L4E_F0_CondExc'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids1), 'sheet_name' : self.parameters.exc_sheet_name1,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[0:6,12:17],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False, 'title' : 'F0 of Exc. Cond. (nS)'})
        dsv = queries.param_filter_query(self.datastore,value_name=['F1_Exc_Cond'])    
        plots['L4E_F1_CondExc'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids1), 'sheet_name' : self.parameters.exc_sheet_name1,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[0:6,18:23],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False, 'title' : 'F1 of Exc. Cond. (nS)'})
        dsv = queries.param_filter_query(self.datastore,value_name=['F0_Inh_Cond-Mean(ICond)'])    
        plots['L4E_F0_CondInh'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids1), 'sheet_name' : self.parameters.exc_sheet_name1,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[0:6,24:29],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False, 'title' : 'F0 of Inh. Cond. (nS)'})
        dsv = queries.param_filter_query(self.datastore,value_name=['F1_Inh_Cond'])    
        plots['L4E_F1_CondInh'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids1), 'sheet_name' : self.parameters.exc_sheet_name1,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[0:6,30:35],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False, 'title' : 'F1 of Inh. Cond. (nS)'})

        # L4 INH
        dsv = queries.param_filter_query(self.datastore,value_name=['F0_Vm-Mean(VM)'])    
        plots['L4I_F0_Vm'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh1), 'sheet_name' : self.parameters.inh_sheet_name1,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[6:12,:5],{'title' : None,'x_label' : None , 'y_label' : 'Layer 4 (INH)','x_axis' : False, 'x_ticks' : False})
        dsv = queries.param_filter_query(self.datastore,value_name=['F1_Vm'])    
        plots['L4I_F1_Vm'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh1), 'sheet_name' : self.parameters.inh_sheet_name1,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[6:12,6:11],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False})
        dsv = queries.param_filter_query(self.datastore,value_name=['F0_Exc_Cond-Mean(ECond)'])    
        plots['L4I_F0_CondExc'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh1), 'sheet_name' : self.parameters.inh_sheet_name1,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[6:12,12:17],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False})
        dsv = queries.param_filter_query(self.datastore,value_name=['F1_Exc_Cond'])    
        plots['L4I_F1_CondExc'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh1), 'sheet_name' : self.parameters.inh_sheet_name1,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[6:12,18:23],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False})
        dsv = queries.param_filter_query(self.datastore,value_name=['F0_Inh_Cond-Mean(ICond)'])    
        plots['L4I_F0_CondInh'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh1), 'sheet_name' : self.parameters.inh_sheet_name1,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[6:12,24:29],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False})
        dsv = queries.param_filter_query(self.datastore,value_name=['F1_Inh_Cond'])    
        plots['L4I_F1_CondInh'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh1), 'sheet_name' : self.parameters.inh_sheet_name1,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[6:12,30:35],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False})
                

        # L2/3 EXC
        dsv = queries.param_filter_query(self.datastore,value_name=['F0_Vm-Mean(VM)'])    
        plots['L23E_F0_Vm'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids2), 'sheet_name' : self.parameters.exc_sheet_name2,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[12:18,:5],{'title' : None,'x_label' : None , 'y_label' : 'Layer 2/3 (EXC)','x_axis' : False, 'x_ticks' : False})
        dsv = queries.param_filter_query(self.datastore,value_name=['F1_Vm'])    
        plots['L23E_F1_Vm'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids2), 'sheet_name' : self.parameters.exc_sheet_name2,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[12:18,6:11],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False})
        dsv = queries.param_filter_query(self.datastore,value_name=['F0_Exc_Cond-Mean(ECond)'])    
        plots['L23E_F0_CondExc'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids2), 'sheet_name' : self.parameters.exc_sheet_name2,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[12:18,12:17],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False})
        dsv = queries.param_filter_query(self.datastore,value_name=['F1_Exc_Cond'])    
        plots['L23E_F1_CondExc'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids2), 'sheet_name' : self.parameters.exc_sheet_name2,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[12:18,18:23],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False})
        dsv = queries.param_filter_query(self.datastore,value_name=['F0_Inh_Cond-Mean(ICond)'])    
        plots['L23E_F0_CondInh'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids2), 'sheet_name' : self.parameters.exc_sheet_name2,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[12:18,24:29],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False})
        dsv = queries.param_filter_query(self.datastore,value_name=['F1_Inh_Cond'])    
        plots['L23E_F1_CondInh'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids2), 'sheet_name' : self.parameters.exc_sheet_name2,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[12:18,30:35],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False})

        # L2/3 INH
        dsv = queries.param_filter_query(self.datastore,value_name=['F0_Vm-Mean(VM)'])    
        plots['L23I_F0_Vm'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh2), 'sheet_name' : self.parameters.inh_sheet_name2,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[18:24,:5],{'title' : None,'x_label' : None , 'y_label' : 'Layer 2/3 (INH)'})
        dsv = queries.param_filter_query(self.datastore,value_name=['F1_Vm'])    
        plots['L23I_F1_Vm'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh2), 'sheet_name' : self.parameters.inh_sheet_name2,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[18:24,6:11],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False})
        dsv = queries.param_filter_query(self.datastore,value_name=['F0_Exc_Cond-Mean(ECond)'])    
        plots['L23I_F0_CondExc'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh2), 'sheet_name' : self.parameters.inh_sheet_name2,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[18:24,12:17],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False})
        dsv = queries.param_filter_query(self.datastore,value_name=['F1_Exc_Cond'])    
        plots['L23I_F1_CondExc'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh2), 'sheet_name' : self.parameters.inh_sheet_name2,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[18:24,18:23],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False})
        dsv = queries.param_filter_query(self.datastore,value_name=['F0_Inh_Cond-Mean(ICond)'])    
        plots['L23I_F0_CondInh'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh2), 'sheet_name' : self.parameters.inh_sheet_name2,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[18:24,24:29],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False})
        dsv = queries.param_filter_query(self.datastore,value_name=['F1_Inh_Cond'])    
        plots['L23I_F1_CondInh'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh2), 'sheet_name' : self.parameters.inh_sheet_name2,'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[18:24,30:35],{'title' : None,'x_label' : None ,'y_label' : None,'x_axis' : False, 'x_ticks' : False})
                
                    
        
        return plots

class TrialToTrialVariabilityComparison(Plotting):

    required_parameters = ParameterSet({
        'sheet_name1' : str, # The name of the sheet in which to do the analysis
        'sheet_name2' : str, # The name of the sheet in which to do the analysis
        'data_ni' : float,
        'data_dg' : float,
    })


    def plot(self):
        self.fig = pylab.figure(facecolor='w', **self.fig_param)
        gs = gridspec.GridSpec(1, 3)
        gs.update(left=0.07, right=0.97, top=0.9, bottom=0.1,wspace=0.1)
        
        orr = list(set([MozaikParametrized.idd(s).orientation for s in queries.param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGrating',st_contrast=100).get_stimuli()]))        
        l4_exc_or = self.datastore.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = self.parameters.sheet_name1)
        l23_exc_or = self.datastore.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = self.parameters.sheet_name2)
        
        
        # lets calculate spont. activity trial to trial variability
        # we assume that the spontaneous activity had already the spikes removed
        
        def calculate_sp(datastore,sheet_name):
            dsv = queries.param_filter_query(datastore,st_name='InternalStimulus',st_direct_stimulation_name='None',sheet_name=sheet_name,analysis_algorithm='ActionPotentialRemoval',ads_unique=True)
            ids = dsv.get_analysis_result()[0].ids
            sp= {}
            for idd in ids:
                assert len(dsv.get_analysis_result()) == 1
                s = dsv.get_analysis_result()[0].get_asl_by_id(idd).magnitude
                sp[idd] = 1/numpy.mean(numpy.std([s[i*int(len(s)/10):(i+1)*int(len(s)/10)] for i in xrange(0,10)],axis=0,ddof=1))

            return sp
        
        sp_l4 = calculate_sp(self.datastore,self.parameters.sheet_name1)
        sp_l23 = calculate_sp(self.datastore,self.parameters.sheet_name2)
        
        def calculate_var_ratio(datastore,sheet_name,sp,ors):
            #lets calculate the mean of trial-to-trial variances across the neurons in the datastore for gratings 
            dsv = queries.param_filter_query(datastore,st_name='FullfieldDriftingSinusoidalGrating',sheet_name=sheet_name,st_contrast=100,analysis_algorithm='TrialVariability',y_axis_name='Vm (no AP) trial-to-trial variance')
            assert queries.equal_ads(dsv, except_params=['stimulus_id'])
            ids = dsv.get_analysis_result()[0].ids
            
            std_gr = 0
            
            for i in ids:
                # find the or pereference of the neuron
                o = orr[numpy.argmin([circular_dist(o,ors[0].get_value_by_id(i),numpy.pi) for o in orr])]
                assert len(queries.param_filter_query(dsv,st_orientation=o,ads_unique=True).get_analysis_result())==1
                a = 1/numpy.mean(numpy.sqrt(queries.param_filter_query(dsv,st_orientation=o,ads_unique=True).get_analysis_result()[0].get_asl_by_id(i).magnitude))
                std_gr = std_gr + a / sp[i]

            std_gr = std_gr / len(ids)

            #lets calculate the mean of trial-to-trial variances across the neurons in the datastore for natural images 
            dsv = queries.param_filter_query(datastore,st_name='NaturalImageWithEyeMovement',sheet_name=sheet_name,y_axis_name='Vm (no AP) trial-to-trial variance',ads_unique=True)
            std_ni_ind = [1/numpy.mean(numpy.sqrt(dsv.get_analysis_result()[0].get_asl_by_id(i).magnitude)) / sp[i] for i in ids]
            std_ni = numpy.mean(std_ni_ind)
            
            return std_gr,std_ni
        
        var_gr_l4,var_ni_l4 = calculate_var_ratio(self.datastore,self.parameters.sheet_name1,sp_l4,l4_exc_or)
        var_gr_l23,var_ni_l23 = calculate_var_ratio(self.datastore,self.parameters.sheet_name2,sp_l23,l23_exc_or)
    
        
        lw = pylab.rcParams['axes.linewidth']
        pylab.rc('axes', linewidth=3)
        width = 0.25
        x = numpy.array([width,1-width])

        def plt(a,b):
            print [a,b]
            rects = pylab.bar(x-width/2.0,[a*100-100,b*100-100],width = width,color='k')
            pylab.xlim(0,1.0)
            pylab.ylim(-30,70)
            pylab.xticks(x,["DG","NI"])
            pylab.yticks([-30,0,70],["70%","100%","170%"])
            pylab.axhline(0.0,color='k',linewidth=3)
            disable_top_right_axis(pylab.gca())    
            disable_xticks(pylab.gca())
            for label in pylab.gca().get_xticklabels() + pylab.gca().get_yticklabels():
                label.set_fontsize(19)
        
        ax = pylab.subplot(gs[0,0])
        plt(self.parameters.data_dg,self.parameters.data_ni)        
        pylab.title("Data",fontsize=19,y=1.05)
        
            
        ax = pylab.subplot(gs[0,1])
        plt(var_gr_l4,var_ni_l4)
        disable_left_axis(ax)
        remove_y_tick_labels()
        pylab.title("Layer 4",fontsize=19,y=1.05)

        ax = pylab.subplot(gs[0,2])
        plt(var_gr_l23,var_ni_l23)
        disable_left_axis(ax)
        remove_y_tick_labels()
        pylab.title("Layer 2/3",fontsize=19,y=1.05)
        
        pylab.rc('axes', linewidth=lw)
        
