import pylab
import numpy

import mozaik.storage.queries as queries
import matplotlib.gridspec as gridspec
from mozaik.visualization.plotting import Plotting
from mozaik.visualization.helper_functions import *
from parameters import ParameterSet

from mozaik.controller import Global

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
          
          
          print len(dsv_simple)
          assert len(dsv_simple) == 1
          assert len(dsv_complex) == 1
          assert len(dsv_simple_v_F0) == 1
          assert len(dsv_complex_v_F0) == 1
          assert len(dsv_simple_v_F1) == 1
          assert len(dsv_complex_v_F1) == 1

          print dsv_simple_v_F0[0].values
          print dsv_simple_v_F1[0].values
          
          simple_v_mr = 2*dsv_simple_v_F1[0].values/abs(dsv_simple_v_F0[0].values)
          print simple_v_mr
          
          complex_v_mr = 2*dsv_complex_v_F1[0].values/abs(dsv_complex_v_F0[0].values)



          dsv_simple = dsv_simple[0]
          dsv_complex = dsv_complex[0]

          gs = gridspec.GridSpecFromSubplotSpec(3, 2,subplot_spec=gs)
          ax = pylab.subplot(gs[0,0])
          ax.hist(dsv_simple.values,bins=numpy.arange(0,2.2,0.2),color='k')
          pylab.ylim(0,450)
          disable_xticks(ax)
          remove_x_tick_labels()
          remove_y_tick_labels()
          #pylab.ylabel('Layer 4',fontsize=15)
          ax = pylab.subplot(gs[1,0])
          ax.hist(dsv_complex.values,bins=numpy.arange(0,2.2,0.2),color='w')
          pylab.ylim(0,450)
          disable_xticks(ax)
          remove_x_tick_labels()
          remove_y_tick_labels()
          #pylab.ylabel('Layer 2/3',fontsize=15)
          ax = pylab.subplot(gs[2,0])
          ax.hist([dsv_simple.values,dsv_complex.values],bins=numpy.arange(0,2.2,0.2),histtype='barstacked',color=['k','w'])
          pylab.ylim(0,450)
          #pylab.ylabel('Pooled',fontsize=15)
          three_tick_axis(ax.xaxis)
          remove_y_tick_labels()
          #pylab.xlabel('Modulation ratio',fontsize=15)
          for label in ax.get_xticklabels() + ax.get_yticklabels(): 
              label.set_fontsize(30) 

          
          if True:
              ax = pylab.subplot(gs[0,1])
              ax.hist(simple_v_mr,bins=numpy.arange(0,2.2,0.2),color='k')
              disable_xticks(ax)
              remove_x_tick_labels()
              remove_y_tick_labels()
              ax = pylab.subplot(gs[1,1])
              ax.hist(complex_v_mr,bins=numpy.arange(0,2.2,0.2),color='k')
              disable_xticks(ax)
              remove_x_tick_labels()
              remove_y_tick_labels()
              ax = pylab.subplot(gs[2,1])
              ax.hist([simple_v_mr,complex_v_mr],bins=numpy.arange(0,2.2,0.2),histtype='barstacked',color=['k','w'])
              three_tick_axis(ax.xaxis)
              remove_y_tick_labels()
              
          for label in ax.get_xticklabels() + ax.get_yticklabels(): 
              label.set_fontsize(30) 
          
          if self.plot_file_name:
                        pylab.savefig(Global.root_directory+self.plot_file_name)
        
