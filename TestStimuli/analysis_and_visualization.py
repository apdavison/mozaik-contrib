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
import sys
sys.path.append('/home/jan/projects/mozaik/contrib')
from Kremkow_plots import *

def perform_analysis_and_visualization(data_store):

    if False:
        import pylab
        pylab.show()


        on_analog = sorted(param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_esyn_ids())
        off_analog = sorted(param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_esyn_ids())

        TrialAveragedFiringRate(param_filter_query(data_store,st_name=['Null'],st_duration=100*147*7),ParameterSet({})).analyse()
        TrialAveragedFiringRate(param_filter_query(data_store,st_name=['FullfieldDriftingSinusoidalGrating']),ParameterSet({})).analyse()

        dsv = param_filter_query(data_store,st_name='Null',analysis_algorithm=['TrialAveragedFiringRate'],st_duration=100*147*7)    
        PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'background_luminance', 'neurons': list(on_analog)[:6], 'sheet_name' : 'X_ON','centered'  : False,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='LuminanceONLGN.png',fig_param={'dpi' : 100,'figsize': (32,15)}).plot()
        PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'background_luminance', 'neurons': list(off_analog)[:6], 'sheet_name' : 'X_OFF','centered'  : False,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='LuminanceOFtLGN.png',fig_param={'dpi' : 100,'figsize': (32,15)}).plot()
        dsv = param_filter_query(data_store,st_name=['Null','InternalStimulus'],st_duration=2*147*7)   
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : on_analog[0], 'sheet_activity' : {},'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="SpontLGN0On.png").plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : off_analog[0], 'sheet_activity' : {},'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="SpontLGN0Off.png").plot()


        dsv = param_filter_query(data_store,st_name=['Null','InternalStimulus'],st_duration=100*147*7)   
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : on_analog[0], 'sheet_activity' : {},'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LuminanceLGN0On.png").plot({'Vm_plot.y_lim' : (-100,-45)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : off_analog[0], 'sheet_activity' : {},'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LuminanceLGN0Off.png").plot({'Vm_plot.y_lim' : (-100,-45)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : on_analog[1], 'sheet_activity' : {},'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LuminanceLGN1On.png").plot({'Vm_plot.y_lim' : (-100,-45)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : off_analog[1], 'sheet_activity' : {},'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LuminanceLGN1Off.png").plot({'Vm_plot.y_lim' : (-100,-45)})

        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate'])    
        PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'contrast', 'neurons': list(on_analog)[:6], 'sheet_name' : 'X_ON','centered'  : False,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='ContrastONLGN.png',fig_param={'dpi' : 100,'figsize': (32,15)}).plot()
        PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'contrast', 'neurons': list(off_analog)[:6], 'sheet_name' : 'X_OFF','centered'  : False,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='ContrasOFtLGN.png',fig_param={'dpi' : 100,'figsize': (32,15)}).plot()

        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=100)   
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : on_analog[0], 'sheet_activity' : {},'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0On.png").plot({'Vm_plot.y_lim' : (-100,-45)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : off_analog[0], 'sheet_activity' : {},'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0Off.png").plot({'Vm_plot.y_lim' : (-100,-45)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : on_analog[1], 'sheet_activity' : {},'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN1On.png").plot({'Vm_plot.y_lim' : (-100,-45)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : off_analog[1], 'sheet_activity' : {},'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN1Off.png").plot({'Vm_plot.y_lim' : (-100,-45)})


        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=0,st_contrast=0)   
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : on_analog[0], 'sheet_activity' : {},'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0On_c0.png").plot({'Vm_plot.y_lim' : (-100,-45)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : off_analog[0], 'sheet_activity' : {},'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN0Off_c0.png").plot({'Vm_plot.y_lim' : (-100,-45)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : on_analog[1], 'sheet_activity' : {},'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN1On_c0.png").plot({'Vm_plot.y_lim' : (-100,-45)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : off_analog[1], 'sheet_activity' : {},'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="LGN1Off_c0.png").plot({'Vm_plot.y_lim' : (-100,-45)})

        dsv = param_filter_query(data_store,st_name='NaturalImageWithEyeMovement')   
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : on_analog[0], 'sheet_activity' : {},'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="NatLGN0On.png").plot({'Vm_plot.y_lim' : (-100,-45)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : off_analog[0], 'sheet_activity' : {},'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="NatLGN0Off.png").plot({'Vm_plot.y_lim' : (-100,-45)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : on_analog[1], 'sheet_activity' : {},'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="NatLGN1On.png").plot({'Vm_plot.y_lim' : (-100,-45)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : off_analog[1], 'sheet_activity' : {},'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (14,12)},plot_file_name="NatLGN1Off.png").plot({'Vm_plot.y_lim' : (-100,-45)})
    else:
        RetinalInputMovie(data_store,ParameterSet({}),plot_file_name="mov").plot()

        import pylab
        pylab.show()

