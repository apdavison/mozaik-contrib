#!/usr/local/bin/ipython -i 
import sys
sys.path.append('/home/jan/projects/mozaik/')
import matplotlib
import time
import os
import pylab
from mozaik.framework.experiment_controller import run_experiments, setup_experiments, setup_logging
from mozaik.visualization.plotting import *
from mozaik.visualization.MRfig import *
from mozaik.analysis.analysis import *
from mozaik.analysis.technical import NeuronAnnotationsToPerNeuronValues
from mozaik.visualization.Kremkow_plots import *
from mozaik.storage.datastore import Hdf5DataStore,PickledDataStore
from NeuroTools.parameters import ParameterSet
from mozaik.storage.queries import *
from mozaik.tools.circ_stat import circular_dist
import mozaik

logger = mozaik.getMozaikLogger("Mozaik")

setup_logging()
data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':'ST'}))
logger.info('Loaded data store')
 
NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
# find neuron with preference closet to 0  

if True:  
    l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')
    l4_exc_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Exc_L4')
    l4_exc = numpy.argmin([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_exc_or[0].values[:20],l4_exc_phase[0].values[:20])])
    l4_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L4')
    l4_inh_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Inh_L4')
    l4_inh = numpy.argmin([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_inh_or[0].values[:20],l4_inh_phase[0].values[:20])])
    l4_exc_or_many = numpy.nonzero(numpy.array([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_exc_or[0].values,l4_exc_phase[0].values)]) < 0.1)[0]

    print "Prefered orientation of plotted exc neurons:"
    print l4_exc_or[0].values[l4_exc]
    print 'index ' + str(l4_exc)
    print "Prefered phase of plotted exc neurons:"
    print l4_exc_phase[0].values[l4_exc]
    print "Prefered orientation of plotted inh neurons:"
    print l4_inh_or[0].values[l4_inh]
    print 'index ' + str(l4_inh)
    print "Prefered phase of plotted inh neurons:"
    print l4_inh_phase[0].values[l4_inh]
    
    
os.mkdir('REPORT')

l4_exc_data = param_filter_query(data_store,sheet_name='V1_Exc_L4')
l4_inh_data = param_filter_query(data_store,sheet_name='V1_Inh_L4')

dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',contrast=100)    
OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),plot_file_name='REPORT/GratingsExc.png',fig_param={'dpi' : 100,'figsize': (14,8)}).plot(title=None)
OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),plot_file_name='REPORT/GratingsInh.png',fig_param={'dpi' : 100,'figsize': (14,8)}).plot(title=None)

dsv = param_filter_query(data_store,st_name='NaturalImageWithEyeMovement')    
OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),plot_file_name='REPORT/NIExc.png',fig_param={'dpi' : 100,'figsize': (14,8)}).plot(title=None)
OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),plot_file_name='REPORT/NIInh.png',fig_param={'dpi' : 100,'figsize': (14,8)}).plot(title=None)

dsv = param_filter_query(data_store,st_name='DriftingGratingWithEyeMovement')    
OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}}),plot_file_name='REPORT/GratingEMExc.png',fig_param={'dpi' : 100,'figsize': (14,8)}).plot(title=None)
OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}}),plot_file_name='REPORT/GratingEMInh.png',fig_param={'dpi' : 100,'figsize': (14,8)}).plot(title=None)

