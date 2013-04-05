import pylab
from mozaik.visualization.plotting import *
from mozaik.analysis.analysis import *
from mozaik.analysis.technical import NeuronAnnotationsToPerNeuronValues
from NeuroTools.parameters import ParameterSet
from mozaik.storage.queries import *
from mozaik.tools.circ_stat import *
from mozaik.tools.misc import *

def verify_push_pull(data_store,connection_name):
    # test the push pull 
    conn = data_store.get_analysis_result(identifier='Connections',proj_name=connection_name)
    assert len(conn) == 1
    conn = conn[0]
    OR_source = data_store.get_analysis_result(sheet_name=conn.source_name,identifier='PerNeuronValue',value_name='LGNAfferentOrientation')
    PHASE_source = data_store.get_analysis_result(sheet_name=conn.source_name,identifier='PerNeuronValue',value_name='LGNAfferentPhase')
    OR_target = data_store.get_analysis_result(sheet_name=conn.target_name,identifier='PerNeuronValue',value_name='LGNAfferentOrientation')
    PHASE_target = data_store.get_analysis_result(sheet_name=conn.target_name,identifier='PerNeuronValue',value_name='LGNAfferentPhase')
    
    assert len(OR_source) == 1
    assert len(PHASE_source) == 1
    assert len(OR_target) == 1
    assert len(PHASE_target) == 1
    OR_source = OR_source[0].values
    PHASE_source = PHASE_source[0].values
    OR_target = OR_target[0].values
    PHASE_target = PHASE_target[0].values

    pylab.figure()
    angles,mag = circ_mean(numpy.tile(OR_source,(conn.weights.shape[1],1)).T,conn.weights,high=numpy.pi,axis=0)
    pylab.subplot(1,2,1)        
    pylab.plot(angles,numpy.array(OR_target),'bo')
    pylab.xlabel('mean orientation of connections')
    pylab.ylabel('set orientation of neuron')
    
    angles,mag = circ_mean(numpy.tile(PHASE_source,(conn.weights.shape[1],1)).T,conn.weights,high=numpy.pi*2,axis=0)
    pylab.subplot(1,2,2)        
    pylab.plot(angles,numpy.array(PHASE_target),'bo')
    pylab.xlabel('mean phase of connections')
    pylab.ylabel('set phase of neuron')
    
    pylab.title(connection_name)


def verify_connectivity(data_store):
    l4_pos = data_store.get_neuron_postions()['V1_Exc_L4']
    l4_inh_pos = data_store.get_neuron_postions()['V1_Inh_L4']
    
    print find_neuron('center',l4_pos)
    print find_neuron('top_right',l4_pos)
    print find_neuron('top_left',l4_pos)
    print find_neuron('bottom_right',l4_pos)
    print find_neuron('bottom_left',l4_pos)
    
    ConnectivityPlot(data_store,ParameterSet({'neuron' : data_store.get_sheet_ids('V1_Exc_L4',find_neuron('center',l4_pos)),'sheet_name' : 'V1_Exc_L4','reversed' : True}),param_filter_query(data_store,identifier='PerNeuronValue',value_name='LGNAfferentOrientation')).plot()
    ConnectivityPlot(data_store,ParameterSet({'neuron' : data_store.get_sheet_ids('V1_Exc_L4',find_neuron('top_right',l4_pos)),'sheet_name' : 'V1_Exc_L4','reversed' : True}),param_filter_query(data_store,identifier='PerNeuronValue',value_name='LGNAfferentOrientation')).plot()
    ConnectivityPlot(data_store,ParameterSet({'neuron' : data_store.get_sheet_ids('V1_Exc_L4',find_neuron('top_left',l4_pos)),'sheet_name' : 'V1_Exc_L4','reversed' : True}),param_filter_query(data_store,identifier='PerNeuronValue',value_name='LGNAfferentOrientation')).plot()
    ConnectivityPlot(data_store,ParameterSet({'neuron' : data_store.get_sheet_ids('V1_Exc_L4',find_neuron('bottom_right',l4_pos)),'sheet_name' : 'V1_Exc_L4','reversed' : True}),param_filter_query(data_store,identifier='PerNeuronValue',value_name='LGNAfferentOrientation')).plot()
    ConnectivityPlot(data_store,ParameterSet({'neuron' : data_store.get_sheet_ids('V1_Exc_L4',find_neuron('bottom_left',l4_pos)),'sheet_name' : 'V1_Exc_L4','reversed' : True}),param_filter_query(data_store,identifier='PerNeuronValue',value_name='LGNAfferentOrientation')).plot()
    
    #ConnectivityPlot(data_store,ParameterSet({'neuron' : data_store.get_sheet_ids('V1_Exc_L4',find_neuron('center',l4_pos)),'sheet_name' : 'V1_Exc_L4','reversed' : True}),param_filter_query(data_store,identifier='PerNeuronValue',value_name='LGNAfferentPhase')).plot()
    #ConnectivityPlot(data_store,ParameterSet({'neuron' : data_store.get_sheet_ids('V1_Exc_L4',find_neuron('center',l4_inh_pos)),'sheet_name' : 'V1_Inh_L4','reversed' : True}),param_filter_query(data_store,identifier='PerNeuronValue',value_name='LGNAfferentOrientation')).plot()
    #ConnectivityPlot(data_store,ParameterSet({'neuron' : data_store.get_sheet_ids('V1_Exc_L4',find_neuron('center',l4_inh_pos)),'sheet_name' : 'V1_Inh_L4','reversed' : True}),param_filter_query(data_store,identifier='PerNeuronValue',value_name='LGNAfferentPhase')).plot()        
    
    verify_push_pull(data_store,'V1L4ExcL4ExcConnection')
    verify_push_pull(data_store,'V1L4ExcL4InhConnection')
    verify_push_pull(data_store,'V1L4InhL4ExcConnection')
    verify_push_pull(data_store,'V1L4InhL4InhConnection')
    
    #verify_push_pull(data_store,'V1ExcL23ExcL23Connection')
    #verify_push_pull(data_store,'V1ExcL23InhL23Connection')
    #verify_push_pull(data_store,'V1InhL23ExcL23Connection')
    #verify_push_pull(data_store,'V1InhL23InhL23Connection')

    #jens_model.connectors['V1ExcL23ExcL23Connection'].store_connections(data_store)    
    #jens_model.connectors['V1ExcL23InhL23Connection'].store_connections(data_store)    
    #jens_model.connectors['V1InhL23ExcL23Connection'].store_connections(data_store)    
    #jens_model.connectors['V1InhL23InhL23Connection'].store_connections(data_store)    


def F0_F1table(data_store,l4_exc):
    # PRINT INFO ON F0 and F1 of Conductances
    print "Excitatory cond, F0, preffered, contrast 100: ", str(param_filter_query(data_store,value_name=['F0_Exc_Cond'],sheet_name='V1_Exc_L4',st_orientation=[numpy.pi/2],st_contrast=100,ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc))
    print "Excitatory cond, F0, null, contrast 100: ", str(param_filter_query(data_store,value_name=['F0_Exc_Cond'],sheet_name='V1_Exc_L4',st_orientation=[0],st_contrast=100,ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc))
    #print "Excitatory cond, F0, preffered, contrast 30: ", str(param_filter_query(data_store,value_name=['F0_Exc_Cond'],sheet_name='V1_Exc_L4',st_orientation=[numpy.pi/2],st_contrast=30,ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc))
    #print "Excitatory cond, F0, null, contrast 30: ", str(param_filter_query(data_store,value_name=['F0_Exc_Cond'],sheet_name='V1_Exc_L4',st_orientation=[0],st_contrast=30,ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc))
    
    print "Excitatory cond, F1, preffered, contrast 100: ", str(param_filter_query(data_store,value_name=['F1_Exc_Cond'],sheet_name='V1_Exc_L4',st_orientation=[numpy.pi/2],st_contrast=100,ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc))
    print "Excitatory cond, F1, null, contrast 100: ", str(param_filter_query(data_store,value_name=['F1_Exc_Cond'],sheet_name='V1_Exc_L4',st_orientation=[0],st_contrast=100,ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc))
    #print "Excitatory cond, F1, preffered, contrast 30: ", str(param_filter_query(data_store,value_name=['F1_Exc_Cond'],sheet_name='V1_Exc_L4',st_orientation=[numpy.pi/2],st_contrast=30,ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc))
    #print "Excitatory cond, F1, null, contrast 30: ", str(param_filter_query(data_store,value_name=['F1_Exc_Cond'],sheet_name='V1_Exc_L4',st_orientation=[0],st_contrast=30,ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc))
      

    #print "Inhibitory cond, F0, preffered, contrast 100: ", str(param_filter_query(data_store,value_name=['F0_Inh_Cond'],sheet_name='V1_Exc_L4',st_orientation=[numpy.pi/2],st_contrast=100,ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc))
    #print "Inhibitory cond, F0, null, contrast 100: ", str(param_filter_query(data_store,value_name=['F0_Inh_Cond'],sheet_name='V1_Exc_L4',st_orientation=[0],st_contrast=100,ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc))
    #print "Inhibitory cond, F0, preffered, contrast 30: ", str(param_filter_query(data_store,value_name=['F0_Inh_Cond'],sheet_name='V1_Exc_L4',st_orientation=[numpy.pi/2],st_contrast=30,ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc))
    #print "Inhibitory cond, F0, null, contrast 30: ", str(param_filter_query(data_store,value_name=['F0_Exc_Cond'],sheet_name='V1_Exc_L4',st_orientation=[0],st_contrast=30,ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc))
    
    #print "Inhibitory cond, F1, preffered, contrast 100: ", str(param_filter_query(data_store,value_name=['F1_Inh_Cond'],sheet_name='V1_Exc_L4',st_orientation=[numpy.pi/2],st_contrast=100,ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc))
    #print "Inhibitory cond, F1, null, contrast 100: ", str(param_filter_query(data_store,value_name=['F1_Inh_Cond'],sheet_name='V1_Exc_L4',st_orientation=[0],st_contrast=100,ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc))
    #print "Inhibitory cond, F1, preffered, contrast 30: ", str(param_filter_query(data_store,value_name=['F1_Inh_Cond'],sheet_name='V1_Exc_L4',st_orientation=[numpy.pi/2],st_contrast=30,ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc))
    #print "Inhibitory cond, F1, null, contrast 30: ", str(param_filter_query(data_store,value_name=['F1_Inh_Cond'],sheet_name='V1_Exc_L4',st_orientation=[0],st_contrast=30,ads_unique=True).get_analysis_result()[0].get_value_by_id(l4_exc))

    
