# -*- coding: utf-8 -*-
"""
"""
import sys
print sys.argv
import matplotlib
matplotlib.use('Agg')
import os
from mozaik.controller import setup_logging
import mozaik
from mozaik.storage.datastore import Hdf5DataStore,PickledDataStore
from mozaik.tools.mozaik_parametrized import colapse, colapse_to_dictionary, MozaikParametrized
from mozaik.analysis.technical import NeuronAnnotationsToPerNeuronValues
from parameters import ParameterSet
import numpy
from mozaik.storage import queries
from mozaik.controller import Global
Global.root_directory = sys.argv[1]+'/'

setup_logging()

data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':sys.argv[1],'store_stimuli' : False}),replace=True)

NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
analog_ids = queries.param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()


dsv = queries.param_filter_query(data_store,st_name='FlashedBar')
for ads in dsv.get_analysis_result():
    sid = MozaikParametrized.idd(ads.stimulus_id)
    sid.x=0
    ads.stimulus_id = str(sid)
for seg in dsv.get_segments():    
    sid = MozaikParametrized.idd(seg.annotations['stimulus'])
    sid.x=0
    seg.annotations['stimulus'] = str(sid)
for seg in dsv.get_segments(null=True):    
    sid = MozaikParametrized.idd(seg.annotations['stimulus'])
    sid.x=0
    seg.annotations['stimulus'] = str(sid)    


def save_data(dirname,dsv,name):

    try:
        os.mkdir(dirname)
    except:
        'all good'
        
    for neuron_id in analog_ids:
        mat_vm = []
        mat_exc = []
        mat_inh = []
        for seg in dsv.get_segments():
            sid = MozaikParametrized.idd(seg.annotations['stimulus'])
            a = seg.get_vm(neuron_id).magnitude
            a= numpy.insert(a,0,sid.trial)
            a= numpy.insert(a,0,sid.y)
            mat_vm.append(a)
            
            a = seg.get_esyn(neuron_id).magnitude
            a= numpy.insert(a,0,sid.trial)
            a= numpy.insert(a,0,sid.y)
            mat_exc.append(a)

            a = seg.get_isyn(neuron_id).magnitude
            a= numpy.insert(a,0,sid.trial)
            a= numpy.insert(a,0,sid.y)
            mat_inh.append(a)

            
        numpy.savetxt(dirname+'/'+'VM_' + name+str(neuron_id)+'.csv',numpy.array(mat_vm))
        numpy.savetxt(dirname+'/'+'ExcC' + name+str(neuron_id)+'.csv',numpy.array(mat_exc))
        numpy.savetxt(dirname+'/'+'InhC' + name+str(neuron_id)+'.csv',numpy.array(mat_inh))

        
d = os.path.split(sys.argv[1])[1]
dsv = queries.param_filter_query(data_store,st_name='FlashedBar',sheet_name='V1_Exc_L4',st_relative_luminance=0)    
save_data('./ForMorgan/'+d,dsv,'Dark_')
dsv = queries.param_filter_query(data_store,st_name='FlashedBar',sheet_name='V1_Exc_L4',st_relative_luminance=1.0)    
save_data('./ForMorgan/'+d,dsv,'Bright_')


