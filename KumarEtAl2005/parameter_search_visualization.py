import sys
from mozaik.meta_workflow.visualization import single_value_visualization
from mozaik.meta_workflow.analysis import collect_results_from_parameter_search
from mozaik.storage.queries import *
from parameters import ParameterSet

assert len(sys.argv) == 2
directory = sys.argv[1]


single_value_visualization("KumarEtAl2007",directory,
                           ParamFilterQuery(ParameterSet({'ads_unique' : False, 'rec_unique' : False, 'params' : ParameterSet({'sheet_name' : 'V1_Exc_L4'})})),
                           value_names=None,filename='Exc.png',resolution=20,treat_nan_as_zero=True,ranges={'Mean(Firing rate)' : (0,60),'Mean(CV of ISI squared)' : (0,1.0), 'Mean(Correlation coefficient(psth (bin=5.0)))' : (0,0.2)})   

single_value_visualization("KumarEtAl2007",directory,
                           ParamFilterQuery(ParameterSet({'ads_unique' : False, 'rec_unique' : False, 'params' : ParameterSet({'sheet_name' : 'V1_Inh_L4'})})),
                           value_names=None,filename='Inh.png',resolution=20,treat_nan_as_zero=True,ranges={'Mean(Firing rate)' : (0,60),'Mean(CV of ISI squared)' : (0,1.0), 'Mean(Correlation coefficient(psth (bin=5.0)))' : (0,0.2)})   


def fexc(datastore):
    a = param_filter_query(datastore,sheet_name='V1_Exc_L4',st_direct_stimulation_name="None")
    assert len(a.get_segments()) == 1
    d  = numpy.array([numpy.array(z) for z in a.get_segments()[0].get_spiketrains()])
    print numpy.shape(d)
    return d

def finh(datastore):
    a = param_filter_query(datastore,sheet_name='V1_Inh_L4',st_direct_stimulation_name="None")
    assert len(a.get_segments()) == 1
    d  = numpy.array([numpy.array(z) for z in a.get_segments()[0].get_spiketrains()])
    print numpy.shape(d)
    return d


collect_results_from_parameter_search("KumarEtAl2007",directory,fexc,"ResExc.pickle")
collect_results_from_parameter_search("KumarEtAl2007",directory,finh,"ResInh.pickle")
