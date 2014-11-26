import sys
from mozaik.meta_workflow.visualization import single_value_visualization
from mozaik.meta_workflow.analysis import collect_results_from_parameter_search
from mozaik.storage.queries import *
from parameters import ParameterSet

assert len(sys.argv) == 2
directory = sys.argv[1]


single_value_visualization("SelfSustainedPushPull",directory,
                           ParamFilterQuery(ParameterSet({'ads_unique' : False, 'rec_unique' : False, 'params' : ParameterSet({'sheet_name' : 'V1_Exc_L4','st_name' : 'InternalStimulus', 'st_direct_stimulation_name' : 'None'})})),
                           value_names=None,filename='Exc.png',resolution=20,treat_nan_as_zero=True,ranges={'Mean(Firing rate)' : (0,60),'Mean(CV of ISI squared)' : (0,1.0)})   

#single_value_visualization("SelfSustainedPushPull",directory,
#                           ParamFilterQuery(ParameterSet({'ads_unique' : False, 'rec_unique' : False, 'params' : ParameterSet({'sheet_name' : 'V1_Inh_L4','st_name' : 'InternalStimulus', 'st_direct_stimulation_name' : 'None'})})),
#                           value_names=None,filename='Inh.png',resolution=20,treat_nan_as_zero=True,ranges={'Mean(Firing rate)' : (0,60),'Mean(CV of ISI squared)' : (0,1.0)})   


