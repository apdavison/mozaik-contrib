import sys
from mozaik.meta_workflow.visualization import single_value_visualization,multi_curve_visualzition
from mozaik.meta_workflow.analysis import export_SingleValues_as_matricies
from mozaik.storage.queries import *
from parameters import ParameterSet

assert len(sys.argv) == 2
directory = sys.argv[1]

#export_SingleValues_as_matricies("KumarEtAl2007",directory,ParamFilterQuery(ParameterSet({'ads_unique' : False, 'rec_unique' : False, 'params' : ParameterSet({'sheet_name' : 'V1_Exc_L4'})})))

#single_value_visualization("KumarEtAl2007",directory,
#                           ParamFilterQuery(ParameterSet({'ads_unique' : False, 'rec_unique' : False, 'params' : ParameterSet({'sheet_name' : 'V1_Exc_L4'})})),
#                           value_names=None,filename='Exc.png',resolution=20,treat_nan_as_zero=True,ranges={'Mean(Firing rate)' : (0,60)})   

multi_curve_visualzition("KumarEtAl2007",directory,'l4_cortex_exc.params.artificial_stimulators.background_act.params.exc_firing_rate',
                         ParamFilterQuery(ParameterSet({'ads_unique' : False, 'rec_unique' : False, 'params' : ParameterSet({'sheet_name' : 'V1_Exc_L4'})})),
                         value_name=None,filename='Exc.png',treat_nan_as_zero=True)   
