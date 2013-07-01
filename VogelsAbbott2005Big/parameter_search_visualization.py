import sys
from mozaik.meta_workflow.visualization import single_value_visualization
from mozaik.storage.queries import *
from parameters import ParameterSet

assert len(sys.argv) == 2
directory = sys.argv[1]


single_value_visualization("VogeslAbbott2005Big",directory,
                           ParamFilterQuery(ParameterSet({'ads_unique' : False, 'rec_unique' : False, 'params' : ParameterSet({'sheet_name' : 'V1_Exc_L4'})})),
                           value_names=None,filename='Exc.png',resolution=30,treat_nan_as_zero=True,ranges={'Mean(Firing rate)' : (0,20), 'Mean(Correlation coefficient(psth (bin=5.0)))' : (0,0.2)})   

single_value_visualization("VogeslAbbott2005Big",directory,
                           ParamFilterQuery(ParameterSet({'ads_unique' : False, 'rec_unique' : False, 'params' : ParameterSet({'sheet_name' : 'V1_Inh_L4'})})),
                           value_names=None,filename='Inh.png',resolution=30,treat_nan_as_zero=True,ranges={'Mean(Firing rate)' : (0,20), 'Mean(Correlation coefficient(psth (bin=5.0)))' : (0,0.2)})   
#'Mean(CV of ISI squared)' : (0,0.1)
