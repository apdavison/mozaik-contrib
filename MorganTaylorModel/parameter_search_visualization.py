import matplotlib
matplotlib.use('Agg')

import sys
from mozaik.meta_workflow.analysis import export_SingleValues_as_matricies
from mozaik.meta_workflow.visualization import single_value_visualization
from mozaik.storage.queries import *
from parameters import ParameterSet

assert len(sys.argv) == 2
directory = sys.argv[1]

#single_value_visualization("MorganTaylorModel",directory,
#                           ParamFilterQuery(ParameterSet({'ads_unique' : False, 'rec_unique' : False, 'params' : ParameterSet({'sheet_name' : 'V1_Exc_L4'})})),value_names=None,filename='Exc.png',resolution=None,treat_nan_as_zero=True,ranges={'Mean HWHH of responsive neurons' : (15,30), 'Mean HWHH difference' : (0,15), 'Faild fits percantage' : (0,0.2)})   

export_SingleValues_as_matricies("MorganTaylorModel",directory,ParamFilterQuery(ParameterSet({'ads_unique' : False, 'rec_unique' : False, 'params' : ParameterSet({'sheet_name' : 'V1_Exc_L4'})})))   