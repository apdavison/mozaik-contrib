import mozaik
import sys
from mozaik.meta_workflow.visualization import single_value_visualization
from mozaik.controller import setup_logging
from mozaik.storage.queries import *
from parameters import ParameterSet

print len(sys.argv)

assert len(sys.argv) == 2
directory = sys.argv[1]

setup_logging()

single_value_visualization("VogeslAbbott2005",directory,
                           ParamFilterQuery(ParameterSet({'ads_unique' : False, 'rec_unique' : False, 'params' : ParameterSet({'sheet_name' : 'Exc1'})})),
                           value_names=None,filename='Exc.png',resolution=None,treat_nan_as_zero=True,ranges={'Mean(Firing rate)' : (0,20),'Mean(CV of ISI squared)' : (0,1.0), 'Mean(Correlation coefficient(psth (bin=5.0)))' : (0,0.2), 'Mean(Var(ECond)/Var(ICond))' : (0,1.0),'Mean(Var(ECond))' : (0,0.002),'Mean(Mean(ECond))' : (0,0.03),'Mean(Mean(ICond))' : (0,0.05),'Mean(Mean(ECond)/Mean(ICond))' : (0,0.8)})   

#single_value_visualization("VogeslAbbott2005",directory,
#                           ParamFilterQuery(ParameterSet({'ads_unique' : False, 'rec_unique' : False, 'params' : ParameterSet({'sheet_name' : 'V1_Inh_L4'})})),
#                           value_names=None,filename='Inh.png',resolution=None,treat_nan_as_zero=True,ranges={'Mean(Firing rate)' : (0,60),'Mean(CV of ISI squared)' : (0,2.0), 'Mean(Correlation coefficient(psth (bin=5.0)))' : (0,0.2), 'Mean(Var(ECond)/Var(ICond))' : (0,1.0),'Mean(Var(Econd))' : (0,0.002)})   
