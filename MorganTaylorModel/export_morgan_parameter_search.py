import sys
import subprocess
from mozaik.meta_workflow.parameter_search import *

master_results_dir = './20171012-170448[param_random_and_fast_afferent_delay_moreaff.defaults]CombinationParamSearch{12}';
run_script = 'export_to_morgan.py'
simulation_name = "MorganTaylorModel"

f = open(master_results_dir+'/parameter_combinations','rb')
combinations = pickle.load(f)
f.close()

# first check whether all parameter combinations contain the same parameter names
assert len(set([tuple(set(comb.keys())) for comb in combinations])) == 1 , "The parameter search didn't occur over a fixed set of parameters"

for i,combination in enumerate(combinations):
	rdn = master_results_dir+'/'+result_directory_name('ParameterSearch',simulation_name,combination)    

	subprocess.call(' '.join(["python",run_script,"'"+rdn+"'"]),shell=True)
    

#+['>']  + ["'"+rdn +'/OUTFILE_analysis'+str(time.time()) + "'"]

