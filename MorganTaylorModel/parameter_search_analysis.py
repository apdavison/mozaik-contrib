import matplotlib
matplotlib.use('Agg')

import sys
from mozaik.meta_workflow.parameter_search import parameter_search_run_script_distributed_slurm_IoV
assert len(sys.argv) == 2, "Wrong number of arguments, required 2 supplied %g" % (len(sys.argv))
directory = sys.argv[1]

parameter_search_run_script_distributed_slurm_IoV("MorganTaylorModel",directory,'run_analysis.py',8)
