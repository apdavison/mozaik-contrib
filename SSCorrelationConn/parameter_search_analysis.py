import sys
from mozaik.meta_workflow.parameter_search import parameter_search_run_script_distributed_slurm
assert len(sys.argv) == 2
directory = sys.argv[1]

parameter_search_run_script_distributed_slurm("SSCorrelationConnectivity",directory,'run_analysis.py',16)
