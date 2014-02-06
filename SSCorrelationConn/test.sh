#!/bin/bash
#SBATCH -J MozaikParamSearchAnalysis
#SBATCH -c 1
source /opt/software/mpi/openmpi-1.6.3-gcc/env
source /home/antolikjan/env/mozaik/bin/activate
cd /home/antolikjan/dev/pkg/mozaik/mozaik/contrib/SelfSustainedPushPull
echo "DSADSA"
mpirun  --mca mtl ^psm python run_analysis.py '20131025-165227[param_AffON_WeakLGN_base_weight=0.0004_StrongerLGN2Inh_Small.defaults]CombinationParamSearch{l4_cortex_inh.rand_struct_ratio:10,l4_cortex_exc.L4ExcL4InhConnection.base_weight:5}//SelfSustainedPushPull_ParameterSearch_____l4_cortex_exc.L4ExcL4InhConnection.base_weight:0.0011_l4_cortex_inh.rand_struct_ratio:0.95' > '20131025-165227[param_AffON_WeakLGN_base_weight=0.0004_StrongerLGN2Inh_Small.defaults]CombinationParamSearch{l4_cortex_inh.rand_struct_ratio:10,l4_cortex_exc.L4ExcL4InhConnection.base_weight:5}//SelfSustainedPushPull_ParameterSearch_____l4_cortex_exc.L4ExcL4InhConnection.base_weight:0.0011_l4_cortex_inh.rand_struct_ratio:0.95/OUTFILE_analysis1382979184.69'
