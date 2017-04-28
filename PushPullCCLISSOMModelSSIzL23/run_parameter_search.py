# -*- coding: utf-8 -*-
import sys
from mozaik.meta_workflow.parameter_search import CombinationParameterSearch,SlurmSequentialBackend
import numpy
import time

if False:
    CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=64),{
									     'retina_lgn.params.gain' : [0.1],
									     'l4_cortex_exc.params.density' : [10],
									     }).run_parameter_search()





if False:
    CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=64),{
									     'l4_cortex_exc.AfferentConnection.base_weight' : [0.0015],
									     'l23_cortex_exc.L23ExcL23InhConnection.base_weight' : [0.003],
									     'l23_cortex_exc.L4ExcL23ExcConnection.base_weight' : [0.003],
									     'l23_cortex_inh.L4ExcL23InhConnection.base_weight' : [0.0001],
									     'l23_cortex_inh.L23InhL23ExcConnection.base_weight' : [0.0025],
									     'l23_cortex_inh.L23InhL23InhConnection.base_weight' : [0.0017],
									     'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0004],
									     'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.002,0.0025,0.003],
									     'l4_cortex_inh.ExcInhAfferentRatio' : [1.7],
									     'l4_cortex_exc.params.density' : [300],
									     'only_afferent' : [False],
									     'l4_cortex_inh.L4InhL4ExcConnection.short_term_plasticity.U' : [0.1,0.13,0.16],
									     }).run_parameter_search()

if False:
    CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=64),{
									     'l23_cortex_exc.L23ExcL23InhConnection.base_weight' : [0.002,0.001],
									     'l23_cortex_exc.L4ExcL23ExcConnection.base_weight' : [0.002,0.001],
									     'l23_cortex_inh.L4ExcL23InhConnection.base_weight' : [0.0001,0.001],
									     'l23_cortex_inh.L23InhL23ExcConnection.base_weight' : [0.0025,0.003,0.0035],
									     'l23_cortex_inh.L23InhL23InhConnection.base_weight' : [0.0017],
									     'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.0005],
									     'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0007,0.00075],
									     'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0018],
									     'l4_cortex_inh.ExcInhAfferentRatio' : [1.4,1.3],
									     'l4_cortex_exc.params.density' : [300],
									     'l4_cortex_inh.L4InhL4ExcConnection.short_term_plasticity.tau_rec' : [25],
									     }).run_parameter_search()


if False:
    CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=64),{
									     'l4_cortex_exc.AfferentConnection.base_weight' : [0.002],
									     'l4_cortex_exc.AfferentConnection.num_samples' : [15,20],
									     'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0007],
									     'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.00065],
									     'l23_cortex_exc.L23ExcL23InhConnection.base_weight' : [0.0015],
									     'l23_cortex_inh.L23InhL23ExcConnection.base_weight' : [0.003],
									     'l23_cortex_exc.L4ExcL23ExcConnection.base_weight' : [0.0015],
									     'l4_cortex_inh.ExcInhAfferentRatio' : [0.7,0.8,0.9],
									     'l4_cortex_exc.params.density' : [300],
									     'l23_cortex_exc.params.density' : [300],
									     'l4_cortex_inh.L4InhL4ExcConnection.short_term_plasticity.tau_fac' : [300,250,350],
									     'l4_cortex_inh.L4InhL4ExcConnection.short_term_plasticity.U' : [0.1,0.11,0.9],
									     }).run_parameter_search()



if False:
    CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=64),{
									     'l4_cortex_exc.AfferentConnection.base_weight' : [0.0005,0.0006],
									     'l4_cortex_exc.AfferentConnection.num_samples' : [45,50],
									     'l23_cortex_exc.L23ExcL23InhConnection.base_weight' : [0.002],
									     'l23_cortex_exc.L4ExcL23ExcConnection.base_weight' : [0.0024],
									     'l23_cortex_inh.L4ExcL23InhConnection.base_weight' : [0.0001],
									     'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0026],
									     'l4_cortex_inh.L4InhL4InhConnection.base_weight' : [0.0002,0.0003],
									     'l4_cortex_inh.ExcInhAfferentRatio' : [1.25,1.2],
									     'l4_cortex_exc.AfferentConnection.short_term_plasticity.tau_rec' : [60],
									     'l4_cortex_exc.params.density' : [300],
									     'l4_cortex_exc.params.sx' : [6000],
									     }).run_parameter_search()



if False:
    CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=64),{
									     'l4_cortex_exc.AfferentConnection.base_weight' : [0.0007],
									     'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.0007],
									     'l4_cortex_exc.AfferentConnection.num_samples' : [45,90],
									     'l4_cortex_inh.L4InhL4InhConnection.base_weight' : [0.0002,0.0003],
									     'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0026,0.002,0.0015],
									     'l4_cortex_inh.ExcInhAfferentRatio' : [1.0,1.1],
									     'l4_cortex_exc.rand_struct_ratio' : [0.85],
									     'l23_cortex_exc.L4ExcL23ExcConnection.base_weight' : [0.0024],
									     'l4_cortex_exc.K' : [500,1000],
									     }).run_parameter_search()


if False:
    CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=64),{
									     'l4_cortex_inh.L4InhL4InhConnection.base_weight' : [0.0004],
									     'l4_cortex_inh.L4InhL4ExcConnection.base_weight' : [0.0007],
									     'l23_cortex_exc.L4ExcL23ExcConnection.base_weight' : [0.0024],
									     'l23_cortex_exc.L4ExcL23ExcConnection.num_samples' : [150],
									     'l23_cortex_inh.L4ExcL23InhConnection.base_weight' : [0.0016],
									     'l4_cortex_inh.ExcInhAfferentRatio' : [1.25],
									     'l4_cortex_exc.params.density' : [600],
									     'l4_cortex_exc.params.cell.params.tau_syn_E' : [3.0],
									     'l4_cortex_inh.L4InhL4ExcConnection.short_term_plasticity' : [None],
									     'l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f1.params.sigma' : [0.3],
									     'retina_lgn.params.gain_control.gain' : [70],
									     }).run_parameter_search()


if False:
    CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=64),{
									     'l23_cortex_exc.L4ExcL23ExcConnection.base_weight' : [0.0007],
									     'l23_cortex_exc.L4ExcL23ExcConnection.num_samples' : [150,300],
									     'l23_cortex_inh.L4ExcL23InhConnection.base_weight' : [0.0006,0.0005],
									     'l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f3.params.arborization_scaler' : [7.0,3.0],
									     'l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f1.params.sigma' : [0.3,1.0],
									     'retina_lgn.params.gain_control.gain' : [50],
									     'l4_cortex_inh.ExcInhAfferentRatio' : [1.0],
									     'retina_lgn.params.noise.stdev' : [3.4],
									     'l4_cortex_exc.rand_struct_ratio' : [0.85],
									     'l4_cortex_inh.L4InhL4InhConnection.base_weight' : [0.0004],
									     'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.00045],
									     'l4_cortex_inh.params.cell.params.tau_m' : [12.0],
									     }).run_parameter_search()

if False:
    CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=64),{
									     'l4_cortex_inh.L4InhL4InhConnection.base_weight' : [0.00067],
									     'l4_cortex_exc.inhibitory_connection_ratio' : [0.4,0.5,0.6],
									     'l4_cortex_exc.params.sx' : [3500],
									     'l4_cortex_inh.params.cell.params.tau_m' : [11.0,12.0,13.0],
									     'l4_cortex_inh.params.cell.params.v_rest' : [-77,-75,-72],
									     'retina_lgn.params.noise.stdev' : [3.4],
									     'l4_cortex_exc.rand_struct_ratio' : [0.85,0.9],
									     }).run_parameter_search()



if False:
    CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=32),{
									     'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0018,0.0019,0.002,0.0021],
									     'l4_cortex_exc.inhibitory_connection_ratio' : [0.65],
									     'l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f3.params.arborization_scaler' : [1.0],
#									     'l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f1.params.sigma' : [0.3],
									     'l4_cortex_exc.rand_struct_ratio' : [0.85],
									     'l4_cortex_inh.params.cell.params.tau_m' : [13.0,12.0],
#									     'l4_cortex_inh.params.cell.params.v_thresh' : [-52.0],
#									     'retina_lgn.params.gain_control.gain' : [50],
#									     'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.0014],
									     'l4_cortex_inh.L4InhL4InhConnection.base_weight' : [0.001],
									     'l4_cortex_exc.AfferentConnection.num_samples' : [40],
									     'l4_cortex_exc.params.sx' : [3500,4500],
									     'l4_cortex_exc.K' : [800],
									     'l4_cortex_exc.params.density' : [800],
#									     'retina_lgn.params.noise.stdev' : [2.7],
#									     'retina_lgn.params.gain_control.gain' : [50],
									     'l4_cortex_exc.l23_aff_ratio' : [0.4],
#									     'l4_cortex_exc.params.cell.params.v_rest' : [-70],
#									     'l4_cortex_inh.params.cell.params.v_rest' : [-70],
									     'l4_cortex_exc.AfferentConnection.short_term_plasticity.tau_rec' : [120],
									     'l4_cortex_exc.AfferentConnection.short_term_plasticity.U' : [0.5],
#									     'l4_cortex_exc.AfferentConnection.or_map_location' : ['./or_map_new','./or_map_new_cropped'],
#									     'l4_cortex_exc.params.cell.params.tau_syn_E' : [2.5],
#									     'l4_cortex_exc.params.cell.params.tau_syn_I' : [6.0],
#									     'l4_cortex_inh.ExcInhAfferentRatio' : [1.0],
#									     'l23_cortex_exc.L4ExcL23ExcConnection.weight_functions.f2.params.sigma' : [0.1],
#									     'retina_lgn.params.gain_control.non_linear_gain.contrast_scaler' : [0.1,0.05],
									     }).run_parameter_search()




if False:
    CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=32),{
									     'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0022,0.0024],
									     'l4_cortex_exc.inhibitory_connection_ratio' : [0.65],
									     'l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f3.params.arborization_scaler' : [1.0],
#									     'l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f1.params.sigma' : [0.3],
									     'l4_cortex_exc.rand_struct_ratio' : [0.9],
									     'l4_cortex_inh.params.cell.params.tau_m' : [11.0],
#									     'l4_cortex_inh.params.cell.params.v_thresh' : [-52.0],
#									     'retina_lgn.params.gain_control.gain' : [50],
									     'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.0014],
									     'l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.tau_rec' : [440,350,250],
#									     'l4_cortex_inh.L4InhL4InhConnection.base_weight' : [0.001],
#									     'l4_cortex_exc.AfferentConnection.num_samples' : [40],
									     'l4_cortex_exc.params.sx' : [4500],
#									     'l4_cortex_exc.K' : [800],
									     'l4_cortex_exc.params.density' : [800],
#									     'retina_lgn.params.noise.stdev' : [2.7],
#									     'retina_lgn.params.gain_control.gain' : [50],
									     'l4_cortex_exc.l23_aff_ratio' : [0.4],
#									     'l4_cortex_exc.params.cell.params.v_rest' : [-70],
#									     'l4_cortex_inh.params.cell.params.v_rest' : [-70],
#									     'l4_cortex_exc.AfferentConnection.short_term_plasticity.tau_rec' : [120],
									     'l4_cortex_exc.params.cell.params.v_rest' : [-70,-73,-76],
#									     'l4_cortex_exc.AfferentConnection.short_term_plasticity.U' : [0.5],
									     'retina_lgn.params.noise.stdev' : [2.9],
#									     'l4_cortex_exc.AfferentConnection.or_map_location' : ['./or_map_new','./or_map_new_cropped'],
#									     'l4_cortex_exc.params.cell.params.tau_syn_E' : [2.5],
#									     'l4_cortex_exc.params.cell.params.tau_syn_I' : [6.0],
#									     'l4_cortex_inh.ExcInhAfferentRatio' : [1.0],
#									     'l23_cortex_exc.L4ExcL23ExcConnection.weight_functions.f2.params.sigma' : [0.1],
#									     'retina_lgn.params.gain_control.non_linear_gain.contrast_scaler' : [0.1,0.05],
									     }).run_parameter_search()



if False:
    CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=32),{
									     'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0023,0.0024],
									     'l4_cortex_exc.inhibitory_connection_ratio' : [0.65],
									     'l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f3.params.arborization_scaler' : [4.0],
									     'l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f1.params.sigma' : [0.3],
									     'l4_cortex_exc.rand_struct_ratio' : [0.9],
									     'l4_cortex_inh.params.cell.params.tau_m' : [11.0],
#									     'l4_cortex_inh.params.cell.params.v_thresh' : [-52.0],
#									     'retina_lgn.params.gain_control.gain' : [50],
									     'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.0014],
									     'l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.tau_rec' : [100,200,400],
#									     'l4_cortex_inh.L4InhL4InhConnection.base_weight' : [0.001],
#									     'l4_cortex_exc.AfferentConnection.num_samples' : [40],
									     'l4_cortex_exc.params.sx' : [5000],
#									     'l4_cortex_exc.K' : [800],
									     'l4_cortex_exc.params.density' : [700,800],
#									     'retina_lgn.params.noise.stdev' : [2.7],
#									     'retina_lgn.params.gain_control.gain' : [50],
									     'l4_cortex_exc.l23_aff_ratio' : [0.4],
#									     'l4_cortex_exc.params.cell.params.v_rest' : [-70],
#									     'l4_cortex_inh.params.cell.params.v_rest' : [-70],
#									     'l4_cortex_exc.AfferentConnection.short_term_plasticity.tau_rec' : [120],
#									     'l4_cortex_exc.AfferentConnection.short_term_plasticity.U' : [0.5],
									     'retina_lgn.params.noise.stdev' : [2.9],
#									     'l4_cortex_exc.AfferentConnection.or_map_location' : ['./or_map_new','./or_map_new_cropped'],
#									     'l4_cortex_exc.params.cell.params.tau_syn_E' : [2.5],
#									     'l4_cortex_exc.params.cell.params.tau_syn_I' : [6.0],
#									     'l4_cortex_inh.ExcInhAfferentRatio' : [1.0],
#									     'l23_cortex_exc.L4ExcL23ExcConnection.weight_functions.f2.params.sigma' : [0.1],
#									     'retina_lgn.params.gain_control.non_linear_gain.contrast_scaler' : [0.1,0.05],
									     }).run_parameter_search()





if False:
    CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{
									     'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0015,0.002,0.0024],
									     'l4_cortex_inh.params.cell.params.tau_m' : [11.0],
									     'l4_cortex_exc.params.cell.params.tau_m' : [26],
									     'l4_cortex_exc.params.cell.params.cm' : [0.1],
									     'l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.tau_rec' : [100,200],
									     'l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.U' : [0.75],
									     'l4_cortex_exc.params.cell.params.tau_syn_E' : [2.5],
									     'l4_cortex_exc.params.cell.params.tau_syn_I' : [8.0,9,10],
									     'l4_cortex_exc.params.cell.params.v_rest' : [-75,-78],
									     'l4_cortex_exc.params.cell.params.delta_T' : [0.8],
									     'l4_cortex_exc.params.cell.params.e_rev_I' : [-80.0,-70.0],
									     'retina_lgn.params.noise.stdev' : [2.7,2.9],
									     'l4_cortex_exc.K' : [800,1600],
									     }).run_parameter_search()





if True:
    CombinationParameterSearch(SlurmSequentialBackend(num_threads=1,num_mpi=16),{
									     'l4_cortex_exc.L4ExcL4InhConnection.base_weight' : [0.0023],
									     'l4_cortex_exc.inhibitory_connection_ratio' : [0.65],
    								    	    'l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f3.params.arborization_scaler' : [4.0],
									     'l23_cortex_exc.L23ExcL23ExcConnection.weight_functions.f1.params.sigma' : [0.3],
									     'l4_cortex_exc.rand_struct_ratio' : [0.9],
									     'l4_cortex_inh.params.cell.params.tau_m' : [11.0],
									     'l4_cortex_exc.L4ExcL4ExcConnection.base_weight' : [0.0014],
									     'l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.tau_rec' : [100],
#									     'l4_cortex_exc.L4ExcL4ExcConnection.short_term_plasticity.U' : [0.75],
									     'l4_cortex_exc.params.density' : [700],
#
									     'l4_cortex_exc.params.sx' : [5000],
								     'l4_cortex_exc.l23_aff_ratio' : [0.4],
#									     'l4_cortex_exc.params.cell.params.tau_syn_I' : [6.0],
									     'retina_lgn.params.noise.stdev' : [2.9],
#									     'l4_cortex_exc.params.cell.params.v_rest' : [-70],
#									     'retina_lgn.params.gain_control.gain' : [70],
#									     'l23_cortex_exc.L23ExcL4ExcConnection.weight_functions.f1.params.arborization_constant' : [50,30,10],
#									     'l4_cortex_inh.L4InhL4ExcConnection.weight_functions.f1.params.sigma' : [3.0],
#									     'l4_cortex_inh.AfferentConnection.delay' : [1.1],
#									     'feedback' : [False],
#								     'l23_cortex_exc.L23ExcL4ExcConnection.num_samples' : [20,30],
#									     'l4_cortex_exc.l23_aff_ratio' : [0.4],
									     }).run_parameter_search()




