#!/usr/local/bin/ipython -i 
from mozaik.experiments import *
from mozaik.experiments.vision import *
from mozaik.sheets.population_selector import RCRandomPercentage
from parameters import ParameterSet
    
def create_experiments(model):

    return  [
                #Lets kick the network up into activation
                PoissonNetworkKick(model,duration=8*8*7,drive_period=200.0,sheet_list=["V1_Exc_L4","V1_Inh_L4","V1_Exc_L2/3","V1_Inh_L2/3"],stimulation_configuration={'component' : 'mozaik.sheets.population_selector.RCRandomPercentage','params' : {'percentage' : 100.0}},lambda_list=[400.0,400.0,400.0,400.0],weight_list=[0.001,0.001,0.001,0.001]),
                #Spontaneous Activity 
                NoStimulation(model,duration=2*5*3*8*7),

                # Measure orientation tuning with full-filed sinusoidal gratins

                MeasureOrientationTuningFullfield(model,num_orientations=12,spatial_frequency=0.8,temporal_frequency=6,grating_duration=147*7,contrasts=[50,100],num_trials=15),

                # Measure response to natural image with simulated eye movement
                MeasureNaturalImagesWithEyeMovement(model,stimulus_duration=3*147*7,num_trials=15),

                #GRATINGS WITH EYEMOVEMENT
                MeasureDriftingSineGratingWithEyeMovement(model,spatial_frequency=0.8,temporal_frequency=6,stimulus_duration=3*147*7,num_trials=15,contrast=100),

            ]



def create_experiments_or(model):
    
    return  [
    
                           #Lets kick the network up into activation
                           PoissonNetworkKick(model,duration=8*8*7,drive_period=200.0,sheet_list=["V1_Exc_L4","V1_Inh_L4","V1_Exc_L2/3","V1_Inh_L2/3"],stimulation_configuration={'component' : 'mozaik.sheets.population_selector.RCRandomPercentage','params' : {'percentage' : 100.0}},lambda_list=[400.0,400.0,400.0,400.0],weight_list=[0.001,0.001,0.001,0.001]),
                           #Spontaneous Activity 
                           NoStimulation(model,duration=2*5*3*8*7),
                           # Measure orientation tuning with full-filed sinusoidal gratins
                           MeasureOrientationTuningFullfield(model,num_orientations=2,spatial_frequency=0.8,temporal_frequency=6,grating_duration=147*7,contrasts=[50,100],num_trials=4),
	       
                           # Measure response to natural image with simulated eye movement
                           MeasureNaturalImagesWithEyeMovement(model,stimulus_duration=147*7,num_trials=4),

            ]
