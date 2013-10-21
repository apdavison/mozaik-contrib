#!/usr/local/bin/ipython -i 
from mozaik.experiments import *
from mozaik.experiments.vision import *
from mozaik.sheets.population_selector import RCRandomPercentage
from parameters import ParameterSet
    
def create_experiments(model):
    
    return  [
                           #Lets kick the network up into activation
                           PoissonNetworkKick(model,duration=8*8*7,drive_period=200.0,sheet_list=["V1_Exc_L4","V1_Inh_L4"],recording_configuration={'component' : 'mozaik.sheets.population_selector.RCRandomPercentage','params' : {'percentage' : 100.0}},lambda_list=[50.0,50.0],weight_list=[0.0006,0.0006]),
                           #Spontaneous Activity 
                           NoStimulation(model,duration=5*3*8*7),
                           # Measure orientation tuning with full-filed sinusoidal gratins
                           #MeasureOrientationTuningFullfield(model,num_orientations=4,spatial_frequency=0.8,temporal_frequency=2,grating_duration=147*7,contrasts=[100],num_trials=10),
                           MeasureOrientationTuningFullfield(model,num_orientations=6,spatial_frequency=0.8,temporal_frequency=2,grating_duration=147*7,contrasts=[30,100],num_trials=5),
            ]
