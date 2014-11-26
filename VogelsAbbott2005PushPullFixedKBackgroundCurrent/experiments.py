#!/usr/local/bin/ipython -i 
from mozaik.experiments import *
from mozaik.experiments.vision import *
from mozaik.sheets.population_selector import RCRandomPercentage
from parameters import ParameterSet
    
def create_experiments(model):

    return  [
                           #Lets kick the network up into activation
                           MeasureSpontaneousActivityWithPoissonStimulation(model,duration=2*2*50*7,sheet_list=["V1_Exc_L4","V1_Inh_L4"],drive_period=100.0,stimulation_configuration={'component' : 'mozaik.sheets.population_selector.RCRandomPercentage','params' : {'percentage' : 100.0}},lambda_list=[10000.0,10000.0],weight_list=[0.002,0.002]),
                           #Spontaneous Activity 
                           MeasureSpontaneousActivity(model,duration=145*7,num_trials=1),
                           #Orientation Tuning
                           #MeasureOrientationTuningFullfield(model,num_orientations=2,spatial_frequency=0.8,temporal_frequency=2,grating_duration=147*7,contrasts=[100],num_trials=2),
            ]
