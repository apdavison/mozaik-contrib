#!/usr/local/bin/ipython -i 
from mozaik.experiments import *
from mozaik.experiments.vision import *
from mozaik.sheets.population_selector import RCRandomPercentage
from parameters import ParameterSet
    
def create_experiments(model):
    l4exc_kick = RCRandomPercentage(model.sheets["V1_Exc_L4"],ParameterSet({'percentage': 20.0}))
    l4inh_kick = RCRandomPercentage(model.sheets["V1_Inh_L4"],ParameterSet({'percentage': 20.0}))

    return  [
                           #Lets kick the network up into activation
                           MeasureSpontaneousActivityWithPoissonStimulation(model,duration=10*7,sheet_list=["V1_Exc_L4","V1_Inh_L4"],recording_configuration_list=[l4exc_kick,l4inh_kick],lambda_list=[100,100]),
                           #Spontaneous Activity 
                           MeasureSpontaneousActivity(model,duration=145*7,num_trials=1),
                           #Orientation Tuning
                           MeasureOrientationTuningFullfield(model,num_orientations=2,spatial_frequency=0.8,temporal_frequency=2,grating_duration=147*7,contrasts=[100],num_trials=2),
            ]
