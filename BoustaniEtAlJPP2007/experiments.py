#!/usr/local/bin/ipython -i 
from mozaik.experiments import *
from mozaik.framework.population_selector import RCRandomPercentage
from parameters import ParameterSet
    
def create_experiments(model):
    l4exc_kick = RCRandomPercentage(model.sheets["V1_Exc_L4"],ParameterSet({'percentage': 10.0}))
    l4inh_kick = RCRandomPercentage(model.sheets["V1_Inh_L4"],ParameterSet({'percentage': 10.0}))

    return  [
                           #Lets kick the network up into activation
                           PoissonNetworkKick(model,duration=70*7,sheet_list=["V1_Exc_L4","V1_Inh_L4"],recording_configuration_list=[l4exc_kick,l4inh_kick],lambda_list=[50,50]),
                           #Spontaneous Activity 
                           NoStimulation(model,duration=145*7),
            ]
