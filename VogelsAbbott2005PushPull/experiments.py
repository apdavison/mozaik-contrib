#!/usr/local/bin/ipython -i 
from mozaik.experiments import *
from mozaik.experiments.vision import *
from mozaik.sheets.population_selector import RCRandomPercentage
from parameters import ParameterSet
    
def create_experiments(model):

    return  [
                           #Lets kick the network up into activation
                           PoissonNetworkKick(model,duration=8*7,sheet_list=["V1_Exc_L4","V1_Inh_L4"],recording_configuration={'component' : 'mozaik.sheets.population_selector.RCRandomPercentage','params' : {'percentage' : 20.0}},lambda_list=[100.0,100.0],weight_list=[0.1,0.1]),
                           #Spontaneous Activity 
                           MeasureSpontaneousActivity(model,duration=145*7,num_trials=1),
            ]
