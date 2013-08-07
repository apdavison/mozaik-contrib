#!/usr/local/bin/ipython -i 
from mozaik.experiments import *
from mozaik.sheets.population_selector import RCRandomPercentage
from parameters import ParameterSet
    
def create_experiments(model):
    
    return  [
                           #Lets kick the network up into activation
                           PoissonNetworkKick(model,duration=10*8*7,drive_period=200.0,sheet_list=["V1_Exc_L4","V1_Inh_L4"],recording_configuration={'component' : 'mozaik.sheets.population_selector.RCRandomPercentage','params' : {'percentage' : 100.0}},lambda_list=[10000.0,10000.0],weight_list=[0.006,0.006]),
                           #Spontaneous Activity 
                           NoStimulation(model,duration=3*8*7),
            ]
