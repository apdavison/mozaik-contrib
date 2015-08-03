#!/usr/local/bin/ipython -i 
from mozaik.experiments import *
from mozaik.sheets.population_selector import RCRandomPercentage
from parameters import ParameterSet
    
def create_experiments(model):
    
    return  [
                           #Lets kick the network up into activation
                           PoissonNetworkKick(model,duration=8*8*7,drive_period=200.0,sheet_list=["Exc1","Exc2","Inh1","Inh2"],stimulation_configuration={'component' : 'mozaik.sheets.population_selector.RCRandomPercentage','params' : {'percentage' : 100.0}},lambda_list=[400.0,400.0,400.0,400.0],weight_list=[0.001,0.001,0.001,0.001]),
                           #Spontaneous Activity 
                           #NoStimulation(model,duration=285*7*50),
			   NoStimulation(model,duration=15*7*50),

            ]
