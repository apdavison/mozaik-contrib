#!/usr/local/bin/ipython -i 
from mozaik.experiments import *
from mozaik.sheets.population_selector import RCRandomPercentage
from parameters import ParameterSet
    
def create_experiments(model):
    
    return  [
                           #Spontaneous Activity 
                           NoStimulation(model,duration=1000*7),
            ]
