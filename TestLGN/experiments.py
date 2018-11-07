#!/usr/local/bin/ipython -i 
from mozaik.experiments import *
from mozaik.experiments.vision import *
from mozaik.sheets.population_selector import RCRandomPercentage
from parameters import ParameterSet
    

def create_experiments(model):

    return  [

                # #Spontaneous Activity 
                # MapPhaseResponseWithBarStimulus(model,ParameterSet({
                #                                                     'x' : 0,
                #                                                     'y' : 0,
                #                                                     'length' : 1/0.8/2.0 * 6.0,
                #                                                     'width' :  1/0.8/4.0,
                #                                                     'orientation' : 0,
                #                                                     'max_offset' : 1/0.8/2.0 * 3.0,
                #                                                     'steps' : 2,
                #                                                     'duration' : 1000,
                #                                                     'flash_duration' : 500, 
                #                                                     'relative_luminance' : 0,
                #                                                     'num_trials' : 1
                #                                                     })),
                                                                    
                # Measure orientation tuning with full-filed sinusoidal gratins
                

            ]



def create_experiments_bar(model):

    return  [

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration' : 2*5*3*8*7})),

                MapPhaseResponseWithBarStimulus(model,ParameterSet({
                                                                    'x' : 0,
                                                                    'y' : 0,
                                                                    'length' : 1/0.8/2.0 * 6.0,
                                                                    'width' :  1/0.8/4.0,
                                                                    'orientation' : 0,
                                                                    'max_offset' : 1/0.8/2.0 * 4.0,
                                                                    'steps' : 11,
                                                                    'duration' : 400,
                                                                    'flash_duration' : 400, 
                                                                    'relative_luminance' : 1.0,
                                                                    'num_trials' : 10
                                                                    })),

                MapPhaseResponseWithBarStimulus(model,ParameterSet({
                                                                    'x' : 0,
                                                                    'y' : 0,
                                                                    'length' : 1/0.8/2.0 * 6.0,
                                                                    'width' :  1/0.8/4.0,
                                                                    'orientation' : 0,
                                                                    'max_offset' : 1/0.8/2.0 * 4.0,
                                                                    'steps' : 11,
                                                                    'duration' : 400,
                                                                    'flash_duration' : 400, 
                                                                    'relative_luminance' : 0.0,
                                                                    'num_trials' : 10
                                                                    })),

            ]



def create_experiments_cs(model):

    return  [

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration' : 2*5*3*8*7})),

		MeasureContrastSensitivity(model,ParameterSet({'orientation' : 0,
                                                                      'spatial_frequency' :0.8,
                                                                      'temporal_frequency' : 2,
                                                                      'grating_duration' : 2*143*7,
                                                                      'contrasts' : [0.1,1,5,10,30,60,90,100],
                                                                      'num_trials' : 5
                                                                      })),
	    ]