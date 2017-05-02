#!/usr/local/bin/ipython -i 
from mozaik.experiments import *
from mozaik.experiments.vision import *
from mozaik.sheets.population_selector import RCRandomPercentage
from parameters import ParameterSet
    
def create_experiments(model):

    return  [
                #Lets kick the network up into activation

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration':8*2*5*3*8*7})),

                # Measure orientation tuning with full-filed sinusoidal gratins
                MeasureOrientationTuningFullfield(model,ParameterSet({'num_orientations' : 8,'spatial_frequency' : 0.8,'temporal_frequency':2,'grating_duration':2*143*7,'contrasts':[5,100],'num_trials':10})),

                # Measure response to natural image with simulated eye movement
                MeasureNaturalImagesWithEyeMovement(model,ParameterSet({'stimulus_duration':2*143*7,'num_trials':10})),

                #GRATINGS WITH EYEMOVEMENT
                #MeasureDriftingSineGratingWithEyeMovement(model,spatial_frequency=0.8,temporal_frequency=6,stimulus_duration=3*143*7,num_trials=15,contrast=100),

            ]


def create_experiments_spont(model):
    
    return  [
                           #Spontaneous Activity 
                            NoStimulation(model,ParameterSet({'duration':8*2*5*3*8*7})),
    ]


def create_experiments_or(model):
    
    return  [
                           #Spontaneous Activity 
                           NoStimulation(model,ParameterSet({'duration':2*5*3*8*7})),
                           # Measure orientation tuning with full-filed sinusoidal gratins
                           MeasureOrientationTuningFullfield(model,ParameterSet({'num_orientations':2,'spatial_frequency':0.8,'temporal_frequency':2,'grating_duration':2*143*7,'contrasts':[5,100],'num_trials':5})),
	       
                           # Measure response to natural image with simulated eye movement
                           MeasureNaturalImagesWithEyeMovement(model,ParameterSet({'stimulus_duration':2*143*7,'num_trials' : 5})),

            ]

def create_experiments_stc(model):
    
    return  [
    
                           #Spontaneous Activity 
                           NoStimulation(model,ParameterSet({'duration':2*5*3*8*7})),
 
                           #Size Tuning  
                           MeasureSizeTuning(model,ParameterSet({'num_sizes':12,'max_size':3.0,'log_spacing' : True,'orientation' : 0,'spatial_frequency' : 0.8,'temporal_frequency' : 2,'grating_duration' : 5*2*143*7,'contrasts' : [5,100],'num_trials' :1,'with_flat': False})),
            ]


def create_experiments_octc(model):
    
    return  [
    
                           #Spontaneous Activity 
                           NoStimulation(model,ParameterSet({'duration':2*5*3*8*7})),
 
                           #OCTC
                           MeasureOrientationContrastTuning(model,ParameterSet({'num_orientations' : 8,'orientation' : 0,'center_radius' : 0.5,'surround_radius' :20.0,'spatial_frequency' : 0.8,'temporal_frequency' : 2,'grating_duration' : 2*143*7,'contrasts' : [100],'num_trials' : 20})),
            ]


def create_experiments_conn(model):
    
    return  [
                           #Spontaneous Activity 
                           NoStimulation(model,ParameterSet({'duration':100})),
            ]



