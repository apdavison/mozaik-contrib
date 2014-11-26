#!/usr/local/bin/ipython -i 
from mozaik.experiments import *
from mozaik.experiments.vision import *

from parameters import ParameterSet
    
def create_experiments(model):
    return             [
                           #OCTC
                           #MeasureOrientationContrastTuning(model,num_orientations=12,orientation=numpy.pi/2,center_radius=1.8,surround_radius=8,spatial_frequency=0.8,temporal_frequency=2,grating_duration=147*7,contrasts=[100],num_trials=2),

                           #Size Tuning  
                           #MeasureSizeTuning(model,num_sizes=15,max_size=4.5,orientation=numpy.pi/2,spatial_frequency=0.8,temporal_frequency=2,grating_duration=147*7,contrasts=[100],num_trials=5),
    
                           #Spontaneous Activity 
                           #MeasureSpontaneousActivity(model,duration=147*7,num_trials=8),
                    
                           #GRATINGS
                           MeasureOrientationTuningFullfield(model,num_orientations=1,spatial_frequency=0.8,temporal_frequency=2,grating_duration=147*7,contrasts=[10,20,30,40,50,60,70,80,90,100],num_trials=15),

                           #IMAGES WITH EYEMOVEMENT
                           #MeasureNaturalImagesWithEyeMovement(model,stimulus_duration=147*7,num_trials=10),

                           #GRATINGS WITH EYEMOVEMENT
                           #MeasureDriftingSineGratingWithEyeMovement(model,spatial_frequency=0.8,temporal_frequency=2,stimulus_duration=147*7,num_trials=10,contrast=100),
                       ]
