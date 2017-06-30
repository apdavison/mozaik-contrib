from mozaik.experiments import *
from mozaik.experiments.vision import *
from mozaik.sheets.population_selector import RCRandomPercentage
from parameters import ParameterSet


def create_experiments_cortical_stimulation(model):

    return  [

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration' : 2*5*3*8*7})),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndHomogeneousOrientedStimulus(model,
                                                                        ParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 10,
                                                                                        'localstimulationarray_parameters' : ParameterSet({   
                                                                                                                            'size': 1200,
                                                                                                                            'spacing' : 50,
                                                                                                                            'itensity_fallof' : 30,
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function',
                                                                                                                            'stimulating_signal_parameters' : ParameterSet({
                                                                                                                                                                                'scale' : 0.02,
                                                                                                                                                                                'sigma' : 30,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 1.0,
                                                                                                                                                                                'duration' : 2*143*7,


                                                                                                                                }),
                                                                                                                            'current_update_interval' : 143*7,
                                                                                                                           })
                                                                                        }))
            ]

def create_experiments_short(model):

    return  [

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration' : 2*5*3*8*7})),

                # Measure orientation tuning with full-filed sinusoidal gratins
                MeasureOrientationTuningFullfield(model,ParameterSet({'num_orientations' : 2,
                                                                      'spatial_frequency' :0.8,
                                                                      'temporal_frequency' : 2,
                                                                      'grating_duration' : 2*143*7,
                                                                      'contrasts' : [100],
                                                                      'num_trials' : 10
                                                                      })),


            ]




def create_experiments(model):

    return  [

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration':8*2*5*3*8*7})),

                # Measure orientation tuning with full-filed sinusoidal gratins
                MeasureOrientationTuningFullfield(model,ParameterSet({'num_orientations' : 8,
                                                                      'spatial_frequency' :0.8,
                                                                      'temporal_frequency' : 2,
                                                                      'grating_duration' : 2*143*7,
                                                                      'contrasts' : [5,100],
                                                                      'num_trials' : 10
                                                                      })),
            ]


def create_experiments_old(model):

    return  [
                #Lets kick the network up into activation

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration':8*2*5*3*8*7})),

                # Measure response to natural image with simulated eye movement
                MeasureNaturalImagesWithEyeMovement(model,ParameterSet({'stimulus_duration':2*143*7,'num_trials':10})),

                # Measure orientation tuning with full-filed sinusoidal gratins
                MeasureOrientationTuningFullfield(model,ParameterSet({'num_orientations' : 8,'spatial_frequency' : 0.8,'temporal_frequency':2,'grating_duration':2*143*7,'contrasts':[5,100],'num_trials':10})),


            ]



def create_experiments_conn(model):
    
    return  [
                           #Spontaneous Activity 
                           NoStimulation(model,ParameterSet({'duration':100})),
            ]



