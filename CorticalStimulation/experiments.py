from mozaik.experiments import *
from mozaik.experiments.vision import *
from mozaik.sheets.population_selector import RCRandomPercentage
from parameters import ParameterSet
from mozaik.tools.distribution_parametrization import MozaikExtendedParameterSet


def create_experiments_tmp(model):

    return  [
                # Measure orientation tuning with full-filed sinusoidal gratins
                MeasureOrientationTuningFullfieldA(model,ParameterSet({'num_orientations' : 2,'spatial_frequency' : 0.8,'temporal_frequency':2,'grating_duration':400,'contrasts':[100],'num_trials':10,'onset_time' : 100, 'offset_time': 300})),

	    ]

def create_experiments_cortical_stimulation_prof(model):

    return  [

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration' : 2*5*3*8*7})),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 10,
											'num_orientations' : 1,
										        'intensities' : [0.01],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3500,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
                                                                                                                                                                                'scale' : 1,
                                                                                                                                                                                'sigma' : 30,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 1.0,
                                                                                                                                                                                'duration' : 400,
																						'onset_time' : 100,
																						'offset_time' : 300,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        }))
	    ]

def create_experiments_cortical_stimulation_luminance(model):
    return  [

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration' : 3*8*7})),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 1,
											'num_orientations' : 1,
										        'intensities' : [0.03125],#[0.0078125,0.015625,0.03125,0.0625,0.125,0.25,0.5,1],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3500,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
                                                                                                                                                                                'scale' : 1,
                                                                                                                                                                                'sigma' : 30,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 400,
																						'onset_time' : 100,
																						'offset_time' : 300,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),






	    ]



def create_experiments_cortical_stimulation_exc(model):

    return  [

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration' : 2*5*3*8*7})),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 10,
											'num_orientations' : 8,
										        'intensities' : [0.0078125,0.015625,0.125,1.0],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3500,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
                                                                                                                                                                                'scale' : 1,
                                                                                                                                                                                'sigma' : 30,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 400,
																						'onset_time' : 100,
																						'offset_time' : 300,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 10,
											'num_orientations' : 1,
										        'intensities' : [0.0078125,0.015625,0.03125,0.0625,0.125,0.25,0.5,1],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3500,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
                                                                                                                                                                                'scale' : 1,
                                                                                                                                                                                'sigma' : 30,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 400,
																						'onset_time' : 100,
																						'offset_time' : 300,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),






	    ]


def create_experiments_cortical_stimulation_excinh(model):

    return  [

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration' : 2*5*3*8*7})),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3','V1_Inh_L2/3'],
                                                                                        'num_trials' : 10,
											'num_orientations' : 8,
										        'intensities' : [0.0078125,0.015625,0.125,1.0],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3500,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
                                                                                                                                                                                'scale' : 1,
                                                                                                                                                                                'sigma' : 30,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 400,
																						'onset_time' : 100,
																						'offset_time' : 300,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3','V1_Inh_L2/3'],
                                                                                        'num_trials' : 10,
											'num_orientations' : 1,
										        'intensities' : [0.0078125,0.015625,0.03125,0.0625,0.125,0.25,0.5,1],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3500,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
                                                                                                                                                                                'scale' : 1,
                                                                                                                                                                                'sigma' : 30,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 400,
																						'onset_time' : 100,
																						'offset_time' : 300,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),






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



            ]



def create_experiments_conn(model):
    
    return  [
                           #Spontaneous Activity 
                           NoStimulation(model,ParameterSet({'duration':100})),
            ]



