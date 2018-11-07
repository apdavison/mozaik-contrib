from mozaik.experiments import *
from mozaik.experiments.vision import *
from mozaik.sheets.population_selector import RCRandomPercentage
from parameters import ParameterSet
from mozaik.tools.distribution_parametrization import MozaikExtendedParameterSet


def create_experiments_tmp(model):

    return  [
                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration':8*2*5*3*8*7})),

                # Measure orientation tuning with full-filed sinusoidal gratins
                MeasureOrientationTuningFullfieldA(model,ParameterSet({'num_orientations' : 8,'spatial_frequency' : 0.8,'temporal_frequency':2,'grating_duration':1000,'contrasts':[3,5,10,100],'num_trials':10,'onset_time' : 200, 'offset_time': 800})),

		MeasureOrientationTuningFullfieldA(model,ParameterSet({'num_orientations' : 1,'spatial_frequency' : 0.8,'temporal_frequency':2,'grating_duration':1000,'contrasts': [0,20,30,40,50,60,70,80,90],'num_trials':10,'onset_time' : 200, 'offset_time': 800})),
	    ]

def create_experiments_contr(model):

    return  [

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration':8*2*5*3*8*7})),

                # Measure orientation tuning with full-filed sinusoidal gratins
                MeasureContrastSensitivityA(model,ParameterSet({'orientation' : 0,
                                                                      'spatial_frequency' :0.8,
                                                                      'temporal_frequency' : 2,
                                                                      'grating_duration' : 1000,
                                                                      'contrasts' : [0,3,5,10,20,30,40,50,60,70,80,90,100],
                                                                      'num_trials' : 10,
								      'onset_time' : 100, 
								      'offset_time': 800
                                                                      })),
            ]


def create_experiments_cortical_stimulation_prof(model):

    return  [

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration' : 2*5*3*8*7})),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 1,
											'num_orientations' : 1,
										        'intensities' : [0.125,0.5],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd10.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
                                                                                                                                                                                'scale' : 1,
                                                                                                                                                                                'sigma' : 30,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 600,
																						'onset_time' : 100,
																						'offset_time' : 500,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),


                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 1,
											'num_orientations' : 1,
										        'intensities' : [0.125,0.25],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd10.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
                                                                                                                                                                                'scale' : 1,
                                                                                                                                                                                'sigma' : 30,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 600,
																						'onset_time' : 100,
																						'offset_time' : 500,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        }))

	    ]

def create_experiments_cortical_stimulation_luminance_exc(model):
    return  [

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration' : 2*5*3*8*7})),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 10,
											'num_orientations' : 1,
										        'intensities' : [0.005,0.01,0.02,0.04,0.08,0.16,0.32,0.64],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd10.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
                                                                                                                                                                                'scale' : 1,
                                                                                                                                                                                'sigma' : 30,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 1000,
																						'onset_time' : 200,
																						'offset_time' : 800,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),


	    ]


def create_experiments_cortical_stimulation_luminance_excinh(model):
    return  [

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration' : 2*5*3*8*7})),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3','V1_Inh_L2/3'],
                                                                                        'num_trials' : 10,
											'num_orientations' : 1,
										        'intensities' : [0.005,0.01,0.02,0.04,0.08,0.16,0.32,0.64],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd10.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
                                                                                                                                                                                'scale' : 1,
                                                                                                                                                                                'sigma' : 30,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 1000,
																						'onset_time' : 200,
																						'offset_time' : 800,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),
	    ]



def create_experiments_cortical_stimulation_or_exc(model):

    return  [

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration' : 2*5*3*8*7})),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol_ContrastBased(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 10,
											'num_orientations' : 8,
										        'contrasts' : [3,5,10,100],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd10.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function_Naka',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
																						'contrast' : 0,
																						'nv_r_max' : 9.099,
																						'nv_c50' : 1.09,
																						'cs_r_max' : 82.77,
																						'cs_c50' : 0.0505,
																						'cs_offset' : 0,
		                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 1000,
																						'onset_time' : 100,
																						'offset_time' : 800,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol_ContrastBased(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 10,
											'num_orientations' : 1,
										        'contrasts' : [0,20,30,40,50,60,70,80,90],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd10.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function_Naka',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
																						'contrast' : 0,
																						'nv_r_max' : 9.099,
																						'nv_c50' : 1.09,
																						'cs_r_max' : 82.77,
																						'cs_c50' : 0.0505,
																						'cs_offset' : 0,
		                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 1000,
																						'onset_time' : 100,
																						'offset_time' : 800,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),


	    ]


def create_experiments_cortical_stimulation_or_excinh(model):

    return  [

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration' : 2*5*3*8*7})),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol_ContrastBased(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3','V1_Inh_L2/3'],
                                                                                        'num_trials' : 10,
											'num_orientations' : 8,
										        'contrasts' : [3,5,10,100],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd10.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function_Naka',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
																						'contrast' : 0,
																						'nv_r_max' : 9.099,
																						'nv_c50' : 1.09,
																						'cs_r_max' :  14000,
																						'cs_c50' : 57.85,
																						'cs_offset' : 0,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 1000,
																						'onset_time' : 100,
																						'offset_time' : 800,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol_ContrastBased(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3','V1_Inh_L2/3'],
                                                                                        'num_trials' : 10,
											'num_orientations' : 1,
										        'contrasts' : [0,20,30,40,50,60,70,80,90],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd10.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function_Naka',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
																						'contrast' : 0,
																						'nv_r_max' : 9.099,
																						'nv_c50' : 1.09,
																						'cs_r_max' :  14000,
																						'cs_c50' : 57.85,
																						'cs_offset' : 0,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 1000,
																						'onset_time' : 100,
																						'offset_time' : 800,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),


	    ]

def create_experiments_cortical_stimulation_or_nolat(model):

    return  [

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration' : 2*5*3*8*7})),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol_ContrastBased(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 1,
											'num_orientations' : 16,
										        'contrasts' : [1,3,5,10,100],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd10.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function_Naka',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
																						'contrast' : 0,
																						'nv_r_max' : 9.099,
																						'nv_c50' : 1.09,
																						'cs_r_max' :  70,
																						'cs_c50' : 0.0042,
																						'cs_exponent' : 2.8,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 1000,
																						'onset_time' : 100,
																						'offset_time' : 800,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol_ContrastBased(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 1,
											'num_orientations' : 1,
										        'contrasts' : [0,20,30,40,50,60,70,80,90],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd10.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function_Naka',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
																						'contrast' : 0,
																						'nv_r_max' : 9.099,
																						'nv_c50' : 1.09,
																						'cs_r_max' :  70,
																						'cs_c50' : 0.0042,
																						'cs_exponent' : 2.8,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 1000,
																						'onset_time' : 100,
																						'offset_time' : 800,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),


	    ]









def create_experiments_cortical_stimulation_sharpness(model):

    return  [

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration' : 2*5*3*8*7})),


                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol_ContrastBased(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 5,
											'num_orientations' : 8,
										        'contrasts' : [3,100],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd10.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function_Naka',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
																						'contrast' : 0,
																						'nv_r_max' : 9.099,
																						'nv_c50' : 1.09,
																						'cs_r_max' : 82.77,
																						'cs_c50' : 0.0505,
																						'cs_exponent' : 1.0,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.1,
                                                                                                                                                                                'duration' : 1000,
																						'onset_time' : 100,
																						'offset_time' : 800,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),


                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol_ContrastBased(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 5,
											'num_orientations' : 8,
										        'contrasts' : [3,100],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd10.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function_Naka',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
																						'contrast' : 0,
																						'nv_r_max' : 9.099,
																						'nv_c50' : 1.09,
																						'cs_r_max' : 82.77,
																						'cs_c50' : 0.0505,
																						'cs_exponent' : 1.0,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.3,
                                                                                                                                                                                'duration' : 1000,
																						'onset_time' : 100,
																						'offset_time' : 800,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),


                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol_ContrastBased(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 5,
											'num_orientations' : 8,
										        'contrasts' : [3,100],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd10.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function_Naka',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
																						'contrast' : 0,
																						'nv_r_max' : 9.099,
																						'nv_c50' : 1.09,
																						'cs_r_max' : 82.77,
																						'cs_c50' : 0.0505,
																						'cs_exponent' : 1.0,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 1000,
																						'onset_time' : 100,
																						'offset_time' : 800,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),


                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol_ContrastBased(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 5,
											'num_orientations' : 8,
										        'contrasts' : [3,100],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd10.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function_Naka',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
																						'contrast' : 0,
																						'nv_r_max' : 9.099,
																						'nv_c50' : 1.09,
																						'cs_r_max' : 82.77,
																						'cs_c50' : 0.0505,
																						'cs_exponent' : 1.0,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 1.0,
                                                                                                                                                                                'duration' : 1000,
																						'onset_time' : 100,
																						'offset_time' : 800,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol_ContrastBased(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 5,
											'num_orientations' : 8,
										        'contrasts' : [3,100],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd10.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function_Naka',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
																						'contrast' : 0,
																						'nv_r_max' : 9.099,
																						'nv_c50' : 1.09,
																						'cs_r_max' : 82.77,
																						'cs_c50' : 0.0505,

																						'cs_exponent' : 1.0,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 2.0,
                                                                                                                                                                                'duration' : 1000,
																						'onset_time' : 100,
																						'offset_time' : 800,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol_ContrastBased(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 5,
											'num_orientations' : 8,
										        'contrasts' : [3,100],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd10.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function_Naka',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
																						'contrast' : 0,
																						'nv_r_max' : 9.099,
																						'nv_c50' : 1.09,
																						'cs_r_max' : 82.77,
																						'cs_c50' : 0.0505,
																						'cs_exponent' : 1.0,
                                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 3.0,
                                                                                                                                                                                'duration' : 1000,
																						'onset_time' : 100,
																						'offset_time' : 800,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),


	    ]


def create_experiments_cortical_stimulation_lee_size(model):

    return  [

                #Spontaneous Activity 
                NoStimulation(model,ParameterSet({'duration' : 2*5*3*8*7})),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol_ContrastBased(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 5,
											'num_orientations' : 8,
										        'contrasts' : [3,100],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 10,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd10.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function_Naka',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
																						'contrast' : 0,
																						'nv_r_max' : 9.099,
																						'nv_c50' : 1.09,
																						'cs_r_max' : 82.77,
																						'cs_c50' : 0.0505,
																						'cs_exponent' : 1.0,
		                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 1000,
																						'onset_time' : 100,
																						'offset_time' : 800,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),


                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol_ContrastBased(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 5,
											'num_orientations' : 8,
										        'contrasts' : [3,100],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 20,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd20.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function_Naka',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
																						'contrast' : 0,
																						'nv_r_max' : 9.099,
																						'nv_c50' : 1.09,
																						'cs_r_max' : 82.77,
																						'cs_c50' : 0.0505,
																						'cs_exponent' : 1.0,
		                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 1000,
																						'onset_time' : 100,
																						'offset_time' : 800,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol_ContrastBased(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 5,
											'num_orientations' : 8,
										        'contrasts' : [3,100],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 50,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd50.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function_Naka',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
																						'contrast' : 0,
																						'nv_r_max' : 9.099,
																						'nv_c50' : 1.09,
																						'cs_r_max' : 82.77,
																						'cs_c50' : 0.0505,
																						'cs_exponent' : 1.0,

		                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 1000,
																						'onset_time' : 100,
																						'offset_time' : 800,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol_ContrastBased(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 5,
											'num_orientations' : 8,
										        'contrasts' : [3,100],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 100,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd100.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function_Naka',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
																						'contrast' : 0,
																						'nv_r_max' : 9.099,
																						'nv_c50' : 1.09,
																						'cs_r_max' : 82.77,
																						'cs_c50' : 0.0505,
																						'cs_exponent' : 1.0,
		                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 1000,
																						'onset_time' : 100,
																						'offset_time' : 800,
                                                                                                                                }),
                                                                                                                            'current_update_interval' : 1,
                                                                                                                           })
                                                                                        })),

                # Measure orientation tuning with full-filed sinusoidal gratins
                CorticalStimulationWithStimulatorArrayAndOrientationTuningProtocol_ContrastBased(model,
                                                                        MozaikExtendedParameterSet({   
                                                                                        'sheet_list' : ['V1_Exc_L2/3'],
                                                                                        'num_trials' : 5,
											'num_orientations' : 8,
										        'contrasts' : [3,100],
                                                                                        'localstimulationarray_parameters' : MozaikExtendedParameterSet({   
                                                                                                                            'size': 3600,
                                                                                                                            'spacing' : 300,
                                                                                                                            'itensity_fallof' : 30,
															    'depth_sampling_step' : 10,
															    'light_source_light_propagation_data' : 'light_scattering_radial_profiles_lsd300.pickle',
                                                                                                                            'stimulating_signal' : 'mozaik.sheets.direct_stimulator.test_stimulating_function_Naka',
                                                                                                                            'stimulating_signal_parameters' : MozaikExtendedParameterSet({
																						'contrast' : 0,
																						'nv_r_max' : 9.099,
																						'nv_c50' : 1.09,
																						'cs_r_max' : 82.77,
																						'cs_c50' : 0.0505,
																						'cs_exponent' : 1.0,
		                                                                                                                                                                'orientation' : 0,
                                                                                                                                                                                'sharpness' : 0.5,
                                                                                                                                                                                'duration' : 1000,
																						'onset_time' : 100,
																						'offset_time' : 800,
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



