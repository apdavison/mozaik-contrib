"""


Usage:  python ffi-experiment.py param/defaults
"""

from pyNN import nest as sim
from model import TestModel
from mozaik.framework.experiment_controller import run_experiments, setup_experiments, Global
from mozaik.framework.experiment import MeasureNaturalImagesWithEyeMovement, \
                                        MeasureSpontaneousActivity, \
                                        MeasureOrientationTuningFullfield, \
                                        MeasureSizeTuning
from mozaik.visualization.plotting import OverviewPlot, RetinalInputMovie
from NeuroTools.parameters import ParameterSet
import mozaik
import os.path

logger = mozaik.getMozaikLogger("Mozaik")


params = setup_experiments('TestStimuli', sim)    
model = TestModel(sim, params)

experiment_list = [
    #Spontaneous Activity 
    #MeasureSpontaneousActivity(model, duration=148*7),

    #IMAGES WITH EYEMOVEMENT
    #MeasureNaturalImagesWithEyeMovement(model, stimulus_duration=200*7, num_trials=1)

    #GRATINGS
    MeasureOrientationTuningFullfield(model, num_orientations=1, spatial_frequency=0.8,
                                      temporal_frequency=2, grating_duration=148*7,
                                      contrasts=[0.5], num_trials=1),
    
    #SIZE TUNING
    #MeasureSizeTuning(model, max_size=4.0, orientation=0.0, spatial_frequency=0.8,
    #                  temporal_frequency=2, grating_duration=148*7, contrasts=[1.0],
    #                  num_trials=1, num_sizes=10),
]

data_store = run_experiments(model, experiment_list)

activity_plot_param = {
    'frame_rate': 5,  
    'bin_width': 50.0, 
    'scatter':  True,
    'resolution': 0
}       
    
OverviewPlot(data_store,
             ParameterSet({'sheet_name': 'X_ON', 'neuron': 0,
                           'sheet_activity': activity_plot_param}),
             plot_file_name=os.path.join(Global.root_directory, "overview_X_ON_0.png"),
             fig_param={'dpi': 100, 'figsize':(14,12)}).plot()
OverviewPlot(data_store,
             ParameterSet({'sheet_name': 'X_OFF', 'neuron': 0,
                           'sheet_activity': {}}),
             plot_file_name=os.path.join(Global.root_directory, "overview_X_OFF_0.png"),
             fig_param={'dpi': 100, 'figsize': (14,12)}).plot()
RetinalInputMovie(data_store,
                  ParameterSet({'frame_rate': 5})).plot()
    
import pylab
pylab.show()
