import sys
sys.path.append('/home/jan/projects/mozaik/')
import matplotlib
from pyNN import nest as sim
from model import PushPullCCModel
from mozaik.framework.experiment_controller import run_experiments, setup_experiments, setup_logging
from mozaik.framework.experiment import *
from mozaik.visualization.plotting import *
from NeuroTools.parameters import ParameterSet
from mozaik.storage.queries import *
import mozaik

logger = mozaik.getMozaikLogger("Mozaik")


params = setup_experiments('FFI',sim)    
jens_model = PushPullCCModel(sim,params)
experiment_list =   [
                       #Spontaneous Activity 
                       #MeasureSpontaneousActivity(jens_model,duration=148*7),

                       #IMAGES WITH EYEMOVEMENT
                       MeasureNaturalImagesWithEyeMovement(jens_model,stimulus_duration=200*7,num_trials=1)

                       #GRATINGS
                       #MeasureOrientationTuningFullfield(jens_model,num_orientations=1,spatial_frequency=0.8,temporal_frequency=2,grating_duration=148*7,contrasts=[0.5],num_trials=1),
                       
                       #SIZE TUNING
                       #MeasureSizeTuning(jens_model,max_size=4.0,orientation=0.0,spatial_frequency=0.8,temporal_frequency=2,grating_duration=148*7,contrasts=[1.0],num_trials=1,num_sizes=10),
                    ]

data_store = run_experiments(jens_model,experiment_list)

activity_plot_param =    {
       'frame_rate' : 5,  
       'bin_width' : 50.0, 
       'scatter' :  True,
       'resolution' : 0
}       
    
OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : 0, 'sheet_activity' : activity_plot_param}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot()
OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : 0, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot()
RetinalInputMovie(data_store,ParameterSet({'frame_rate' : 5})).plot()
    
import pylab
pylab.show()
