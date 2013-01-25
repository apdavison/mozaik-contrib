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
from mozaik.storage.datastore import PickledDataStore
import mozaik

logger = mozaik.getMozaikLogger("Mozaik")

if True:
    params = setup_experiments('FFI',sim)    
    jens_model = PushPullCCModel(sim,params)
    experiment_list =   [
                       #Spontaneous Activity 
                       #MeasureSpontaneousActivity(jens_model,duration=50*7,num_trials=1),

                       #IMAGES WITH EYEMOVEMENT
                       MeasureNaturalImagesWithEyeMovement(jens_model,stimulus_duration=200*7,num_trials=1),

                       #GRATINGS
                       MeasureOrientationTuningFullfield(jens_model,num_orientations=1,spatial_frequency=0.8,temporal_frequency=2,grating_duration=50*7,contrasts=[100],num_trials=1),
                       
                       #GRATINGS WITH EYEMOVEMENT
                       MeasureDriftingSineGratingWithEyeMovement(jens_model,spatial_frequency=0.8,temporal_frequency=2,stimulus_duration=200*7,num_trials=1,contrast=100),
                       
                       #SIZE TUNING
                       #MeasureSizeTuning(jens_model,max_size=1.0,orientation=0.0,spatial_frequency=0.8,temporal_frequency=2,grating_duration=50*7,contrasts=[20,50,100],num_trials=1,num_sizes=3),
                    ]

    data_store = run_experiments(jens_model,experiment_list)
    data_store.save()
else:
    data_store = PickledDataStore(load=True,parameters=ParameterSet({'root_directory':'A'}))

activity_plot_param =    {
       'frame_rate' : 5,  
       'bin_width' : 5.0, 
       'scatter' :  True,
       'resolution' : 0
}       
OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_ON', 'neuron' : 0, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot()
OverviewPlot(data_store,ParameterSet({'sheet_name' : 'X_OFF', 'neuron' : 0, 'sheet_activity' : {}}),fig_param={'dpi' : 100,'figsize': (14,12)}).plot()
RetinalInputMovie(data_store,ParameterSet({'frame_rate' : 5})).plot()
    
import pylab
pylab.show()
