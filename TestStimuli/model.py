"""

"""

import sys
import numpy
from NeuroTools.parameters import ParameterSet
from mozaik.models.model import Model
from mozaik.framework import load_component
from mozaik.framework.space import VisualRegion

class TestModel(Model):
    
    required_parameters = ParameterSet({
        'retina_lgn': ParameterSet,
        'visual_field': ParameterSet 
    })
    
    def __init__(self, simulator, parameters):
        Model.__init__(self, simulator, parameters)        
        
        RetinaLGN = load_component(self.parameters.retina_lgn.component)
      
        # Build and instrument the network
        self.visual_field = VisualRegion(location_x=self.parameters.visual_field.centre[0],
                                         location_y=self.parameters.visual_field.centre[1],
                                         size_x=self.parameters.visual_field.size[0],
                                         size_y=self.parameters.visual_field.size[1])
        self.input_layer = RetinaLGN(self, self.parameters.retina_lgn.params)

        # which neurons to record
        tr = {
            'spikes': 'all', 
            'v': numpy.arange(0, 60, 1),
            'gsyn_exc': numpy.arange(0, 60, 1),
            'gsyn_inh': numpy.arange(0, 60, 1),
        }

        self.input_layer.sheets['X_ON'].to_record = tr #'all'
        self.input_layer.sheets['X_OFF'].to_record = tr #'all'
