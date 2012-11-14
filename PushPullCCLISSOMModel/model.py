import sys
sys.path.append('/home/jan/projects/mozaik/')
import numpy
from NeuroTools.parameters import ParameterSet
from mozaik.models.model import Model
from mozaik.connectors.meta_connectors import GaborConnector
from mozaik.connectors.modular_connectors import ModularProbabilisticConnector
from mozaik.framework import load_component
from mozaik.framework.space import VisualRegion

class PushPullCCModel(Model):
    
    required_parameters = ParameterSet({
        'l4_cortex_exc' : ParameterSet, 
        'l4_cortex_inh' : ParameterSet, 
        'l23_cortex_exc' : ParameterSet, 
        'l23_cortex_inh' : ParameterSet, 
        'retina_lgn' : ParameterSet ,
        'visual_field' : ParameterSet 
    })
    
    def __init__(self,simulator,parameters):
        Model.__init__(self,simulator,parameters)        
        # Load components
        CortexExcL4 = load_component(self.parameters.l4_cortex_exc.component)
        CortexInhL4 = load_component(self.parameters.l4_cortex_inh.component)
        #CortexExcL23 = load_component(self.parameters.l23_cortex_exc.component)
        #CortexInhL23 = load_component(self.parameters.l23_cortex_inh.component)
        
        RetinaLGN = load_component(self.parameters.retina_lgn.component)
      
        # Build and instrument the network
        self.visual_field = VisualRegion(location_x=self.parameters.visual_field.centre[0],location_y=self.parameters.visual_field.centre[1],size_x=self.parameters.visual_field.size[0],size_y=self.parameters.visual_field.size[1])
        self.input_layer = RetinaLGN(self, self.parameters.retina_lgn.params)
        cortex_exc_l4 = CortexExcL4(self, self.parameters.l4_cortex_exc.params)
        cortex_inh_l4 = CortexInhL4(self, self.parameters.l4_cortex_inh.params)
        #cortex_exc_l23 = CortexExcL23(self, self.parameters.l23_cortex_exc.params)
        #cortex_inh_l23 = CortexInhL23(self, self.parameters.l23_cortex_inh.params)

        # which neurons to record
        tr = {'spikes' : 'all', 
              'v' : numpy.arange(0,21,1),
              'gsyn_exc' :numpy.arange(0,21,1),
              'gsyn_inh' : numpy.arange(0,21,1),
        }
        
        cortex_exc_l4.to_record = tr #'all'
        cortex_inh_l4.to_record = tr #'all'
        #cortex_exc_l23.to_record = tr #'all'
        #cortex_inh_l23.to_record = tr #'all'
        self.input_layer.sheets['X_ON'].to_record = tr #'all'
        self.input_layer.sheets['X_OFF'].to_record = tr #'all'

        # initialize projections
        GaborConnector(self,self.input_layer.sheets['X_ON'],self.input_layer.sheets['X_OFF'],cortex_exc_l4,self.parameters.l4_cortex_exc.AfferentConnection,'V1AffConnection')
        GaborConnector(self,self.input_layer.sheets['X_ON'],self.input_layer.sheets['X_OFF'],cortex_inh_l4,self.parameters.l4_cortex_inh.AfferentConnection,'V1AffInhConnection')
        
        ModularProbabilisticConnector(self,'V1L4ExcL4ExcConnection',cortex_exc_l4,cortex_exc_l4,self.parameters.l4_cortex_exc.L4ExcL4ExcConnection).connect()
        ModularProbabilisticConnector(self,'V1L4ExcL4InhConnection',cortex_exc_l4,cortex_inh_l4,self.parameters.l4_cortex_exc.L4ExcL4InhConnection).connect()
        ModularProbabilisticConnector(self,'V1L4InhL4ExcConnection',cortex_inh_l4,cortex_exc_l4,self.parameters.l4_cortex_inh.L4InhL4ExcConnection).connect()
        ModularProbabilisticConnector(self,'V1L4InhL4InhConnection',cortex_inh_l4,cortex_inh_l4,self.parameters.l4_cortex_inh.L4InhL4InhConnection).connect()
 
        #ModularProbabilisticConnector(self,'V1ExcL23ExcL23Connection',cortex_exc_l23,cortex_exc_l23,self.parameters.l23_cortex_exc.L23ExcL23ExcConnection).connect()
        #ModularProbabilisticConnector(self,'V1ExcL23InhL23Connection',cortex_exc_l23,cortex_inh_l23,self.parameters.l23_cortex_exc.L23ExcL23InhConnection).connect()
        #ModularProbabilisticConnector(self,'V1InhL23ExcL23Connection',cortex_inh_l23,cortex_exc_l23,self.parameters.l23_cortex_inh.L23InhL23ExcConnection).connect()
        #ModularProbabilisticConnector(self,'V1InhL23InhL23Connection',cortex_inh_l23,cortex_inh_l23,self.parameters.l23_cortex_inh.L23InhL23InhConnection).connect()
        #ModularProbabilisticConnector(self,'V1ExcL4ExcL23Connection',cortex_exc_l4,cortex_exc_l23,self.parameters.l23_cortex_exc.L4ExcL23ExcConnection).connect()
        #ModularProbabilisticConnector(self,'V1ExcL23ExcL4Connection',cortex_exc_l23,cortex_exc_l4,self.parameters.l23_cortex_exc.L23ExcL4ExcConnection).connect()
        #ModularProbabilisticConnector(self,'V1ExcL23InhL4Connection',cortex_exc_l23,cortex_inh_l4,self.parameters.l23_cortex_exc.L23ExcL4InhConnection).connect()
        
        
        
