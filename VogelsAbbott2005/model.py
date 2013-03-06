import sys
sys.path.append('/home/jan/projects/mozaik/')
from NeuroTools.parameters import ParameterSet
from mozaik.models.model import Model
from mozaik.connectors.fast_connectors import UniformProbabilisticArborization
from mozaik.framework import load_component

class VogelsAbbott(Model):
    
    required_parameters = ParameterSet({
        'l4_cortex_exc' : ParameterSet, 
        'l4_cortex_inh' : ParameterSet, 
    })
    
    def __init__(self,simulator,parameters):
        Model.__init__(self,simulator,parameters)        
        # Load components
        CortexExcL4 = load_component(self.parameters.l4_cortex_exc.component)
        CortexInhL4 = load_component(self.parameters.l4_cortex_inh.component)
        
        cortex_exc_l4 = CortexExcL4(self, self.parameters.l4_cortex_exc.params)
        cortex_inh_l4 = CortexInhL4(self, self.parameters.l4_cortex_inh.params)

        # initialize projections
        UniformProbabilisticArborization(self,'V1L4ExcL4ExcConnection',cortex_exc_l4,cortex_exc_l4,self.parameters.l4_cortex_exc.L4ExcL4ExcConnection).connect()
        UniformProbabilisticArborization(self,'V1L4ExcL4InhConnection',cortex_exc_l4,cortex_inh_l4,self.parameters.l4_cortex_exc.L4ExcL4InhConnection).connect()
        UniformProbabilisticArborization(self,'V1L4InhL4ExcConnection',cortex_inh_l4,cortex_exc_l4,self.parameters.l4_cortex_inh.L4InhL4ExcConnection).connect()
        UniformProbabilisticArborization(self,'V1L4InhL4InhConnection',cortex_inh_l4,cortex_inh_l4,self.parameters.l4_cortex_inh.L4InhL4InhConnection).connect()
