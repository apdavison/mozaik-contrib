from parameters import ParameterSet
from mozaik.models import Model
from mozaik.connectors.fast import FixedKConnector
from mozaik import load_component

class VogelsAbbott(Model):
    
    required_parameters = ParameterSet({
        'l4_cortex_exc' : ParameterSet, 
        'l4_cortex_inh' : ParameterSet, 
    })
    
    def __init__(self, sim, num_threads, parameters):
        Model.__init__(self, sim, num_threads, parameters)
        # Load components
        CortexExcL4 = load_component(self.parameters.l4_cortex_exc.component)
        CortexInhL4 = load_component(self.parameters.l4_cortex_inh.component)
        
        cortex_exc_l4 = CortexExcL4(self, self.parameters.l4_cortex_exc.params)
        cortex_inh_l4 = CortexInhL4(self, self.parameters.l4_cortex_inh.params)

        # initialize projections
        FixedKConnector(self,'V1L4ExcL4ExcConnection',cortex_exc_l4,cortex_exc_l4,self.parameters.l4_cortex_exc.L4ExcL4ExcConnection).connect()
        FixedKConnector(self,'V1L4ExcL4InhConnection',cortex_exc_l4,cortex_inh_l4,self.parameters.l4_cortex_exc.L4ExcL4InhConnection).connect()
        FixedKConnector(self,'V1L4InhL4ExcConnection',cortex_inh_l4,cortex_exc_l4,self.parameters.l4_cortex_inh.L4InhL4ExcConnection).connect()
        FixedKConnector(self,'V1L4InhL4InhConnection',cortex_inh_l4,cortex_inh_l4,self.parameters.l4_cortex_inh.L4InhL4InhConnection).connect()

