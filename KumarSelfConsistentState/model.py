from parameters import ParameterSet
from mozaik.models import Model
from mozaik.connectors.fast import FixedKConnector
from mozaik import load_component

class KumarEtAl2007(Model):
    
    required_parameters = ParameterSet({
        'l4_cortex_exc' : ParameterSet, 
    })
    
    def __init__(self, sim, num_threads, parameters):
        Model.__init__(self, sim, num_threads, parameters)
        # Load components
        CortexExcL4 = load_component(self.parameters.l4_cortex_exc.component)
        cortex_exc_l4 = CortexExcL4(self, self.parameters.l4_cortex_exc.params)


