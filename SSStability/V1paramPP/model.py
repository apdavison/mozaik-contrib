from parameters import ParameterSet
from mozaik.models import Model
from mozaik.connectors.fast import UniformProbabilisticArborization
from mozaik import load_component

class VogelsAbbott(Model):
    
    required_parameters = ParameterSet({
        'exc1' : ParameterSet, 
        'inh1' : ParameterSet, 
        'exc2' : ParameterSet, 
        'inh2' : ParameterSet, 
    })
    
    def __init__(self, sim, num_threads, parameters):
        Model.__init__(self, sim, num_threads, parameters)
        # Load components
        CortexExc1 = load_component(self.parameters.exc1.component)
	CortexExc2 = load_component(self.parameters.exc2.component)
	CortexInh1 = load_component(self.parameters.inh1.component)
	CortexInh2 = load_component(self.parameters.inh2.component)

        cortex_exc1 = CortexExc1(self, self.parameters.exc1.params)
        cortex_exc2 = CortexExc2(self, self.parameters.exc2.params)
        cortex_inh1 = CortexInh1(self, self.parameters.inh1.params)
        cortex_inh2 = CortexInh2(self, self.parameters.inh2.params)

        # initialize projections
        UniformProbabilisticArborization(self,'V1Exc1Exc1Connection',cortex_exc1,cortex_exc1,self.parameters.exc1.Exc1Exc1Connection).connect()
        UniformProbabilisticArborization(self,'V1Exc1Exc2Connection',cortex_exc1,cortex_exc2,self.parameters.exc1.Exc1Exc2Connection).connect()
        UniformProbabilisticArborization(self,'V1Exc1Inh1Connection',cortex_exc1,cortex_inh1,self.parameters.exc1.Exc1Inh1Connection).connect()
        UniformProbabilisticArborization(self,'V1Exc1Inh2Connection',cortex_exc1,cortex_inh2,self.parameters.exc1.Exc1Inh2Connection).connect()

        UniformProbabilisticArborization(self,'V1Exc2Exc1Connection',cortex_exc2,cortex_exc1,self.parameters.exc2.Exc2Exc1Connection).connect()
        UniformProbabilisticArborization(self,'V1Exc2Exc2Connection',cortex_exc2,cortex_exc2,self.parameters.exc2.Exc2Exc2Connection).connect()
        UniformProbabilisticArborization(self,'V1Exc2Inh1Connection',cortex_exc2,cortex_inh1,self.parameters.exc2.Exc2Inh1Connection).connect()
        UniformProbabilisticArborization(self,'V1Exc2Inh2Connection',cortex_exc2,cortex_inh2,self.parameters.exc2.Exc2Inh2Connection).connect()

        UniformProbabilisticArborization(self,'V1Inh1Exc1Connection',cortex_inh1,cortex_exc1,self.parameters.inh1.Inh1Exc1Connection).connect()
        UniformProbabilisticArborization(self,'V1Inh1Exc2Connection',cortex_inh1,cortex_exc2,self.parameters.inh1.Inh1Exc2Connection).connect()
        UniformProbabilisticArborization(self,'V1Inh1Inh1Connection',cortex_inh1,cortex_inh1,self.parameters.inh1.Inh1Inh1Connection).connect()
        UniformProbabilisticArborization(self,'V1Inh1Inh2Connection',cortex_inh1,cortex_inh2,self.parameters.inh1.Inh1Inh2Connection).connect()

        UniformProbabilisticArborization(self,'V1Inh2Exc1Connection',cortex_inh2,cortex_exc1,self.parameters.inh2.Inh2Exc1Connection).connect()
        UniformProbabilisticArborization(self,'V1Inh2Exc2Connection',cortex_inh2,cortex_exc2,self.parameters.inh2.Inh2Exc2Connection).connect()
        UniformProbabilisticArborization(self,'V1Inh2Inh1Connection',cortex_inh2,cortex_inh1,self.parameters.inh2.Inh2Inh1Connection).connect()
        UniformProbabilisticArborization(self,'V1Inh2Inh2Connection',cortex_inh2,cortex_inh2,self.parameters.inh2.Inh2Inh2Connection).connect()

