# -*- coding: utf-8 -*-
from parameters import ParameterSet
from mozaik.models import Model
from mozaik.connectors.fast import UniformProbabilisticArborization
from mozaik import load_component

class Boustani2007(Model):
    
    """
    This model is a re-implementation of a model by Jens Kremkow presented in following paper:
    
    S. El Boustani, M. Pospischil, M. Rudolph-Lilith, & A. Destexhe

    [Activated cortical states: experiments, analyses and models.](http://www.sciencedirect.com/science/article/pii/S0928425707000319)

    Journal of Physiology, 2007, Paris, 101(1–3), 99–109. 

    DOI: https://doi.org/10.1016/j.jphysparis.2007.10.001
    """

    required_parameters = ParameterSet({
	'sheets' : ParameterSet({
            'l4_cortex_exc' : ParameterSet, 
	    'l4_cortex_inh' : ParameterSet, 
		})
    })
    
    def __init__(self, sim, num_threads, parameters):
        Model.__init__(self, sim, num_threads, parameters)        
        # Load components
        CortexExcL4 = load_component(self.parameters.sheets.l4_cortex_exc.component)
        CortexInhL4 = load_component(self.parameters.sheets.l4_cortex_inh.component)
        
        cortex_exc_l4 = CortexExcL4(self, self.parameters.sheets.l4_cortex_exc.params)
        cortex_inh_l4 = CortexInhL4(self, self.parameters.sheets.l4_cortex_inh.params)

        # initialize projections
        UniformProbabilisticArborization(self,'V1L4ExcL4ExcConnection',cortex_exc_l4,cortex_exc_l4,self.parameters.sheets.l4_cortex_exc.L4ExcL4ExcConnection).connect()
        UniformProbabilisticArborization(self,'V1L4ExcL4InhConnection',cortex_exc_l4,cortex_inh_l4,self.parameters.sheets.l4_cortex_exc.L4ExcL4InhConnection).connect()
        UniformProbabilisticArborization(self,'V1L4InhL4ExcConnection',cortex_inh_l4,cortex_exc_l4,self.parameters.sheets.l4_cortex_inh.L4InhL4ExcConnection).connect()
        UniformProbabilisticArborization(self,'V1L4InhL4InhConnection',cortex_inh_l4,cortex_inh_l4,self.parameters.sheets.l4_cortex_inh.L4InhL4InhConnection).connect()

