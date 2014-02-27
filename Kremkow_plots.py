"""
docstring goes here
"""
import numpy
import numpy.random
import mozaik.storage.queries as queries
from mozaik.visualization.plotting import (Plotting, GSynPlot,
                                           VmPlot, ConductanceSignalListPlot,
                                           AnalogSignalListPlot,OverviewPlot,PerNeuronValueScatterPlot,PlotTuningCurve,PerNeuronValuePlot)
from parameters import ParameterSet
import matplotlib.gridspec as gridspec
from mozaik.visualization.simple_plot import SpikeRasterPlot, SpikeHistogramPlot



class KremkowOverviewFigure(Plotting):
    required_parameters = ParameterSet({
        'sheet_name': str,  # the name of the sheet for which to plot
        'neuron': int,  # which neuron to show
    })

    def subplot(self, subplotspec):
        plots = {}
        gs = gridspec.GridSpecFromSubplotSpec(12, 18, subplot_spec=subplotspec,
                                              hspace=1.0, wspace=1.0)
        
        assert queries.equal_stimulus(self.datastore,['trial'])        
        
        lgn_on_dsv = queries.param_filter_query(self.datastore, sheet_name='X_ON')
        lgn_off_dsv = queries.param_filter_query(self.datastore, sheet_name='X_OFF')
        lgn_spikes_a = [[s.spiketrains[0:2] for s in lgn_on_dsv.get_segments()],
                      [s.spiketrains[0:2] for s in lgn_off_dsv.get_segments()]]
        lgn_spikes_b = [[s.spiketrains[2:4] for s in lgn_on_dsv.get_segments()],
                      [s.spiketrains[2:4] for s in lgn_off_dsv.get_segments()]]
        
        plots['LGN_SRP1'] = (SpikeRasterPlot(lgn_spikes_a),gs[1:4, 0:5],{'x_axis' : False, 'x_label': None,'colors':['#FACC2E', '#0080FF']})
        plots['LGN_SHP1'] = (SpikeHistogramPlot(lgn_spikes_a),gs[4:5, 0:5],{'x_axis' : False, 'x_label': None,'colors':['#FACC2E', '#0080FF']})
        plots['LGN_SRP2'] = (SpikeRasterPlot(lgn_spikes_b),gs[7:10, 0:5],{'x_axis' : False, 'x_label': None,'colors':['#FACC2E', '#0080FF']})
        plots['LGN_SHP2'] = (SpikeHistogramPlot(lgn_spikes_b),gs[10:11, 0:5],{'colors':['#FACC2E', '#0080FF']})
                     
        dsv1 = queries.param_filter_query(self.datastore,sheet_name=self.parameters.sheet_name)
        sp = [[[s.get_spiketrain(self.parameters.neuron)] for s in dsv1.get_segments()]]
        
        plots['V1_SRP1'] = (SpikeRasterPlot(sp),gs[:3, 6:14],{'x_axis' : False, 'x_label': None})
        plots['V1_SHP1'] = (SpikeRasterPlot(sp),gs[:3, 6:14],{'x_axis' : False, 'x_label': None})

        p = {}
        p['title']=None
        p['x_axis']=None
        p['x_label']=None
        plots['Vm_Plot'] = (VmPlot(dsv1, ParameterSet({'sheet_name': self.parameters.sheet_name,'neuron': self.parameters.neuron})),gs[4:8, 6:14],p)                                  
        p = {}
        p['title']=None
        plots['Gsyn_Plot'] = (GSynPlot(dsv1, ParameterSet({'sheet_name': self.parameters.sheet_name,'neuron': self.parameters.neuron})),gs[8:12, 6:14],p)                                  
        
        
        
        plots['GSTA_Plot'] = (ConductanceSignalListPlot(queries.TagBasedQuery(ParameterSet({'tags': ['GSTA']})).query(dsv1),ParameterSet({'normalize_individually': True})),gs[7:10, 15:],{})                                  
        
        #p = {}
        #p['mean'] = False
        #AnalogSignalListPlot(dsv, ParameterSet({'sheet_name': self.parameters.sheet_name,
        #                                        'ylabel': 'AC (norm)'})).subplot(gs[2:5, 15:], p)

        return plots

class StimulusResponseComparison(Plotting):
    required_parameters = ParameterSet({
        'sheet_name': str,  # the name of the sheet for which to plot
        'neuron': int,  # which neuron to show
    })

    def subplot(self, subplotspec):
        plots = {}
        gs = gridspec.GridSpecFromSubplotSpec(6, 18, subplot_spec=subplotspec,
                                              hspace=1.0, wspace=1.0)
        
        dsv = queries.param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=numpy.pi/2,contrast=100)            
        plots['Gratings'] = (OverviewPlot(dsv, ParameterSet({'sheet_name': self.parameters.sheet_name,'neuron': self.parameters.neuron, 'sheet_activity' : {}})),gs[0:2,:],{'x_label': None})
        dsv = queries.param_filter_query(self.datastore,st_name='NaturalImageWithEyeMovement')            
        plots['NIwEM'] = (OverviewPlot(dsv, ParameterSet({'sheet_name': self.parameters.sheet_name,'neuron': self.parameters.neuron, 'sheet_activity' : {}})),gs[2:4,:],{'x_label': None})
        dsv = queries.param_filter_query(self.datastore,st_name='DriftingGratingWithEyeMovement')            
        plots['GratingsWithEM'] = (OverviewPlot(dsv, ParameterSet({'sheet_name': self.parameters.sheet_name,'neuron': self.parameters.neuron, 'sheet_activity' : {}})),gs[4:6,:],{})
        
        
        return plots
        

class OrientationTuningSummary(Plotting):
    required_parameters = ParameterSet({})

    def subplot(self, subplotspec):
        plots = {}
        gs = gridspec.GridSpecFromSubplotSpec(20, 14, subplot_spec=subplotspec,
                                              hspace=1.0, wspace=1.0)
        
        
        analog_ids = sorted(numpy.random.permutation(queries.param_filter_query(self.datastore,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()))
        analog_ids_inh = sorted(numpy.random.permutation(queries.param_filter_query(self.datastore,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()))
        
        dsv1 = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['LGNAfferentOrientation'],sheet_name='V1_Exc_L4')
        dsv2 = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['orientation preference'],sheet_name='V1_Exc_L4',st_contrast=100,analysis_algorithm='GaussianTuningCurveFit')
        plots['Or_corr_exc'] = (PerNeuronValueScatterPlot(dsv1+dsv2, ParameterSet({})),gs[0:4,0:2],{'x_label' : None,'y_label' : None,'x_lim': (-0.1,numpy.pi),'y_lim' : (-0.1,numpy.pi)})
        dsv1 = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['LGNAfferentOrientation'],sheet_name='V1_Inh_L4')
        dsv2 = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['orientation preference'],sheet_name='V1_Inh_L4',st_contrast=100,analysis_algorithm='GaussianTuningCurveFit')
        plots['Or_corr_ing'] = (PerNeuronValueScatterPlot(dsv1+dsv2, ParameterSet({})),gs[0:4,2:4],{'x_label' : None,'y_label' : None,'x_lim': (-0.1,numpy.pi),'y_lim' : (-0.1,numpy.pi)})
        
        
                
        dsv = queries.param_filter_query(self.datastore,value_name=['orientation HWHH'],sheet_name=['V1_Exc_L4','V1_Inh_L4'])    
        plots['HWHH'] = (PerNeuronValueScatterPlot(dsv, ParameterSet({})),gs[0:4,5:9],{'x_lim': (-5,60),'y_lim' : (-5,60),'identity_line' : True})
        
        dsv = queries.param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate'])
        plots['ExcORTCMean'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids), 'sheet_name' : 'V1_Exc_L4','centered'  : True,'mean' : True})),gs[6:10,:3],{'title' : None,'x_label' : None})
        plots['ExcORTC1'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[0:7]), 'sheet_name' : 'V1_Exc_L4','centered'  : True,'mean' : False})),gs[6:8,3:],{'title' : None,'x_label' : None})
        plots['ExcORTC2'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[7:14]), 'sheet_name' : 'V1_Exc_L4','centered'  : True,'mean' : False})),gs[8:10,3:],{'title' : None,'x_label' : None})

        plots['InhORTCMean'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh), 'sheet_name' : 'V1_Inh_L4','centered'  : True,'mean' : True})),gs[11:15,:3],{'title' : None,'x_label' : None})
        plots['InhORTC1'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh[0:7]), 'sheet_name' : 'V1_Inh_L4','centered'  : True,'mean' : False})),gs[11:13,3:],{'title' : None,'x_label' : None})
        plots['InhORTC2'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh[7:14]), 'sheet_name' : 'V1_Inh_L4','centered'  : True,'mean' : False})),gs[13:15,3:],{'title' : None})

        
        dsv = queries.param_filter_query(self.datastore,value_name=['orientation HWHH'])    
        plots['HWHHHistogram'] = (PerNeuronValuePlot(dsv, ParameterSet({'cortical_view' : False})),gs[17:,:],{})


        return plots
        
