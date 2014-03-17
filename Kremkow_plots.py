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
        plots['GSTA_Plot'] = (ConductanceSignalListPlot(queries.TagBasedQuery(ParameterSet({'tags': ['GSTA']})).query(dsv1),ParameterSet({'normalize_individually': True, 'neurons' : [self.parameters.neuron]})),gs[7:10, 15:],{})                                  
        
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
                
        dsv = param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGrating',st_orientation=orr[numpy.argmin([circular_dist(o,l4_exc_or[0].get_value_by_id(self.parameters.neuron),numpy.pi)  for o in orr])],st_contrast=100)             
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
        gs = gridspec.GridSpecFromSubplotSpec(20, 15, subplot_spec=subplotspec,
                                              hspace=1.0, wspace=1.0)
        
        
        analog_ids = sorted(numpy.random.permutation(queries.param_filter_query(self.datastore,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()))
        analog_ids_inh = sorted(numpy.random.permutation(queries.param_filter_query(self.datastore,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()))
        
        dsv1 = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['LGNAfferentOrientation'],sheet_name='V1_Exc_L4')
        dsv2 = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['orientation preference'],sheet_name='V1_Exc_L4',st_contrast=100,analysis_algorithm='GaussianTuningCurveFit')
        plots['Or_corr_exc'] = (PerNeuronValueScatterPlot(dsv1+dsv2, ParameterSet({})),gs[0:4,3:5],{'x_label' : 'OR measured','y_label' : 'OR set','x_lim': (0.0,numpy.pi),'y_lim' : (0.0,numpy.pi)})
        dsv1 = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['LGNAfferentOrientation'],sheet_name='V1_Inh_L4')
        dsv2 = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['orientation preference'],sheet_name='V1_Inh_L4',st_contrast=100,analysis_algorithm='GaussianTuningCurveFit')
        plots['Or_corr_ing'] = (PerNeuronValueScatterPlot(dsv1+dsv2, ParameterSet({})),gs[0:4,5:7],{'x_label' : 'OR measured','y_label' : None,'x_lim': (0.0,numpy.pi),'y_lim' : (0.0,numpy.pi)})
        
        
                
        dsv = queries.param_filter_query(self.datastore,value_name=['orientation HWHH'],sheet_name=['V1_Exc_L4','V1_Inh_L4'])    
        plots['HWHH'] = (PerNeuronValueScatterPlot(dsv, ParameterSet({})),gs[0:4,8:12],{'x_lim': (0,50),'y_lim' : (0,50),'identity_line' : True, 'x_label' : 'HWHH Cont. 100%','y_label' : 'HWHH Cont. 50%'})
        
        dsv = queries.param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate'])
        plots['ExcORTCMean'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids), 'sheet_name' : 'V1_Exc_L4','centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[6:10,:3],{'title' : None,'x_label' : None , 'y_label' : 'EXC\nfiring rate (sp/s)'})
        plots['ExcORTC1'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[0:7]), 'sheet_name' : 'V1_Exc_L4','centered'  : True,'mean' : False,'pool' : False,'polar' : False})),gs[6:8,3:],{'title' : None,'x_label' : None,'y_axis' : False,'x_axis' : False, 'x_ticks' : False})
        plots['ExcORTC2'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[7:14]), 'sheet_name' : 'V1_Exc_L4','centered'  : True,'mean' : False,'pool' : False,'polar': False})),gs[8:10,3:],{'title' : None,'x_label' : None,'y_axis' : False,'x_axis' : False})

        plots['InhORTCMean'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh), 'sheet_name' : 'V1_Inh_L4','centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[11:15,:3],{'title' : None, 'y_label' : 'INH\nfiring rate (sp/s)'})
        plots['InhORTC1'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh[0:7]), 'sheet_name' : 'V1_Inh_L4','centered'  : True,'mean' : False,'pool' : False,'polar' : False})),gs[11:13,3:],{'title' : None,'x_label' : None,'y_axis' : False,'x_axis' : False})
        plots['InhORTC2'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh[7:14]), 'sheet_name' : 'V1_Inh_L4','centered'  : True,'mean' : False,'pool' : False,'polar' : False})),gs[13:15,3:],{'title' : None,'y_axis' : None})

        
        dsv = queries.param_filter_query(self.datastore,value_name=['orientation HWHH'],sheet_name=['V1_Exc_L4'])    
        plots['HWHHHistogramExc'] = (PerNeuronValuePlot(dsv, ParameterSet({'cortical_view' : False})),gs[17:,1:7],{'title' : 'Excitatory' , 'x_lim' : (0.0,50.0), 'x_label' : 'HWHH'})
        dsv = queries.param_filter_query(self.datastore,value_name=['orientation HWHH'],sheet_name=['V1_Inh_L4'])    
        plots['HWHHHistogramInh'] = (PerNeuronValuePlot(dsv, ParameterSet({'cortical_view' : False})),gs[17:,8:14],{'title' : 'Inhibitory' , 'x_lim' : (0.0,50.0), 'x_label' : 'HWHH'})

        return plots
        

class ConductanceAndVmTuningSummary(Plotting):
    required_parameters = ParameterSet({
        'many' : bool, # If true it will show 4 times as many (but twice as small) example neurons
    })

    def subplot(self, subplotspec):
        plots = {}
        gs = gridspec.GridSpecFromSubplotSpec(16, 14, subplot_spec=subplotspec,
                                              hspace=1.0, wspace=1.0)
        
        analog_ids = sorted(numpy.random.permutation(queries.param_filter_query(self.datastore,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()))
        analog_ids_inh = sorted(numpy.random.permutation(queries.param_filter_query(self.datastore,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()))
        
        dsv = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['F0_Exc_Cond','F0_Inh_Cond'],sheet_name='V1_Exc_L4')
        plots['MeanF0'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids), 'sheet_name' : 'V1_Exc_L4','centered'  : True,'mean' : True,'pool' : True,'polar' : True})),gs[:4,:3],{'legend' : True, 'y_label': 'F0(Cond)' ,'title' : None, 'x_ticks' : None, 'x_label' : None,'colors': {'F0_Exc_Cond contrast : 100' : '#FF0000' , 'F0_Exc_Cond contrast : 50' : '#FFACAC','F0_Inh_Cond contrast : 100' : '#0000FF' , 'F0_Inh_Cond contrast : 50' : '#ACACFF'}})

        dsv = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['F1_Exc_Cond','F1_Inh_Cond'],sheet_name='V1_Exc_L4')
        plots['MeanF1'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids), 'sheet_name' : 'V1_Exc_L4','centered'  : True,'mean' : True,'pool' : True,'polar' : True})),gs[4:8,:3],{'y_label': 'F1(Cond)','title' : None, 'x_ticks' : None, 'x_label' : None,'colors': {'F1_Exc_Cond contrast : 100' : '#FF0000' , 'F1_Exc_Cond contrast : 50' : '#FFACAC','F1_Inh_Cond contrast : 100' : '#0000FF' , 'F1_Inh_Cond contrast : 50' : '#ACACFF'}})

        dsv = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['F0_Vm'],sheet_name='V1_Exc_L4')
        plots['MeanVMF0'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids), 'sheet_name' : 'V1_Exc_L4','centered'  : True,'mean' : True,'pool' : True,'polar' : True})),gs[8:12,:3],{'y_label': 'F0(Vm)' ,'title' : None, 'x_ticks' : None, 'x_label' : None,'colors': {'F0_Vm contrast : 100' : '#000000' , 'F0_Vm contrast : 50' : '#ACACAC'}})

        dsv = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['F1_Vm'],sheet_name='V1_Exc_L4')
        plots['MeanVMF1'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids), 'sheet_name' : 'V1_Exc_L4','centered'  : True,'mean' : True,'pool' : True,'polar' : True})),gs[12:16,:3],{'y_label': 'F1(Vm)','title' : None, 'x_ticks' : None, 'x_label' : None,'colors': {'F1_Vm contrast : 100' : '#000000' , 'F1_Vm contrast : 50' : '#ACACAC'}})

        if self.parameters.many:

            dsv = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['F0_Exc_Cond','F0_Inh_Cond'],sheet_name='V1_Exc_L4')
            plots['F0a'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[20:30]), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : False,'pool' : True,'polar' : True})),gs[:2,3:],{'y_label': None,'title' : None, 'x_ticks' : None, 'x_label' : None,'colors': {'F0_Exc_Cond contrast : 100' : '#FF0000' , 'F0_Exc_Cond contrast : 50' : '#FFACAC','F0_Inh_Cond contrast : 100' : '#0000FF' , 'F0_Inh_Cond contrast : 50' : '#ACACFF'}})

            dsv = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['F0_Exc_Cond','F0_Inh_Cond'],sheet_name='V1_Exc_L4')
            plots['F0b'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[30:40]), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : False,'pool' : True,'polar' : True})),gs[2:4,3:],{'y_label': None,'title' : None, 'x_ticks' : None, 'x_label' : None,'colors': {'F0_Exc_Cond contrast : 100' : '#FF0000' , 'F0_Exc_Cond contrast : 50' : '#FFACAC','F0_Inh_Cond contrast : 100' : '#0000FF' , 'F0_Inh_Cond contrast : 50' : '#ACACFF'}})


            dsv = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['F1_Exc_Cond','F1_Inh_Cond'],sheet_name='V1_Exc_L4')
            plots['F1a'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[20:30]), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : False,'pool' : True,'polar' : True})),gs[4:6,3:],{'y_label': None,'title' : None, 'x_ticks' : None, 'x_label' : None,'colors': {'F1_Exc_Cond contrast : 100' : '#FF0000' , 'F1_Exc_Cond contrast : 50' : '#FFACAC','F1_Inh_Cond contrast : 100' : '#0000FF' , 'F1_Inh_Cond contrast : 50' : '#ACACFF'}})

            dsv = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['F1_Exc_Cond','F1_Inh_Cond'],sheet_name='V1_Exc_L4')
            plots['F1b'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[30:40]), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : False,'pool' : True,'polar' : True})),gs[6:8,3:],{'y_label': None,'title' : None, 'x_ticks' : None, 'x_label' : None,'colors': {'F1_Exc_Cond contrast : 100' : '#FF0000' , 'F1_Exc_Cond contrast : 50' : '#FFACAC','F1_Inh_Cond contrast : 100' : '#0000FF' , 'F1_Inh_Cond contrast : 50' : '#ACACFF'}})


            dsv = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['F0_Vm'],sheet_name='V1_Exc_L4')
            plots['VMF0a'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[20:30]), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : False,'pool' : True,'polar' : True})),gs[8:10,3:],{'y_label': None ,'title' : None, 'x_ticks' : None, 'x_label' : None,'colors': {'F0_Vm contrast : 100' : '#000000' , 'F0_Vm contrast : 50' : '#ACACAC'}})

            dsv = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['F0_Vm'],sheet_name='V1_Exc_L4')
            plots['VMF0b'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[30:40]), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : False,'pool' : True,'polar' : True})),gs[10:12,3:],{'y_label': None ,'title' : None, 'x_ticks' : None, 'x_label' : None,'colors': {'F0_Vm contrast : 100' : '#000000' , 'F0_Vm contrast : 50' : '#ACACAC'}})


            dsv = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['F1_Vm'],sheet_name='V1_Exc_L4')
            plots['VMF1a'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[20:30]), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : False,'pool' : True,'polar' : True})),gs[12:14,3:],{'y_label': None,'title' : None, 'x_ticks' : None, 'x_label' : None,'colors': {'F1_Vm contrast : 100' : '#000000' , 'F1_Vm contrast : 50' : '#ACACAC'}})

            dsv = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['F1_Vm'],sheet_name='V1_Exc_L4')
            plots['VMF1b'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[30:40]), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : False,'pool' : True,'polar' : True})),gs[14:16,3:],{'y_label': None,'title' : None, 'x_ticks' : None, 'x_label' : None,'colors': {'F1_Vm contrast : 100' : '#000000' , 'F1_Vm contrast : 50' : '#ACACAC'}})

        else:
            dsv = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['F0_Exc_Cond'],sheet_name='V1_Exc_L4')
            plots['F0E'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[:5]), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : False,'pool' : False,'polar' : False})),gs[:4,3:],{'title' : None, 'x_ticks' : None, 'x_label' : None,'colors': {'contrast : 100' : '#FF0000' , 'contrast : 50' : '#FFACAC'}})


            dsv = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['F1_Exc_Cond'],sheet_name='V1_Exc_L4')
            plots['F1E'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[:5]), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : False,'pool' : False,'polar' : False})),gs[5:9,3:],{'title' : None, 'x_ticks' : None, 'x_label' : None,'colors': {'contrast : 100' : '#FF0000' , 'contrast : 50' : '#FFACAC'}})

            dsv = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['F0_Inh_Cond'],sheet_name='V1_Exc_L4')
            plots['F0I'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[:5]), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : False,'pool' : False,'polar' : False})),gs[10:14,3:],{'title' : None, 'x_ticks' : None, 'x_label' : None,'colors': {'contrast : 100' : '#0000FF' , 'contrast : 50' : '#ACACFF'}})

            dsv = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name=['F1_Inh_Cond'],sheet_name='V1_Exc_L4')
            plots['F1I'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[:5]), 'sheet_name' : 'V1_Exc_L4','centered'  : False,'mean' : False,'pool' : False,'polar' : False})),gs[15:19,3:],{'title' : None, 'x_ticks' : None, 'x_label' : None,'colors': {'contrast : 100' : '#0000FF' , 'contrast : 50' : '#ACACFF'}})


        return plots
        

