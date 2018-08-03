import pylab
import numpy
import pickle
import math
import mozaik.storage.queries as queries
import matplotlib.gridspec as gridspec
from mozaik.visualization.plotting import Plotting
from mozaik.visualization.helper_functions import *
from parameters import ParameterSet
from mozaik.storage.queries import *
from mozaik.analysis.analysis import *
from mozaik.controller import Global
from mozaik.visualization.plotting import (Plotting, GSynPlot,RasterPlot,PerNeuronAnalogSignalScatterPlot,
                                           VmPlot, ConductanceSignalListPlot,ScatterPlot,
                                           AnalogSignalListPlot,OverviewPlot,PerNeuronValueScatterPlot,PlotTuningCurve,PerNeuronValuePlot,CorticalColumnRasterPlot)
from mozaik.visualization.simple_plot import *
                                           
from mozaik.visualization.plot_constructors import LinePlot, PerStimulusPlot, PerStimulusADSPlot, ADSGridPlot                                           
from mozaik.tools.circ_stat import circular_dist
import mozaik.visualization.helper_functions as phf

from mozaik.controller import Global
logger = mozaik.getMozaikLogger()

class ResponseOverview(Plotting):

    required_parameters = ParameterSet({
        'sheet_name': str,  # the name of the sheet for which to plot
        'neuron': int,
    })
    
    def subplot(self, subplotspec):
        dsv = queries.param_filter_query(self.datastore,sheet_name=self.parameters.sheet_name)
        return PerStimulusPlot(dsv, function=self._ploter, title_style="Clever").make_line_plot(subplotspec)

    def _ploter(self, dsv,subplotspec):
        d = []
        gs = gridspec.GridSpecFromSubplotSpec(1, 7, subplot_spec=subplotspec,wspace=0.7)
        params = {}
        params['x_label']  = None
        params['x_tick_labels'] = None
        params['x_tick_style'] ='Custom'
        params['y_label'] = 'trial'

        scale = 0.08
        #lams = self.datastore.full_datastore.get_analysis_result(identifier='PerNeuronValue',analysis_algorithm ='NeuronAnnotationsToPerNeuronValues', value_name = 'Light activation magnitude(V1_Exc_L2/3,'+str(scale) + ',0.0,0.0,1.0,10)', sheet_name = 'V1_Exc_L2/3')
        for z in self.datastore.full_datastore.get_analysis_result(identifier='PerNeuronValue',analysis_algorithm ='NeuronAnnotationsToPerNeuronValues', sheet_name = 'V1_Exc_L2/3'):
            print z.value_name            

        lams = self.datastore.full_datastore.get_analysis_result(identifier='PerNeuronValue',analysis_algorithm ='NeuronAnnotationsToPerNeuronValues', value_name = 'Light activation magnitude(V1_Exc_L2/3,'+str(scale) + ',0.0,0.5,10)', sheet_name = 'V1_Exc_L2/3')
        fontsize = 15
        
        if len(lams) == 1:
            lams = lams[0]
            spike_ids = param_filter_query(self.datastore,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_spike_train_ids()
            dsv1 = param_filter_query(self.datastore,value_name=['Firing rate'],st_name='InternalStimulus',sheet_name='V1_Exc_L2/3',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_orientation=0,st_stimulating_signal_parameters_scale=scale)

            axis = pylab.subplot(gs[0,6])
            pylab.plot(lams.get_value_by_id(list(spike_ids)),numpy.array(dsv1.get_analysis_result()[0].get_value_by_id(list(spike_ids)))*2,'o',markersize=2)
            pylab.xlabel(r'photon flux $(ph/cm^2/s)$',fontsize=fontsize)
            #pylab.yticks([0,50,100],fontsize=fontsize)
            pylab.xticks([0,5e15],[0.0,r'$5 \times 10^{15}$'],fontsize=fontsize)
            pylab.xlim(0,5e15)
            pylab.ylabel('spikes/s',fontsize=fontsize)
            phf.disable_top_right_axis(axis)
            for label in axis.get_xticklabels() + axis.get_yticklabels():
                label.set_fontsize(fontsize)
        else:
            axis = pylab.subplot(gs[0,6])
            pylab.plot([0,1],[0,1],'o',markersize=2)
            pylab.xlabel(r'photon flux $(ph/cm^2/s)$',fontsize=fontsize)
            pylab.yticks([0,50,100],fontsize=fontsize)
            pylab.xticks([0,5e15],[0.0,r'$5 \times 10^{15}$'],fontsize=fontsize)
            pylab.xlim(0,5e15)
            pylab.ylabel('(spikes/s)',fontsize=fontsize)
            phf.disable_top_right_axis(axis)
            for label in axis.get_xticklabels() + axis.get_yticklabels():
                label.set_fontsize(fontsize)

        d.extend([ ("Spike_plot",RasterPlot(dsv,
                   ParameterSet({'sheet_name': self.parameters.sheet_name,
                                 'trial_averaged_histogram': False,'spontaneous' : False,
                                 'neurons': [self.parameters.neuron]})
                   ),gs[0, 0:2],{ 'title' : None}),
                
                 ("Conductance_plot",GSynPlot(dsv,
                   ParameterSet({'sheet_name': self.parameters.sheet_name,'spontaneous' : False,
                               'neuron': self.parameters.neuron})
                 ),gs[0, 2:4], { 'title' : None}),

                 ("Vm_plot",VmPlot(dsv,
                   ParameterSet({'sheet_name': self.parameters.sheet_name,'spontaneous' : False,
                             'neuron': self.parameters.neuron})
                 ),gs[0, 4:6], {'title' : None})
              ])
        return d




class LightStimulationOverview(Plotting):
    required_parameters = ParameterSet({})

    def subplot(self, subplotspec):
        plots = {}

        scale = 0.08
        scales = [0.005,0.01,0.04,0.08,0.16,0.32][::-1]

        gs = gridspec.GridSpecFromSubplotSpec(23,60, subplot_spec=subplotspec,hspace=1.0, wspace=100.0)

        analog_ids = param_filter_query(self.datastore,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_esyn_ids()
        spike_ids = param_filter_query(self.datastore,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_spike_train_ids()

        l23_exc_or = self.datastore.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L2/3')[0]
        x = self.datastore.get_neuron_postions()['V1_Exc_L2/3'][0]
        y = self.datastore.get_neuron_postions()['V1_Exc_L2/3'][1]
        depth = self.datastore.get_neuron_postions()['V1_Exc_L2/3'][2]
        ors = numpy.array(l23_exc_or.ids)[numpy.nonzero(numpy.array([circular_dist(o,0,numpy.pi) for o in l23_exc_or.values]) < 0.15)[0]]
        close_spikes = numpy.array([a for a in ors if (a in spike_ids)])
        close_spikes_sheet_indexes = self.datastore.get_sheet_indexes(sheet_name='V1_Exc_L2/3',neuron_ids=close_spikes)
        close_spikes = close_spikes[numpy.sqrt(x[close_spikes_sheet_indexes]**2 + y[close_spikes_sheet_indexes]**2)<1.5]

        depth_of_close_spikes = depth[self.datastore.get_sheet_indexes(sheet_name='V1_Exc_L2/3',neuron_ids=close_spikes)]
        depth_of_close_spikes ,close_spikes = zip(*sorted(zip(depth_of_close_spikes,close_spikes),reverse=True))
        close_analogs = list(set(close_spikes) & set(analog_ids))

        depth_of_all_spikes = depth[self.datastore.get_sheet_indexes(sheet_name='V1_Exc_L2/3',neuron_ids=spike_ids)]


        lams = self.datastore.get_analysis_result(identifier='PerNeuronValue',analysis_algorithm ='NeuronAnnotationsToPerNeuronValues', value_name = 'Light activation magnitude(V1_Exc_L2/3,'+str(scale) + ',0.0,0.5,10)', sheet_name = 'V1_Exc_L2/3')[0]
        lamsChR = self.datastore.get_analysis_result(identifier='PerNeuronValue',analysis_algorithm ='NeuronAnnotationsToPerNeuronValues', value_name = 'Light activation magnitude ChR('+str(scale) + ',0.0_0.5_10)', sheet_name = 'V1_Exc_L2/3')[0]

        dsv = param_filter_query(self.datastore,st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_orientation=0,st_stimulating_signal_parameters_scale=scale,st_trial=0)
        plots['ExampleRaster'] = (RasterPlot(dsv, ParameterSet({'sheet_name' : 'V1_Exc_L2/3', 'neurons' : list(close_spikes[:100]), 'trial_averaged_histogram' : True, 'spontaneous' : False})),gs[:,19:34],{})
        
       
        fontsize = 15
        axis = pylab.subplot(gs[:10,0:16])
        for sc in scales:
            segs = param_filter_query(self.datastore,st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',sheet_name='V1_Exc_L2/3',st_stimulating_signal_parameters_orientation=0,st_stimulating_signal_parameters_scale=sc,st_trial=0).get_segments()
            assert len(segs) == 1
            pylab.plot(segs[0].get_vm(close_analogs[0]))

        phf.disable_top_right_axis(axis)
        pylab.xlim(0,1000)
        pylab.ylim(-80,-40)
        pylab.yticks([-80,-60,-40],fontsize=fontsize)
        pylab.xticks([0,200,800,1000],[0,0.2,0.8,1.0],fontsize=fontsize)
        pylab.ylabel('Vm (mV)',fontsize=fontsize)
        for label in axis.get_xticklabels() + axis.get_yticklabels():
            label.set_fontsize(fontsize)
        
        
        axis = pylab.subplot(gs[13:,0:16])
        ids = self.datastore.get_sheet_indexes(sheet_name='V1_Exc_L2/3',neuron_ids=close_analogs)
        for sc in scales:
            f = open(Global.root_directory+'/mixed_signalsV1_Exc_L2_3_'+str(sc)+'_0.0_0.5_10.pickle')
            mixed_signals = pickle.load(f)
            pylab.plot(mixed_signals[ids[0]],label=('%.3f' % sc) + '$e^{16} (ph/cm^2/s)$'  )
        pylab.legend(loc='upper right',ncol=2)
        phf.disable_top_right_axis(axis)
        pylab.xlim(0,1000)
        pylab.xticks([0,200,800,1000],[0,0.2,0.8,1.0])
        pylab.xlabel('time (s)',fontsize=fontsize)
        pylab.yticks([0,0.15,0.3],fontsize=fontsize)
        pylab.ylim(-0.02,0.4)
        pylab.ylabel('inward current (nA)',fontsize=fontsize)
        for label in axis.get_xticklabels() + axis.get_yticklabels():
            label.set_fontsize(fontsize)
        
        dsv = param_filter_query(self.datastore,value_name=['Firing rate'],st_name='InternalStimulus',sheet_name='V1_Exc_L2/3',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_orientation=0,st_stimulating_signal_parameters_scale=scale)

        axis = pylab.subplot(gs[:10,37:47])
        pylab.plot(lams.get_value_by_id(list(spike_ids)),numpy.array(dsv.get_analysis_result()[0].get_value_by_id(list(spike_ids)))*2,'o')
        pylab.xlabel(r'photon flux $(ph/cm^2/s)$',fontsize=fontsize)
        pylab.yticks([0,50,100],fontsize=fontsize)
        pylab.xticks([0,1e16],[0.0,r'$10^{16}$'],fontsize=fontsize)
        pylab.xlim(0,1.0e16)
        pylab.ylabel('firing rate (spk/s)',fontsize=fontsize)
        phf.disable_top_right_axis(axis)
        for label in axis.get_xticklabels() + axis.get_yticklabels():
            label.set_fontsize(fontsize)


        axis = pylab.subplot(gs[13:,37:47])
        pylab.plot(lams.get_value_by_id(list(spike_ids)),lamsChR.get_value_by_id(list(spike_ids)),'o')
        pylab.xlabel(r'photon flux $(ph/cm^2/s)$',fontsize=fontsize)
        pylab.xticks([0,1e16],[0.0,r'$10^{16}$'],fontsize=fontsize)
        pylab.xlim(0,1.0e16)
        pylab.yticks([0,0.05,0.1],fontsize=fontsize)
        pylab.ylabel('inward current (nA)',fontsize=fontsize)
        phf.disable_top_right_axis(axis)
        for label in axis.get_xticklabels() + axis.get_yticklabels():
            label.set_fontsize(fontsize)


        axis = pylab.subplot(gs[:10,50:60])
        pylab.plot(depth_of_close_spikes,numpy.array(dsv.get_analysis_result()[0].get_value_by_id(list(close_spikes)))*2,'o')
        pylab.xticks([100,250,400],[100,250,400])
        pylab.xlabel(r'depth ($\mu m$)',fontsize=fontsize)
        pylab.yticks([0,50,100],fontsize=fontsize)
        pylab.ylabel('firing rate (spk/s)',fontsize=fontsize)
        pylab.xlim(100,400)
        phf.disable_top_right_axis(axis)
        for label in axis.get_xticklabels() + axis.get_yticklabels():
            label.set_fontsize(fontsize)



        axis = pylab.subplot(gs[13:,50:60])
        pylab.scatter(l23_exc_or.values-numpy.pi/2,lams.values,s=3)
        pylab.xticks([-numpy.pi/2, 0, numpy.pi/2],["-$\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$"])
        pylab.xlabel('orientation',fontsize=fontsize)
        pylab.ylabel(r'photon flux $(ph/cm^2/s)$',fontsize=fontsize)
        pylab.yticks([0,1e16],[0.0,r'$10^{16}$'],fontsize=fontsize)
        pylab.ylim(0,1.0e16)
        pylab.xlim([-numpy.pi/2, numpy.pi/2])
        phf.disable_top_right_axis(axis)
        for label in axis.get_xticklabels() + axis.get_yticklabels():
            label.set_fontsize(fontsize)

        return plots






class OrientationTuningSummaryFiringRates(Plotting):
    required_parameters = ParameterSet({
              'cortical_stimulation'   : bool,
              'scale'                  : bool
        })

    def subplot(self, subplotspec):
        plots = {}
        gs = gridspec.GridSpecFromSubplotSpec(6, 29, subplot_spec=subplotspec,
                                              hspace=1.0, wspace=5.0)
        
        spike_ids = sorted(numpy.random.permutation(queries.param_filter_query(self.datastore,sheet_name='V1_Exc_L2/3').get_segments()[0].get_stored_spike_train_ids()))

        scales = [0.01,0.07,0.14]
        contrasts = [3,5,10,100]
        fontsize = 15
        
        if self.parameters.cortical_stimulation:
            if self.parameters.scale:
                base = queries.param_filter_query(self.datastore,sheet_name='V1_Exc_L2/3',st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_scale = scales[1],value_name=['stimulating_signal_parameters_orientation baseline of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(spike_ids)
                mmax = queries.param_filter_query(self.datastore,sheet_name='V1_Exc_L2/3',st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_scale = scales[1],value_name=['stimulating_signal_parameters_orientation max of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(spike_ids)
            else:
                base = queries.param_filter_query(self.datastore,sheet_name='V1_Exc_L2/3',st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_contrast = contrasts[-1],value_name=['stimulating_signal_parameters_orientation baseline of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(spike_ids)
                mmax = queries.param_filter_query(self.datastore,sheet_name='V1_Exc_L2/3',st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_contrast = contrasts[-1],value_name=['stimulating_signal_parameters_orientation max of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(spike_ids)                    
        else:
            base = queries.param_filter_query(self.datastore,sheet_name='V1_Exc_L2/3',st_name=['FullfieldDriftingSinusoidalGratingA'],st_contrast=contrasts[-1],value_name=['orientation baseline of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(spike_ids)
            mmax = queries.param_filter_query(self.datastore,sheet_name='V1_Exc_L2/3',st_name=['FullfieldDriftingSinusoidalGratingA'],st_contrast=contrasts[-1],value_name=['orientation max of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(spike_ids)
        #responsive_spike_ids = numpy.array(spike_ids)[numpy.array(base)+numpy.array(mmax) > 0.0]
        responsive_spike_ids = spike_ids
                
        if self.parameters.cortical_stimulation:
            if self.parameters.scale:
                dsv = queries.param_filter_query(self.datastore,st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',value_name='Firing rate',analysis_algorithm=['TrialAveragedFiringRate'],st_stimulating_signal_parameters_scale = scales)
            else:
                dsv = queries.param_filter_query(self.datastore,st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',value_name='Firing rate',analysis_algorithm=['TrialAveragedFiringRate'],st_stimulating_signal_parameters_contrast = contrasts)

            plots['ORTCMean'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'stimulating_signal_parameters_orientation', 'neurons': list(responsive_spike_ids), 'sheet_name' : "V1_Exc_L2/3",'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[:,:6],{'x_label' : 'orientation','y_lim' : (0,None),'title' : None, 'y_label' : 'firing rate (sp/s)'})
            plots['ORTC1'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'stimulating_signal_parameters_orientation', 'neurons': list(responsive_spike_ids[0:3]), 'sheet_name' : "V1_Exc_L2/3",'centered'  : True,'mean' : False,'pool' : False,'polar' : False})),gs[0:3,6:15],{'y_lim' : (0,None),'title' : None,'left_border' : None, 'x_label' : None,'y_axis' : False,'x_axis' : False})
            plots['ORTC2'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'stimulating_signal_parameters_orientation', 'neurons': list(responsive_spike_ids[3:6]), 'sheet_name' : "V1_Exc_L2/3",'centered'  : True,'mean' : False,'pool' : False,'polar' : False})),gs[3:,6:15],{'y_lim' : (0,None),'title' : None,'left_border' : None, 'y_axis' : None,'x_axis' : False})
        else:
            dsv = queries.param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGratingA',value_name='Firing rate',analysis_algorithm=['TrialAveragedFiringRate'])
            plots['ORTCMean'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(responsive_spike_ids), 'sheet_name' : "V1_Exc_L2/3",'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[:,:6],{'y_lim' : (0,None),'title' : None, 'y_label' : 'firing rate (sp/s)'})
            plots['ORTC1'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(responsive_spike_ids[0:3]), 'sheet_name' : "V1_Exc_L2/3",'centered'  : True,'mean' : False,'pool' : False,'polar' : False})),gs[0:3,6:15],{'y_lim' : (0,None),'title' : None,'left_border' : None, 'x_label' : None,'y_axis' : False,'x_axis' : False})
            plots['ORTC2'] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'orientation', 'neurons': list(responsive_spike_ids[3:6]), 'sheet_name' : "V1_Exc_L2/3",'centered'  : True,'mean' : False,'pool' : False,'polar' : False})),gs[3:,6:15],{'y_lim' : (0,None),'title' : None,'left_border' : None, 'y_axis' : None,'x_axis' : False})

        if self.parameters.cortical_stimulation:   
            if self.parameters.scale:   
                low_c = queries.param_filter_query(self.datastore,value_name=['stimulating_signal_parameters_orientation HWHH of Firing rate'],sheet_name=['V1_Exc_L2/3'],st_stimulating_signal_parameters_scale=scales[1]).get_analysis_result()[0].get_value_by_id(responsive_spike_ids)
                high_c = queries.param_filter_query(self.datastore,value_name=['stimulating_signal_parameters_orientation HWHH of Firing rate'],sheet_name=['V1_Exc_L2/3'],st_stimulating_signal_parameters_scale=scales[2]).get_analysis_result()[0].get_value_by_id(responsive_spike_ids)
            else:
                low_c = queries.param_filter_query(self.datastore,value_name=['stimulating_signal_parameters_orientation HWHH of Firing rate'],sheet_name=['V1_Exc_L2/3'],st_stimulating_signal_parameters_contrast=contrasts[0]).get_analysis_result()[0].get_value_by_id(responsive_spike_ids)
                high_c = queries.param_filter_query(self.datastore,value_name=['stimulating_signal_parameters_orientation HWHH of Firing rate'],sheet_name=['V1_Exc_L2/3'],st_stimulating_signal_parameters_contrast=contrasts[3]).get_analysis_result()[0].get_value_by_id(responsive_spike_ids)
        else:
                low_c = queries.param_filter_query(self.datastore,value_name=['orientation HWHH of Firing rate'],sheet_name=['V1_Exc_L2/3'],st_contrast=contrasts[0]).get_analysis_result()[0].get_value_by_id(responsive_spike_ids)
                high_c = queries.param_filter_query(self.datastore,value_name=['orientation HWHH of Firing rate'],sheet_name=['V1_Exc_L2/3'],st_contrast=contrasts[3]).get_analysis_result()[0].get_value_by_id(responsive_spike_ids)

        axis = pylab.subplot(gs[:,16:22])
        pylab.scatter(low_c,high_c,color='k')
        pylab.plot([0,40],[0,40],color='#AAAAAA')
        pylab.xlim(0,40)
        pylab.xticks([0,20,40])
        pylab.xlabel('HWHH low contrast',fontsize=fontsize)
        pylab.ylim(0,40)
        pylab.yticks([0,20,40])
        pylab.ylabel('HWHH  high contrast',fontsize=fontsize)
        for label in axis.get_xticklabels() + axis.get_yticklabels():
            label.set_fontsize(fontsize)
        phf.disable_top_right_axis(pylab.gca())

        if self.parameters.cortical_stimulation:
            if self.parameters.scale:
                dsv = queries.param_filter_query(self.datastore,value_name=['stimulating_signal_parameters_orientation HWHH of Firing rate'],sheet_name=['V1_Exc_L2/3'],st_stimulating_signal_parameters_scale=scales[-1])    
            else:
                dsv = queries.param_filter_query(self.datastore,value_name=['stimulating_signal_parameters_orientation HWHH of Firing rate'],sheet_name=['V1_Exc_L2/3'],st_stimulating_signal_parameters_contrast=contrasts[-1])    
        else:
            dsv = queries.param_filter_query(self.datastore,value_name=['orientation HWHH of Firing rate'],sheet_name=['V1_Exc_L2/3'],st_contrast=[100])    
        
        y = dsv.get_analysis_result()[0].get_value_by_id(responsive_spike_ids)
        axis = pylab.subplot(gs[:,23:])
        pylab.hist(y,color='k',bins=numpy.arange(0,50,3))
        pylab.xlim(0,50)
        pylab.xticks([0,25,50])
        pylab.ylim(0,300)
        pylab.yticks([0,150,300])
        pylab.ylabel('# neurons',fontsize=fontsize)
        for label in axis.get_xticklabels() + axis.get_yticklabels():
            label.set_fontsize(fontsize)
        phf.disable_top_right_axis(pylab.gca())
        pylab.xlabel('HWHH high contrast',fontsize=fontsize)

        return plots











class StatisticsOverview(Plotting):
    required_parameters = ParameterSet({
        'cortical_stimulation'   : bool,
        'type' : int,
        'window_length' : float,
    })
    
    def calculate_cc(self,sheet_name,neurons):
                if self.parameters.type == 0:
                    orr = list(set([MozaikParametrized.idd(s).orientation for s in queries.param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGratingA',st_contrast=100).get_stimuli()]))        
                    dsv1 =  queries.param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGratingA',st_contrast=100,sheet_name=sheet_name,analysis_algorithm='TrialToTrialCrossCorrelationOfAnalogSignalList')
                    oor = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name = 'orientation preference', sheet_name = sheet_name,st_contrast=100).get_analysis_result()

                else:
                    orr = list(set([MozaikParametrized.idd(s).direct_stimulation_parameters.stimulating_signal_parameters.orientation.value for s in queries.param_filter_query(self.datastore,st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_contrast=100).get_stimuli()]))        
                    dsv1 =  queries.param_filter_query(self.datastore,st_stimulating_signal_parameters_contrast=100,sheet_name=sheet_name,analysis_algorithm='TrialToTrialCrossCorrelationOfAnalogSignalList')
                    oor = queries.param_filter_query(self.datastore,identifier='PerNeuronValue',value_name = 'stimulating_signal_parameters_orientation preference', sheet_name = sheet_name,st_stimulating_signal_parameters_contrast=100).get_analysis_result()
                
                print oor
                vm_gr_asls = []
                psth_gr_asls = []
                
                if True:
                    for neuron_idd in neurons:
                        col = orr[numpy.argmin([circular_dist(o,oor[0].get_value_by_id(neuron_idd),numpy.pi)  for o in orr])]
                        if self.parameters.type == 0:
                            dsv =  queries.param_filter_query(dsv1,y_axis_name='trial-trial cross-correlation of Vm (no AP)',st_orientation=col,ads_unique=True)
                        else:
                            dsv =  queries.param_filter_query(dsv1,y_axis_name='trial-trial cross-correlation of Vm (no AP)',st_stimulating_signal_parameters_orientation=col,ads_unique=True)

                        vm_gr_asls.append(dsv.get_analysis_result()[0].get_asl_by_id(neuron_idd))
                        if self.parameters.type == 0:
                            dsv =  queries.param_filter_query(dsv1,y_axis_name='trial-trial cross-correlation of psth (bin=10.0)',st_orientation=col,ads_unique=True)
                        else:
                            dsv =  queries.param_filter_query(dsv1,y_axis_name='trial-trial cross-correlation of psth (bin=10.0)',st_stimulating_signal_parameters_orientation=col,ads_unique=True)

                        psth_gr_asls.append(dsv.get_analysis_result()[0].get_asl_by_id(neuron_idd))

                vm_cc_gr = numpy.mean(numpy.array(vm_gr_asls),axis=0)
                psth_cc_gr = numpy.mean(numpy.array(psth_gr_asls),axis=0)
                return vm_cc_gr,psth_cc_gr

    def subplot(self, subplotspec):
        plots = {}
        gs = gridspec.GridSpecFromSubplotSpec(2,3, subplot_spec=subplotspec,hspace=0.25, wspace=0.3)

        fontsize=15
        
        analog_ids = param_filter_query(self.datastore,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_esyn_ids()
        spike_ids = param_filter_query(self.datastore,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_spike_train_ids()

        l23_exc_or = self.datastore.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L2/3')[0]
        x = self.datastore.get_neuron_postions()['V1_Exc_L2/3'][0]
        y = self.datastore.get_neuron_postions()['V1_Exc_L2/3'][1]
        ors = numpy.array(l23_exc_or.ids)[numpy.nonzero(numpy.array([circular_dist(o,0,numpy.pi) for o in l23_exc_or.values]) < 0.15)[0]]
        close_spikes = numpy.array([a for a in ors if (a in spike_ids)])
        close_spikes_sheet_indexes = self.datastore.get_sheet_indexes(sheet_name='V1_Exc_L2/3',neuron_ids=close_spikes)
        close_spikes = close_spikes[numpy.sqrt(x[close_spikes_sheet_indexes]**2 + y[close_spikes_sheet_indexes]**2)<1.5]
        close_analogs = list(set(close_spikes) & set(analog_ids))

        if self.parameters.cortical_stimulation:
           dsv_pref = queries.param_filter_query(self.datastore,st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_orientation=0,st_stimulating_signal_parameters_contrast=100)
           dsv_nonpref = queries.param_filter_query(self.datastore,st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_orientation=1.57079632679,st_stimulating_signal_parameters_contrast=100)
        else:
           dsv_pref = queries.param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGratingA',st_orientation=0,st_contrast=100)
           dsv_nonpref = queries.param_filter_query(self.datastore,st_name='FullfieldDriftingSinusoidalGratingA',st_orientation=numpy.pi/2,st_contrast=100)

        def autolabel(rects,offset=0.16):
            # attach some text labels
            for rect in rects:
                pylab.gca().text(rect.get_x()+0.043,rect.get_y()  + rect.get_height() + abs(pylab.gca().get_ylim()[0] - pylab.gca().get_ylim()[1])*offset,
                        '%.3f' % float(rect.get_height()),ha='center', va='bottom',fontsize=13,rotation=90)


        color_cycle = {
                            'Bl':(0,0,0),
                            'Or':(.9,.6,0),
                            'SB':(.35,.7,.9),
                            'bG': (0,.6,.5),
                            'Ye':(.95,.9,.25),
                            'Bu':(0,.45,.7),
                            'Ve':(.8,.4,0),
                            'rP':(.8,.6,.7),
                        }
        colors = [color_cycle[c] for c in sorted(color_cycle.keys())]

        mean_CC_pref = numpy.mean(param_filter_query(dsv_pref,analysis_algorithm='TrialMean',sheet_name='V1_Exc_L2/3',identifier='PerNeuronValue',value_name='Trial-to-trial-mean(Correlation coefficient(psth (bin=10.0)))',ads_unique=True).get_analysis_result()[0].get_value_by_ids(close_analogs,close_analogs)[numpy.triu_indices(len(close_analogs),1)])
        mean_CC_nonpref = numpy.mean(param_filter_query(dsv_nonpref,analysis_algorithm='TrialMean',sheet_name='V1_Exc_L2/3',identifier='PerNeuronValue',value_name='Trial-to-trial-mean(Correlation coefficient(psth (bin=10.0)))',ads_unique=True).get_analysis_result()[0].get_value_by_ids(close_analogs,close_analogs)[numpy.triu_indices(len(close_analogs),1)])
        std_CC_pref = numpy.std(param_filter_query(dsv_pref,analysis_algorithm='TrialMean',sheet_name='V1_Exc_L2/3',identifier='PerNeuronValue',value_name='Trial-to-trial-mean(Correlation coefficient(psth (bin=10.0)))',ads_unique=True).get_analysis_result()[0].get_value_by_ids(close_analogs,close_analogs)[numpy.triu_indices(len(close_analogs),1)])
        std_CC_nonpref = numpy.std(param_filter_query(dsv_nonpref,analysis_algorithm='TrialMean',sheet_name='V1_Exc_L2/3',identifier='PerNeuronValue',value_name='Trial-to-trial-mean(Correlation coefficient(psth (bin=10.0)))',ads_unique=True).get_analysis_result()[0].get_value_by_ids(close_analogs,close_analogs)[numpy.triu_indices(len(close_analogs),1)])

        pylab.subplot(gs[0,0])
        r1 = pylab.bar(numpy.array([0.2,0.6])+self.parameters.type*0.1,[mean_CC_pref,mean_CC_nonpref],width=0.08,color=colors[self.parameters.type],yerr=[std_CC_pref,std_CC_nonpref],error_kw=dict(ecolor='gray', lw=2, capsize=5, capthick=2))
        pylab.xlim(0,1.0)
        pylab.ylim(0,0.4)

        pylab.ylabel('Corr. Coef.',fontsize=fontsize)
        pylab.yticks([0,0.4])
        pylab.xticks([0.3,0.7],['preferred','orthogonal'])
        for label in pylab.gca().get_xticklabels() + pylab.gca().get_yticklabels():
            label.set_fontsize(fontsize)
        phf.disable_top_right_axis(pylab.gca())
        autolabel(r1,offset=0.3)

        mean_Exc_pref = numpy.squeeze([a.magnitude for a in param_filter_query(dsv_pref,analysis_algorithm='TrialMean',sheet_name='V1_Exc_L2/3',y_axis_name='exc. conductance trial-to-trial mean',ads_unique=True).get_analysis_result()[0].get_asl_by_id(close_analogs)])
        mean_Inh_pref = numpy.squeeze([a.magnitude for a in param_filter_query(dsv_pref,analysis_algorithm='TrialMean',sheet_name='V1_Exc_L2/3',y_axis_name='inh. conductance trial-to-trial mean',ads_unique=True).get_analysis_result()[0].get_asl_by_id(close_analogs)])
        mean_Exc_nonpref = numpy.squeeze([a.magnitude for a in param_filter_query(dsv_nonpref,analysis_algorithm='TrialMean',sheet_name='V1_Exc_L2/3',y_axis_name='exc. conductance trial-to-trial mean',ads_unique=True).get_analysis_result()[0].get_asl_by_id(close_analogs)])
        mean_Inh_nonpref = numpy.squeeze([a.magnitude for a in param_filter_query(dsv_nonpref,analysis_algorithm='TrialMean',sheet_name='V1_Exc_L2/3',y_axis_name='inh. conductance trial-to-trial mean',ads_unique=True).get_analysis_result()[0].get_asl_by_id(close_analogs)])



        mean_Con_Ratio_pref = numpy.mean(numpy.mean(mean_Exc_pref,axis=1)/numpy.mean(mean_Inh_pref,axis=1))
        mean_Con_Ratio_nonpref = numpy.mean(numpy.mean(mean_Exc_nonpref,axis=1)/numpy.mean(mean_Inh_nonpref,axis=1))
        std_Con_Ratio_pref = numpy.std(numpy.mean(mean_Exc_pref,axis=1)/numpy.mean(mean_Inh_pref,axis=1))
        std_Con_Ratio_nonpref = numpy.std(numpy.mean(mean_Exc_nonpref,axis=1)/numpy.mean(mean_Inh_nonpref,axis=1))

        pylab.subplot(gs[1,0])
        r1 = pylab.bar(numpy.array([0.2,0.6])+self.parameters.type*0.1,[mean_Con_Ratio_pref,mean_Con_Ratio_nonpref],width=0.08,yerr=[std_Con_Ratio_pref,std_Con_Ratio_nonpref],color=colors[self.parameters.type],error_kw=dict(ecolor='gray', lw=2, capsize=5, capthick=2))
        pylab.xlim(0,1.0)
        pylab.ylim(0,1.0)
        pylab.ylabel('Mean Exc:Inh cond. ratio',fontsize=fontsize)
        pylab.yticks([0,1.0])
        pylab.xticks([0.3,0.7],['preferred','orthogonal'])
        for label in pylab.gca().get_xticklabels() + pylab.gca().get_yticklabels():
            label.set_fontsize(fontsize)
        phf.disable_top_right_axis(pylab.gca())
        autolabel(r1)



        mean_ToTVar_pref = numpy.mean(param_filter_query(dsv_pref,analysis_algorithm='AnalogSignal_PerNeuronMeanVar',sheet_name='V1_Exc_L2/3',identifier='PerNeuronValue',value_name='PerNeuronMean(Vm (no AP) trial-to-trial variance)',ads_unique=True).get_analysis_result()[0].get_value_by_id(close_analogs))
        mean_ToTVar_nonpref = numpy.mean(param_filter_query(dsv_nonpref,analysis_algorithm='AnalogSignal_PerNeuronMeanVar',sheet_name='V1_Exc_L2/3',identifier='PerNeuronValue',value_name='PerNeuronMean(Vm (no AP) trial-to-trial variance)',ads_unique=True).get_analysis_result()[0].get_value_by_id(close_analogs))
        std_ToTVar_pref = numpy.std(param_filter_query(dsv_pref,analysis_algorithm='AnalogSignal_PerNeuronMeanVar',sheet_name='V1_Exc_L2/3',identifier='PerNeuronValue',value_name='PerNeuronMean(Vm (no AP) trial-to-trial variance)',ads_unique=True).get_analysis_result()[0].get_value_by_id(close_analogs))
        std_ToTVar_nonpref = numpy.std(param_filter_query(dsv_nonpref,analysis_algorithm='AnalogSignal_PerNeuronMeanVar',sheet_name='V1_Exc_L2/3',identifier='PerNeuronValue',value_name='PerNeuronMean(Vm (no AP) trial-to-trial variance)',ads_unique=True).get_analysis_result()[0].get_value_by_id(close_analogs))

        pylab.subplot(gs[0,1])
        r1 = pylab.bar(numpy.array([0.2,0.6])+self.parameters.type*0.1,[mean_ToTVar_pref,mean_ToTVar_nonpref],width=0.08,yerr=[std_ToTVar_pref,std_ToTVar_nonpref],color=colors[self.parameters.type],error_kw=dict(ecolor='gray', lw=2, capsize=5, capthick=2))
        pylab.xlim(0,1.0)
        pylab.ylim(0,25.0)

        pylab.ylabel('Trial-to-trial Vm variance (mV)',fontsize=fontsize)
        pylab.yticks([0,25.0])
        pylab.xticks([0.3,0.7],['preferred','orthogonal'])
        for label in pylab.gca().get_xticklabels() + pylab.gca().get_yticklabels():
            label.set_fontsize(fontsize)
        phf.disable_top_right_axis(pylab.gca())
        autolabel(r1)



        mean_Var_pref = numpy.mean(param_filter_query(dsv_pref,analysis_algorithm='TrialMean',sheet_name='V1_Exc_L2/3',identifier='PerNeuronValue',value_name='Trial-to-trial-mean(PerNeuronVar(Vm (no AP)))',ads_unique=True).get_analysis_result()[0].get_value_by_id(close_analogs))
        mean_Var_nonpref = numpy.mean(param_filter_query(dsv_nonpref,analysis_algorithm='TrialMean',sheet_name='V1_Exc_L2/3',identifier='PerNeuronValue',value_name='Trial-to-trial-mean(PerNeuronVar(Vm (no AP)))',ads_unique=True).get_analysis_result()[0].get_value_by_id(close_analogs))
        std_Var_pref = numpy.std(param_filter_query(dsv_pref,analysis_algorithm='TrialMean',sheet_name='V1_Exc_L2/3',identifier='PerNeuronValue',value_name='Trial-to-trial-mean(PerNeuronVar(Vm (no AP)))',ads_unique=True).get_analysis_result()[0].get_value_by_id(close_analogs))
        std_Var_nonpref = numpy.std(param_filter_query(dsv_nonpref,analysis_algorithm='TrialMean',sheet_name='V1_Exc_L2/3',identifier='PerNeuronValue',value_name='Trial-to-trial-mean(PerNeuronVar(Vm (no AP)))',ads_unique=True).get_analysis_result()[0].get_value_by_id(close_analogs))

        pylab.subplot(gs[1,1])
        r1 = pylab.bar(numpy.array([0.2,0.6])+self.parameters.type*0.1,[mean_Var_pref,mean_Var_nonpref],width=0.08,yerr=[std_Var_pref,std_Var_nonpref],color=colors[self.parameters.type],error_kw=dict(ecolor='gray', lw=2, capsize=5, capthick=2))
        pylab.xlim(0,1.0)
        pylab.ylim(0,40.0)

        pylab.ylabel('Trial averaged Vm variance (mV)',fontsize=fontsize)
        pylab.yticks([0,40.0])
        pylab.xticks([0.3,0.7],['preferred','orthogonal'])
        for label in pylab.gca().get_xticklabels() + pylab.gca().get_yticklabels():
            label.set_fontsize(fontsize)
        phf.disable_top_right_axis(pylab.gca())
        autolabel(r1)

        if False:
            logger.info('DSADSADFSAFSDAFSDFA ' + str(close_analogs))
            vm_cc_gr_s1,psth_cc_gr_s1 = self.calculate_cc('V1_Exc_L2/3',close_analogs)
            logger.info('DSADSADFSAFSDAFSDFA ' + str(len(vm_cc_gr_s1)) + '  ' + str(len(psth_cc_gr_s1)))
            z = int(min(self.parameters.window_length,(len(vm_cc_gr_s1)-1)/2)/2)*2
            logger.info('DSADSADFSAFSDAFSDFA ' + str(z))
            bin_size=10.0
            plots["Spike_sheet_1"] = (StandardStyleLinePlot([numpy.linspace(-z,z,2*int(z/bin_size)+1)], [psth_cc_gr_s1[int(len(psth_cc_gr_s1)/2)-int(z/bin_size):int(len(psth_cc_gr_s1)/2)+int(z/bin_size)+1]]),gs[0,2],{'colors': [colors[self.parameters.type]], 'x_tick_style' : 'Custom', 'x_ticks' : [],'y_tick_style' : 'Custom', 'y_ticks' : [0,0.2], 'y_tick_labels' : [0.0,0.1], 'linewidth' : 2.0, 'y_lim' : (-0.02,0.2),'y_label' : 'spikes'})
            plots["Vm_sheet_1"] = (StandardStyleLinePlot([numpy.linspace(-z,z,2*z+1)], [vm_cc_gr_s1.flatten()[int(len(vm_cc_gr_s1)/2)-z:int(len(vm_cc_gr_s1)/2)+z+1]]),gs[1,2],{'x_label' : 'time(ms)', 'colors':[colors[self.parameters.type]], 'x_tick_style' : 'Custom', 'x_ticks' : [-z,0,z], 'x_tick_labels' : [-self.parameters.window_length,0,self.parameters.window_length],'y_tick_style' : 'Custom', 'y_ticks' : [-0.05,0,0.2], 'y_tick_labels' : [-0.05,0.0,0.2], 'linewidth' : 2.0, 'y_lim' : (-0.05,0.1),'y_label' : 'Vm'})


        x = param_filter_query(dsv_pref,analysis_algorithm='TrialMean',sheet_name='V1_Exc_L2/3',identifier='PerNeuronValue',value_name='Trial-to-trial-mean(PerNeuronVar(Vm (no AP)))',ads_unique=True).get_analysis_result()[0].get_value_by_id(close_analogs)
        y = param_filter_query(dsv_pref,sheet_name='V1_Exc_L2/3',value_name=['Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(close_analogs)

        pylab.subplot(gs[0,2])
        pylab.plot(x,y,'o',color=colors[self.parameters.type])
        pylab.xlim(8,22.0)
        pylab.ylim(0,20.0)
        pylab.xticks([10,15,20])
        pylab.yticks([0,10,20])
        pylab.xlabel('Vm variance (mV)',fontsize=fontsize)
        pylab.ylabel('firing rate (sp/s)',fontsize=fontsize)

        for label in pylab.gca().get_xticklabels() + pylab.gca().get_yticklabels():
            label.set_fontsize(fontsize)
        phf.disable_top_right_axis(pylab.gca())
        

        x = param_filter_query(dsv_pref,analysis_algorithm='Analog_MeanSTDAndFanoFactor',sheet_name='V1_Exc_L2/3',identifier='PerNeuronValue',value_name='Mean(VM)',ads_unique=True).get_analysis_result()[0].get_value_by_id(close_analogs)
        #x = param_filter_query(dsv_pref,analysis_algorithm='SubtractPNVfromPNVS',sheet_name='V1_Exc_L2/3',identifier='PerNeuronValue',value_name='Mean(VM)-Mean(VM)',ads_unique=True).get_analysis_result()[0].get_value_by_id(close_analogs)
        y = param_filter_query(dsv_pref,sheet_name='V1_Exc_L2/3',value_name=['Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(close_analogs)

        pylab.subplot(gs[1,2])
        pylab.plot(x,y,'o',color=colors[self.parameters.type])
        pylab.xlim(-61.5,-58.5)
        pylab.ylim(0,20.0)
        pylab.xticks([-61,-58.0])
        pylab.yticks([0,20,20])
        pylab.xlabel('Vm mean (mV)',fontsize=fontsize)
        pylab.ylabel('firing rate (sp/s)',fontsize=fontsize)

        for label in pylab.gca().get_xticklabels() + pylab.gca().get_yticklabels():
            label.set_fontsize(fontsize)
        phf.disable_top_right_axis(pylab.gca())
        
        return plots


class ContrastResponse(Plotting):
    required_parameters = ParameterSet({
        'cortical_stimulation' : int,
    })

    def get_vals(self,dsv,neuron,variable,mean=False):
        assert queries.ads_with_equal_stimulus_type(dsv)
        assert queries.equal_ads(dsv,except_params=['stimulus_id'])
        pnvs = dsv.get_analysis_result()

        st = [MozaikParametrized.idd(s.stimulus_id) for s in pnvs]
        if mean:
            tc_dict = colapse_to_dictionary([numpy.mean(z.get_value_by_id(neuron)) for z in pnvs],st,variable)
        else:
            tc_dict = colapse_to_dictionary([z.get_value_by_id(neuron) for z in pnvs],st,variable)

        rads = tc_dict.values()[0][0]
        values = tc_dict.values()[0][1]
        a, b = zip(*sorted(zip(rads,values)))
        return numpy.array(a),numpy.array(b)

    def subplot(self, subplotspec):
        plots = {}
        gs = gridspec.GridSpecFromSubplotSpec(2,2, subplot_spec=subplotspec,hspace=0.25, wspace=0.3)

        scale=0.16

        fontsize=15
        spike_ids = param_filter_query(self.datastore,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_spike_train_ids()
        l23_exc_or = self.datastore.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L2/3')[0]
        l23_exc_or_many = numpy.array(spike_ids)[numpy.nonzero(numpy.array([circular_dist(l23_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids]) < 0.1)[0]]

        dsv = queries.param_filter_query(self.datastore,sheet_name='V1_Exc_L2/3',st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_orientation=0,analysis_algorithm='TrialAveragedFiringRate',value_name='Firing rate')
        dsv_var = queries.param_filter_query(self.datastore,sheet_name='V1_Exc_L2/3',st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_orientation=0,analysis_algorithm='TrialAveragedFiringRate',value_name='Tria-to-trial Var of Firing rate')
        variable = 'stimulating_signal_parameters_scale'  

        fitfunc = lambda p,x: p[1]*numpy.power(x,p[0])/(numpy.power(x,p[0])+p[2])

        x,y = self.get_vals(dsv,l23_exc_or_many[0],variable)
        x,yvar = self.get_vals(dsv_var,l23_exc_or_many[0],variable)
        x = numpy.concatenate(([0],x))
        y = numpy.concatenate(([0],y))
        yvar = numpy.concatenate(([0],yvar))
        
        c50 = param_filter_query(self.datastore,analysis_algorithm='NakaRushtonTuningCurveFit',value_name='stimulating_signal_parameters_scale Naka-Rushton c50 of Firing rate').get_analysis_result()[0].get_value_by_id(l23_exc_or_many[0])
        scale = param_filter_query(self.datastore,analysis_algorithm='NakaRushtonTuningCurveFit',value_name='stimulating_signal_parameters_scale Naka-Rushton scaler of Firing rate').get_analysis_result()[0].get_value_by_id(l23_exc_or_many[0])
        exponent = param_filter_query(self.datastore,analysis_algorithm='NakaRushtonTuningCurveFit',value_name='stimulating_signal_parameters_scale Naka-Rushton exponent of Firing rate').get_analysis_result()[0].get_value_by_id(l23_exc_or_many[0])
        yfit = fitfunc([exponent,scale,c50],numpy.arange(0,numpy.max(x),0.01))

        ax = pylab.subplot(gs[0,0])
        ax.errorbar(x,y,yerr=yvar,fmt='o',color='k')         
        ax.plot(numpy.arange(0,numpy.max(x),0.01),yfit,color='gray')   
        disable_top_right_axis(pylab.gca())  
        three_tick_axis(ax.yaxis)  
        three_tick_axis(ax.xaxis)
        pylab.ylim(ymin=0,ymax=int(math.ceil(numpy.max(y+yvar))+1)/2*2)
        pylab.xlim(0,0.33)
        pylab.xticks([0,0.33],["0",r"$0.33 \times 10^{16}$"])

        for label in pylab.gca().get_xticklabels() + pylab.gca().get_yticklabels():
                label.set_fontsize(fontsize)
        pylab.ylabel('spikes/s',fontsize=fontsize)
        remove_x_tick_labels()


        x,y = self.get_vals(dsv,l23_exc_or_many[3],variable)
        x,yvar = self.get_vals(dsv_var,l23_exc_or_many[3],variable)
        x = numpy.concatenate(([0],x))
        y = numpy.concatenate(([0],y))
        yvar = numpy.concatenate(([0],yvar))

        c50 = param_filter_query(self.datastore,analysis_algorithm='NakaRushtonTuningCurveFit',value_name='stimulating_signal_parameters_scale Naka-Rushton c50 of Firing rate').get_analysis_result()[0].get_value_by_id(l23_exc_or_many[3])
        scale = param_filter_query(self.datastore,analysis_algorithm='NakaRushtonTuningCurveFit',value_name='stimulating_signal_parameters_scale Naka-Rushton scaler of Firing rate').get_analysis_result()[0].get_value_by_id(l23_exc_or_many[3])
        exponent = param_filter_query(self.datastore,analysis_algorithm='NakaRushtonTuningCurveFit',value_name='stimulating_signal_parameters_scale Naka-Rushton exponent of Firing rate').get_analysis_result()[0].get_value_by_id(l23_exc_or_many[3])
        yfit = fitfunc([exponent,scale,c50],numpy.arange(0,numpy.max(x),0.01))


        ax = pylab.subplot(gs[1,0])
        ax.errorbar(x,y,yerr=yvar,fmt='o',color='k')      
        ax.plot(numpy.arange(0,numpy.max(x),0.01),yfit,color='gray')            
        disable_top_right_axis(pylab.gca())  
        three_tick_axis(ax.yaxis)  
        three_tick_axis(ax.xaxis)
        pylab.ylim(ymin=0,ymax=int(math.ceil(numpy.max(y+yvar))+1)/2*2)
        pylab.xlim(0,0.33)
        pylab.xticks([0,0.33],["0",r"$0.33 \times 10^{16}$"])
        
        for label in pylab.gca().get_xticklabels() + pylab.gca().get_yticklabels():
                label.set_fontsize(fontsize)
        pylab.xlabel('scale',fontsize=fontsize)
        pylab.ylabel('spikes/s',fontsize=fontsize)


        x,y = self.get_vals(dsv,l23_exc_or_many,variable,mean=True)
        x,yvar = self.get_vals(dsv_var,l23_exc_or_many,variable)
        x = numpy.concatenate(([0],x))
        y = numpy.concatenate(([0],y))
        

        c50 = param_filter_query(self.datastore,analysis_algorithm='NakaRushtonTuningCurveFit',value_name='stimulating_signal_parameters_scale Mean Naka-Rushton c50 of Firing rate').get_analysis_result()[0].value
        scale = param_filter_query(self.datastore,analysis_algorithm='NakaRushtonTuningCurveFit',value_name='stimulating_signal_parameters_scale Mean Naka-Rushton scaler of Firing rate').get_analysis_result()[0].value
        exponent = param_filter_query(self.datastore,analysis_algorithm='NakaRushtonTuningCurveFit',value_name='stimulating_signal_parameters_scale Mean Naka-Rushton exponent of Firing rate').get_analysis_result()[0].value
        yfit = fitfunc([exponent,scale,c50],numpy.arange(0,numpy.max(x),0.01))


        ax = pylab.subplot(gs[:,1])
        ax.plot(x,y,'ok')      
        ax.errorbar(x,y,yerr=yvar,fmt='o',color='k')      
        ax.plot(numpy.arange(0,numpy.max(x),0.01),yfit,color='gray')            
        disable_top_right_axis(pylab.gca())  
        three_tick_axis(ax.yaxis)  
        three_tick_axis(ax.xaxis)
        pylab.ylim(ymin=0,ymax=int(math.ceil(numpy.max(y))+1)/2*2)
        pylab.xlim(0,0.33)
        pylab.xticks([0,0.33],["0",r"$0.33 \times 10^{16}$"])
        
        for label in pylab.gca().get_xticklabels() + pylab.gca().get_yticklabels():
                label.set_fontsize(fontsize)
        pylab.xlabel('scale',fontsize=fontsize)
        pylab.ylabel('spikes/s',fontsize=fontsize)

        return plots


class ContrastResponseTransformed(Plotting):
    required_parameters = ParameterSet({
        'cortical_stimulation' : int,
        'type' : int,
    })

    def get_vals(self,dsv,neuron,variable,mean=False):
        assert queries.ads_with_equal_stimulus_type(dsv)
        assert queries.equal_ads(dsv,except_params=['stimulus_id'])
        pnvs = dsv.get_analysis_result()

        st = [MozaikParametrized.idd(s.stimulus_id) for s in pnvs]
        if mean:
            tc_dict = colapse_to_dictionary([numpy.mean(z.get_value_by_id(neuron)) for z in pnvs],st,variable)
        else:
            tc_dict = colapse_to_dictionary([z.get_value_by_id(neuron) for z in pnvs],st,variable)

        rads = tc_dict.values()[0][0]
        values = tc_dict.values()[0][1]
        a, b = zip(*sorted(zip(rads,values)))
        return numpy.array(a),numpy.array(b)

    def subplot(self, subplotspec):

        color_cycle = {
                            'Bl':(0,0,0),
                            'Or':(.9,.6,0),
                            'SB':(.35,.7,.9),
                            'bG': (0,.6,.5),
                            'Ye':(.95,.9,.25),
                            'Bu':(0,.45,.7),
                            'Ve':(.8,.4,0),
                            'rP':(.8,.6,.7),
                        }
        colors = [color_cycle[c] for c in sorted(color_cycle.keys())]
        
        plots = {}
        gs = gridspec.GridSpecFromSubplotSpec(1,1, subplot_spec=subplotspec,hspace=0.25, wspace=0.3)
        
        spike_ids = param_filter_query(self.datastore,sheet_name="V1_Exc_L2/3").get_segments()[0].get_stored_spike_train_ids()
        l23_exc_or = self.datastore.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L2/3')[0]
        l23_exc_or_many = numpy.array(spike_ids)[numpy.nonzero(numpy.array([circular_dist(l23_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids]) < 0.1)[0]]

        if self.parameters.cortical_stimulation:
            dsv = queries.param_filter_query(self.datastore,sheet_name='V1_Exc_L2/3',st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_orientation=0,analysis_algorithm='TrialAveragedFiringRate',value_name='Firing rate')
            dsv_var = queries.param_filter_query(self.datastore,sheet_name='V1_Exc_L2/3',st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_orientation=0,analysis_algorithm='TrialAveragedFiringRate',value_name='Tria-to-trial Var of Firing rate')
            variable = 'stimulating_signal_parameters_contrast'  
        else:
            dsv = queries.param_filter_query(self.datastore,sheet_name='V1_Exc_L2/3',st_name='FullfieldDriftingSinusoidalGratingA',st_orientation=0,analysis_algorithm='TrialAveragedFiringRate',value_name='Firing rate')
            dsv_var = queries.param_filter_query(self.datastore,sheet_name='V1_Exc_L2/3',st_name='FullfieldDriftingSinusoidalGratingA',st_orientation=0,analysis_algorithm='TrialAveragedFiringRate',value_name='Tria-to-trial Var of  Firing rate')
            variable = 'contrast'  

        fontsize = 15
       
        fitfunc = lambda p,x: p[1]*numpy.power(x,p[0])/(numpy.power(x,p[0])+p[2])

        x,y = self.get_vals(dsv,l23_exc_or_many,variable,mean=True)
        c50 = param_filter_query(self.datastore,analysis_algorithm='NakaRushtonTuningCurveFit',value_name=variable + ' Mean Naka-Rushton c50 of Firing rate').get_analysis_result()[0].value
        scale = param_filter_query(self.datastore,analysis_algorithm='NakaRushtonTuningCurveFit',value_name=variable + ' Mean Naka-Rushton scaler of Firing rate').get_analysis_result()[0].value
        exponent = param_filter_query(self.datastore,analysis_algorithm='NakaRushtonTuningCurveFit',value_name=variable + ' Mean Naka-Rushton exponent of Firing rate').get_analysis_result()[0].value
        yfit = fitfunc([exponent,scale,c50],numpy.arange(0,100,1))

        ax = pylab.subplot(gs[0,0])
        ax.plot(x,y,'o',color=colors[self.parameters.type])         
        ax.plot(numpy.arange(0,100,1),yfit,color=colors[self.parameters.type])         
        disable_top_right_axis(pylab.gca())  
        three_tick_axis(ax.yaxis)  
        three_tick_axis(ax.xaxis)
        pylab.ylim(0,10)
        pylab.xlim(0,100)
        
        for label in pylab.gca().get_xticklabels() + pylab.gca().get_yticklabels():
                label.set_fontsize(19)
        pylab.ylabel('spikes/s',fontsize=fontsize)
        #pylab.xlabel('contrast',fontsize=fontsize)
        remove_x_tick_labels()

        return plots

class OrientationTuningSharpness(Plotting):
    required_parameters = ParameterSet({
        })

    def subplot(self, subplotspec):
        plots = {}
        gs = gridspec.GridSpecFromSubplotSpec(4, 6, subplot_spec=subplotspec,
                                              hspace=0.5    , wspace=0.1) 
        
        spike_ids = sorted(numpy.random.permutation(queries.param_filter_query(self.datastore,sheet_name='V1_Exc_L2/3').get_segments()[0].get_stored_spike_train_ids()))
        l23_exc_or = self.datastore.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L2/3')[0]
        l23_exc_or_many = numpy.array(spike_ids)[numpy.nonzero(numpy.array([circular_dist(l23_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids]) < 0.1)[0]]


        contrasts = [3,100]
        sharpness = [0.1,0.3,0.5,1.0,2.0,3.0]
        fontsize = 15
        
        for i in xrange(0,6):    

            #lamsChR = self.datastore.get_analysis_result(identifier='PerNeuronValue',analysis_algorithm ='NeuronAnnotationsToPerNeuronValues', value_name = 'Light activation magnitude ChR(V1_Exc_L2/3,0.00616172489758,0.0,'+ str(sharpness[i]) +',10)', sheet_name = 'V1_Exc_L2/3')[0]
            queries.param_filter_query(self.datastore,identifier='PerNeuronValue',analysis_algorithm ='NeuronAnnotationsToPerNeuronValues',sheet_name = 'V1_Exc_L2/3').print_content(full_ADS=True)
            lams = self.datastore.get_analysis_result(identifier='PerNeuronValue',analysis_algorithm ='NeuronAnnotationsToPerNeuronValues', value_name = 'Light activation ma13gnitude(V1_Exc_L2/3,0.0061617248975781415,0.0,'+ str(sharpness[i]) +',10)', sheet_name = 'V1_Exc_L2/3')[0]
 
            axis = pylab.subplot(gs[0,i])
            pylab.scatter((l23_exc_or.values+numpy.pi/2)%numpy.pi-numpy.pi/2,lams.values,s=3)
            pylab.xticks([-numpy.pi/2, 0, numpy.pi/2],["-$\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$"])
            if i==0:
                pylab.ylabel(r'$photons/cm^2/s$',fontsize=fontsize)
                pylab.xlabel('orientation',fontsize=fontsize)
                pylab.yticks([0,5e14],[0.0,r'$5 {\times} 10^{14}$'],fontsize=fontsize)
            else:
                phf.disable_left_axis(axis)
                pylab.yticks([])
            
            pylab.ylim(0,5.0e14)
            pylab.xlim([-numpy.pi/2, numpy.pi/2])
            phf.disable_top_right_axis(axis)
            for label in axis.get_xticklabels() + axis.get_yticklabels():
                label.set_fontsize(fontsize)
            
            pylab.title('$\\sigma = ' + str(sharpness[i]) + '$',fontsize=fontsize)


            base = queries.param_filter_query(self.datastore,sheet_name='V1_Exc_L2/3',st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_sharpness=sharpness[i],st_stimulating_signal_parameters_contrast = 100,value_name=['stimulating_signal_parameters_orientation baseline of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(spike_ids)
            mmax = queries.param_filter_query(self.datastore,sheet_name='V1_Exc_L2/3',st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_stimulating_signal_parameters_sharpness=sharpness[i],st_stimulating_signal_parameters_contrast = 100,value_name=['stimulating_signal_parameters_orientation max of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(spike_ids)                    
            responsive_spike_ids = numpy.array(spike_ids)[numpy.array(base)+numpy.array(mmax) > 1.0]
                    
            dsv = queries.param_filter_query(self.datastore,st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',value_name='Firing rate',analysis_algorithm=['TrialAveragedFiringRate'],st_stimulating_signal_parameters_contrast = contrasts,st_stimulating_signal_parameters_sharpness = sharpness[i])
            if i == 0:
                plots['ORTCMean'+str(i)] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'stimulating_signal_parameters_orientation', 'neurons': list(responsive_spike_ids), 'sheet_name' : "V1_Exc_L2/3",'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[1,i],{'x_label' : 'orientation','y_lim' : (0,10), 'y_label' : 'firing rate (sp/s)', 'title' : None})
            else:
                plots['ORTCMean'+str(i)] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'stimulating_signal_parameters_orientation', 'neurons': list(responsive_spike_ids), 'sheet_name' : "V1_Exc_L2/3",'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[1,i],{'x_label' : 'orientation','y_lim' : (0,10),'title' : None, 'x_label' : None, 'left_border' : False, 'y_axis' : False, 'title' : None})
            

            low_c = queries.param_filter_query(self.datastore,value_name=['stimulating_signal_parameters_orientation HWHH of Firing rate'],sheet_name=['V1_Exc_L2/3'],st_stimulating_signal_parameters_contrast=contrasts[0],st_stimulating_signal_parameters_sharpness=sharpness[i]).get_analysis_result()[0].get_value_by_id(responsive_spike_ids)
            high_c = queries.param_filter_query(self.datastore,value_name=['stimulating_signal_parameters_orientation HWHH of Firing rate'],sheet_name=['V1_Exc_L2/3'],st_stimulating_signal_parameters_contrast=contrasts[1],st_stimulating_signal_parameters_sharpness=sharpness[i]).get_analysis_result()[0].get_value_by_id(responsive_spike_ids)

            if True:
                axis = pylab.subplot(gs[3,i])
                logger.info(low_c)
                logger.info(high_c)
                pylab.plot([0,60],[0,60],color='#AAAAAA')
                pylab.plot(low_c,high_c,'ok')
                pylab.xlim(0,60)
                pylab.ylim(0,60)
                if i == 0:
                    pylab.xlabel('HWHH low cont.',fontsize=fontsize)
                    pylab.ylabel('HWHH high cont.',fontsize=fontsize)
                    pylab.yticks([0,30,60])
                else:
                    phf.disable_left_axis(axis)
                    pylab.yticks([])

                pylab.xticks([0,30,60])
                
                
                phf.disable_top_right_axis(axis)
                for label in axis.get_xticklabels() + axis.get_yticklabels():
                    label.set_fontsize(fontsize)

            dsv = queries.param_filter_query(self.datastore,value_name=['stimulating_signal_parameters_orientation HWHH of Firing rate'],sheet_name=['V1_Exc_L2/3'],st_stimulating_signal_parameters_contrast=100,st_stimulating_signal_parameters_sharpness = sharpness[i])    
            if i == 0:
                plots['HWHHHistogramExc'+str(i)] = (PerNeuronValuePlot(dsv, ParameterSet({'cortical_view' : False})),gs[2,i],{'mark_mean' : True, 'x_lim' : (0.0,50.0), 'y_lim' : (0.0,400.0), 'title' : None,'y_label' : '# neurons', 'x_label' : 'HWHH', 'mark_value' : 27})
            else:
                plots['HWHHHistogramExc'+str(i)] = (PerNeuronValuePlot(dsv, ParameterSet({'cortical_view' : False})),gs[2,i],{'mark_mean' : True, 'x_lim' : (0.0,50.0),'y_lim' : (0.0,400.0), 'x_label' : None,'title' : None,'y_label' : '# neurons','left_border' : False, 'y_axis' : False, 'mark_value' : 27})

        return plots



class OrientationTuningLEE_size(Plotting):
    required_parameters = ParameterSet({
        })

    def subplot(self, subplotspec):
        plots = {}
        gs = gridspec.GridSpecFromSubplotSpec(4, 5, subplot_spec=subplotspec,
                                              hspace=0.5, wspace=0.1)
        
        spike_ids = sorted(numpy.random.permutation(queries.param_filter_query(self.datastore,sheet_name='V1_Exc_L2/3').get_segments()[0].get_stored_spike_train_ids()))
        l23_exc_or = self.datastore.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L2/3')[0]
        l23_exc_or_many = numpy.array(spike_ids)[numpy.nonzero(numpy.array([circular_dist(l23_exc_or.get_value_by_id(i),0,numpy.pi)  for i in spike_ids]) < 0.1)[0]]


        contrasts = [3,100]
        sizes = [10,20,50,100,300]
        fontsize = 15
        
        for i in xrange(0,5):    

            if True:
                lams = self.datastore.get_analysis_result(identifier='PerNeuronValue',analysis_algorithm ='NeuronAnnotationsToPerNeuronValues', value_name = 'Light activation magnitude(V1_Exc_L2/3,0.00616172489758,0.0,0.5,'+ str(sizes[i]) +')', sheet_name = 'V1_Exc_L2/3')[0]
     
                axis = pylab.subplot(gs[0,i])
                pylab.scatter((l23_exc_or.values+numpy.pi/2)%numpy.pi-numpy.pi/2,lams.values,s=3)
                pylab.xticks([-numpy.pi/2, 0, numpy.pi/2],["-$\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$"])
                if i==0:
                    pylab.ylabel(r'$photons/cm^2/s$',fontsize=fontsize)
                    pylab.xlabel('orientation',fontsize=fontsize)
                    pylab.yticks([0,6e14],[0.0,r'$6 {\times} 10^{14}$'],fontsize=fontsize)
                else:
                    phf.disable_left_axis(axis)
                    pylab.yticks([])

                pylab.ylim(0,6.0e14)
                pylab.xlim([-numpy.pi/2, numpy.pi/2])
                pylab.title('$\\sigma = ' + str(sizes[i]) + '$')
                phf.disable_top_right_axis(axis)
                for label in axis.get_xticklabels() + axis.get_yticklabels():
                    label.set_fontsize(fontsize)


            base = queries.param_filter_query(self.datastore,sheet_name='V1_Exc_L2/3',st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_spacing=sizes[i],st_stimulating_signal_parameters_contrast = 100,value_name=['stimulating_signal_parameters_orientation baseline of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(spike_ids)
            mmax = queries.param_filter_query(self.datastore,sheet_name='V1_Exc_L2/3',st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',st_spacing=sizes[i],st_stimulating_signal_parameters_contrast = 100,value_name=['stimulating_signal_parameters_orientation max of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(spike_ids)                    
            responsive_spike_ids = numpy.array(spike_ids)[numpy.array(base)+numpy.array(mmax) > 1.0]
                    
            dsv = queries.param_filter_query(self.datastore,st_name='InternalStimulus',st_direct_stimulation_name='LocalStimulatorArray',value_name='Firing rate',analysis_algorithm=['TrialAveragedFiringRate'],st_stimulating_signal_parameters_contrast = contrasts,st_spacing = sizes[i])
            if i == 0:
                plots['ORTCMean'+str(i)] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'stimulating_signal_parameters_orientation', 'neurons': list(responsive_spike_ids), 'sheet_name' : "V1_Exc_L2/3",'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[1,i],{'x_label' : 'orientation','y_lim' : (0,12), 'y_label' : 'firing rate (sp/s)', 'title' : None})
            else:
                plots['ORTCMean'+str(i)] = (PlotTuningCurve(dsv, ParameterSet({'parameter_name' : 'stimulating_signal_parameters_orientation', 'neurons': list(responsive_spike_ids), 'sheet_name' : "V1_Exc_L2/3",'centered'  : True,'mean' : True,'pool' : False,'polar' : False})),gs[1,i],{'x_label' : 'orientation','y_lim' : (0,12),'title' : None, 'x_label' : None, 'left_border' : False, 'y_axis' : False, 'title' : None})
            

            low_c = queries.param_filter_query(self.datastore,value_name=['stimulating_signal_parameters_orientation HWHH of Firing rate'],sheet_name=['V1_Exc_L2/3'],st_stimulating_signal_parameters_contrast=contrasts[0],st_spacing = sizes[i]).get_analysis_result()[0].get_value_by_id(responsive_spike_ids)
            high_c = queries.param_filter_query(self.datastore,value_name=['stimulating_signal_parameters_orientation HWHH of Firing rate'],sheet_name=['V1_Exc_L2/3'],st_stimulating_signal_parameters_contrast=contrasts[1],st_spacing = sizes[i]).get_analysis_result()[0].get_value_by_id(responsive_spike_ids)

            axis = pylab.subplot(gs[2,i])
            pylab.plot([0,60],[0,60],color='#AAAAAA')
            pylab.plot(low_c,high_c,'ok')
            pylab.xlim(0,60)
            pylab.xticks([0,30,60])
            pylab.ylim(0,60)
            if i == 0:
                pylab.xlabel('HWHH low contrast',fontsize=fontsize)
                pylab.yticks([0,30,60])
                pylab.ylabel('HWHH  high contrast',fontsize=fontsize)
            else:
                phf.disable_left_axis(axis)
                pylab.yticks([])
            for label in axis.get_xticklabels() + axis.get_yticklabels():
                label.set_fontsize(fontsize)
            phf.disable_top_right_axis(axis)

            dsv = queries.param_filter_query(self.datastore,value_name=['stimulating_signal_parameters_orientation HWHH of Firing rate'],sheet_name=['V1_Exc_L2/3'],st_stimulating_signal_parameters_contrast=100,st_spacing = sizes[i])    
            if i == 0:
                plots['HWHHHistogramExc'+str(i)] = (PerNeuronValuePlot(dsv, ParameterSet({'cortical_view' : False})),gs[3,i],{'mark_mean' : True, 'x_lim' : (0.0,50.0), 'y_lim' : (0.0,400.0), 'title' : None,'y_label' : '# neurons', 'x_label' : 'HWHH', 'mark_value' : 27})
            else:
                plots['HWHHHistogramExc'+str(i)] = (PerNeuronValuePlot(dsv, ParameterSet({'cortical_view' : False})),gs[3,i],{'mark_mean' : True, 'x_lim' : (0.0,50.0),'y_lim' : (0.0,400.0), 'x_label' : None,'title' : None,'y_label' : '# neurons','left_border' : False, 'y_axis' : False, 'mark_value' : 27})

        return plots

