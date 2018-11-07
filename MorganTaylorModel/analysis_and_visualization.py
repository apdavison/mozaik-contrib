import mozaik
from mozaik.visualization.plotting import *
from mozaik.analysis.technical import NeuronAnnotationsToPerNeuronValues
from mozaik.analysis.analysis import *
from mozaik.analysis.vision import *
from mozaik.storage.queries import *
from mozaik.storage.datastore import PickledDataStore
from mozaik.controller import Global

import sys
sys.path.append('/home/jantolik/projects/mozaik/contrib')
import Kremkow_plots
from Kremkow_plots import *
from lsv1m_paper import *


logger = mozaik.getMozaikLogger()

import os


class AAA(Analysis):
      """
      """
      def perform_analysis(self):

	    spike_ids = param_filter_query(self.datastore,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
	    spike_ids_inh = param_filter_query(self.datastore,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()

            mmax5 = param_filter_query(self.datastore,sheet_name='V1_Exc_L4',st_direct_stimulation_name="None",st_name=['FullfieldDriftingSinusoidalGrating'],st_contrast=5,value_name=['orientation max of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(spike_ids)
            mmax100 = param_filter_query(self.datastore,sheet_name='V1_Exc_L4',st_direct_stimulation_name="None",st_name=['FullfieldDriftingSinusoidalGrating'],st_contrast=100,value_name=['orientation max of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(spike_ids)
            responsive_spike_ids_l4E = numpy.array(spike_ids)[numpy.logical_and(numpy.array(mmax5) > 2.0,numpy.array(mmax100) > 2.0)]

	    self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.mean(mmax5),period=None,value_name = 'Tuning Height LC',sheet_name='V1_Exc_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        
	    self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.mean(mmax100),period=None,value_name = 'Tuning Height HC',sheet_name='V1_Exc_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        

            mmax5 = param_filter_query(self.datastore,sheet_name='V1_Inh_L4',st_direct_stimulation_name="None",st_name=['FullfieldDriftingSinusoidalGrating'],st_contrast=5,value_name=['orientation max of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(spike_ids_inh)
            mmax100 = param_filter_query(self.datastore,sheet_name='V1_Inh_L4',st_direct_stimulation_name="None",st_name=['FullfieldDriftingSinusoidalGrating'],st_contrast=100,value_name=['orientation max of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(spike_ids_inh)
            responsive_spike_ids_l4I = numpy.array(spike_ids_inh)[numpy.logical_and(numpy.array(mmax5) > 2.0,numpy.array(mmax100) > 2.0)]

	    self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.mean(mmax5),period=None,value_name = 'Tuning Height LC',sheet_name='V1_Inh_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        
	    self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.mean(mmax100),period=None,value_name = 'Tuning Height HC',sheet_name='V1_Inh_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        


            hwhh_hc = numpy.array(param_filter_query(self.datastore,sheet_name='V1_Exc_L4',st_direct_stimulation_name="None",st_name=['FullfieldDriftingSinusoidalGrating'],st_contrast=100,value_name=['orientation HWHH of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(responsive_spike_ids_l4E))
            hwhh_lc = numpy.array(param_filter_query(self.datastore,sheet_name='V1_Exc_L4',st_direct_stimulation_name="None",st_name=['FullfieldDriftingSinusoidalGrating'],st_contrast=5,value_name=['orientation HWHH of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(responsive_spike_ids_l4E))

            self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.mean(hwhh_hc[hwhh_hc<1000]),period=None,value_name = 'Mean HWHH of responsive neurons',sheet_name='V1_Exc_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        
            self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.mean(abs(hwhh_hc-hwhh_lc)[abs(hwhh_hc-hwhh_lc)<1000]),period=None,value_name = 'Mean HWHH difference',sheet_name='V1_Exc_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        
	    print numpy.mean(abs(hwhh_hc-hwhh_lc)[abs(hwhh_hc-hwhh_lc)<1000])
	    self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.sum(numpy.logical_or(hwhh_hc<=0,hwhh_hc>=100)),period=None,value_name = 'Faild fits HC',sheet_name='V1_Exc_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        
	    self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.sum(numpy.logical_or(hwhh_lc<=0,hwhh_lc>=100)),period=None,value_name = 'Faild fits LC',sheet_name='V1_Exc_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        
	    
            hwhh_hc = numpy.array(param_filter_query(self.datastore,sheet_name='V1_Inh_L4',st_direct_stimulation_name="None",st_name=['FullfieldDriftingSinusoidalGrating'],st_contrast=100,value_name=['orientation HWHH of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(responsive_spike_ids_l4I))
            hwhh_lc = numpy.array(param_filter_query(self.datastore,sheet_name='V1_Inh_L4',st_direct_stimulation_name="None",st_name=['FullfieldDriftingSinusoidalGrating'],st_contrast=5,value_name=['orientation HWHH of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(responsive_spike_ids_l4I))

            self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.mean(hwhh_hc[hwhh_hc<1000]),period=None,value_name = 'Mean HWHH of responsive neurons',sheet_name='V1_Inh_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        
            self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.mean(abs(hwhh_hc-hwhh_lc)[abs(hwhh_hc-hwhh_lc)<1000]),period=None,value_name = 'Mean HWHH difference',sheet_name='V1_Inh_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        
	    print numpy.mean(abs(hwhh_hc-hwhh_lc)[abs(hwhh_hc-hwhh_lc)<1000])
	    self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.sum(numpy.logical_or(hwhh_hc<=0,hwhh_hc>=100)),period=None,value_name = 'Faild fits HC',sheet_name='V1_Inh_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        
	    self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.sum(numpy.logical_or(hwhh_lc<=0,hwhh_lc>=100)),period=None,value_name = 'Faild fits LC',sheet_name='V1_Inh_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        


class AAA1(Analysis):
      """
      """
      def perform_analysis(self):
	    TrialMean(self.datastore,ParameterSet({'vm': False,  'cond_exc': True, 'cond_inh': True})).analyse()
	    dsv = param_filter_query(self.datastore,analysis_algorithm='TrialMean',st_name='FullfieldDriftingSinusoidalGrating',sheet_name='V1_Exc_L4',y_axis_name=['inh. conductance trial-to-trial mean','exc. conductance trial-to-trial mean'],st_orientation=0,st_contrast=100)
	    mozaik.analysis.analysis.AnalogSignal_PerNeuronBetweenSignalCorrelation(dsv,ParameterSet({'value_name1':'inh. conductance trial-to-trial mean','value_name2' : 'exc. conductance trial-to-trial mean'})).analyse()

	    dsv = param_filter_query(self.datastore,analysis_algorithm='AnalogSignal_PerNeuronBetweenSignalCorrelation')	    
	    PopulationMeanAndVar(dsv,ParameterSet({})).analyse()
	    PopulationMedian(dsv,ParameterSet({})).analyse()

	    spike_ids = param_filter_query(self.datastore,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
	    spike_ids_inh = param_filter_query(self.datastore,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()

	    wellfit5 = param_filter_query(self.datastore,sheet_name='V1_Exc_L4',st_direct_stimulation_name=None,st_name=['FullfieldDriftingSinusoidalGrating'],st_contrast=5,value_name=['orientation fitting error of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(spike_ids)
	    wellfit100 = param_filter_query(self.datastore,sheet_name='V1_Exc_L4',st_direct_stimulation_name=None,st_name=['FullfieldDriftingSinusoidalGrating'],st_contrast=100,value_name=['orientation fitting error of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(spike_ids)
	    print wellfit5
	    print wellfit100
            wellfit_spike_ids_l4E = numpy.array(spike_ids)[numpy.logical_and(numpy.array(wellfit5) < 0.2,numpy.array(wellfit100) < 0.2)]

	    wellfit5 = param_filter_query(self.datastore,sheet_name='V1_Inh_L4',st_direct_stimulation_name=None,st_name=['FullfieldDriftingSinusoidalGrating'],st_contrast=5,value_name=['orientation fitting error of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(spike_ids_inh)
	    wellfit100 = param_filter_query(self.datastore,sheet_name='V1_Inh_L4',st_direct_stimulation_name=None,st_name=['FullfieldDriftingSinusoidalGrating'],st_contrast=100,value_name=['orientation fitting error of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(spike_ids_inh)
	    print wellfit5
	    print wellfit100
            wellfit_spike_ids_l4I = numpy.array(spike_ids_inh)[numpy.logical_and(numpy.array(wellfit5) < 0.2,numpy.array(wellfit100) < 0.2)]


            hwhh_hc = numpy.array(param_filter_query(self.datastore,sheet_name='V1_Exc_L4',st_direct_stimulation_name=None,st_name=['FullfieldDriftingSinusoidalGrating'],st_contrast=100,value_name=['orientation HWHH of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(wellfit_spike_ids_l4E))
            hwhh_lc = numpy.array(param_filter_query(self.datastore,sheet_name='V1_Exc_L4',st_direct_stimulation_name=None,st_name=['FullfieldDriftingSinusoidalGrating'],st_contrast=5,value_name=['orientation HWHH of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(wellfit_spike_ids_l4E))

            self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.mean(hwhh_hc),period=None,value_name = 'Mean HWHH of responsive neurons',sheet_name='V1_Exc_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        
            self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.mean(abs(hwhh_hc-hwhh_lc)),period=None,value_name = 'Mean HWHH difference',sheet_name='V1_Exc_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        
            self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.median(hwhh_hc),period=None,value_name = 'Median HWHH of responsive neurons',sheet_name='V1_Exc_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        
            self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.median(abs(hwhh_hc-hwhh_lc)),period=None,value_name = 'Median HWHH difference',sheet_name='V1_Exc_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        
	    self.datastore.full_datastore.add_analysis_result(SingleValue(value=(len(spike_ids) - len(wellfit_spike_ids_l4E))/len(spike_ids),period=None,value_name = 'Faild fits percantage',sheet_name='V1_Exc_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        
	    
            hwhh_hc = numpy.array(param_filter_query(self.datastore,sheet_name='V1_Inh_L4',st_direct_stimulation_name=None,st_name=['FullfieldDriftingSinusoidalGrating'],st_contrast=100,value_name=['orientation HWHH of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(wellfit_spike_ids_l4I))
            hwhh_lc = numpy.array(param_filter_query(self.datastore,sheet_name='V1_Inh_L4',st_direct_stimulation_name=None,st_name=['FullfieldDriftingSinusoidalGrating'],st_contrast=5,value_name=['orientation HWHH of Firing rate'],ads_unique=True).get_analysis_result()[0].get_value_by_id(wellfit_spike_ids_l4I))

            self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.mean(hwhh_hc),period=None,value_name = 'Mean HWHH of responsive neurons',sheet_name='V1_Inh_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        
            self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.mean(abs(hwhh_hc-hwhh_lc)),period=None,value_name = 'Mean HWHH difference',sheet_name='V1_Inh_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        
            self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.median(hwhh_hc),period=None,value_name = 'Median HWHH of responsive neurons',sheet_name='V1_Inh_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        
            self.datastore.full_datastore.add_analysis_result(SingleValue(value=numpy.median(abs(hwhh_hc-hwhh_lc)),period=None,value_name = 'Median HWHH difference',sheet_name='V1_Inh_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        
	    self.datastore.full_datastore.add_analysis_result(SingleValue(value=(len(spike_ids_inh) - len(wellfit_spike_ids_l4E))/len(spike_ids_inh),period=None,value_name = 'Faild fits percentage',sheet_name='V1_Inh_L4',tags=self.tags,analysis_algorithm=self.__class__.__name__,stimulus_id=None))        
	    
	    param_filter_query(self.datastore,sheet_name=["V1_Exc_L4","V1_Exc_L4"],identifier=['SingleValue','PerNeuronValue']).remove_ads_outside_of_dsv()





def ana1(data_store):
    print "DSADSA"
    sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4','V1_Inh_L4']))
    exc_sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4']))
    TrialAveragedFiringRate(param_filter_query(data_store,st_name="FullfieldDriftingSinusoidalGrating"),ParameterSet({})).analyse()

    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name='V1_Exc_L4',value_name='Firing rate')    
    GaussianTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()

    dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name='V1_Inh_L4',value_name='Firing rate')
    GaussianTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()

    PSTH(param_filter_query(data_store),ParameterSet({'bin_length' : 10.0 })).analyse()
    NeuronToNeuronAnalogSignalCorrelations(param_filter_query(data_store,analysis_algorithm='PSTH',st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({'convert_nan_to_zero' : True})).analyse()

    dsv = param_filter_query(data_store,analysis_algorithm='NeuronToNeuronAnalogSignalCorrelations')
    PopulationMeanAndVar(dsv,ParameterSet({})).analyse()

    #dsv = queries.param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name = 'V1_Exc_L4',st_orientation=0,st_trial=1)
    Irregularity(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()
    PopulationMeanAndVar(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()


    AAA1(data_store,ParameterSet({})).analyse()

    data_store.save()


def analysis(data_store,analog_ids,analog_ids_inh,gratings,bars):
    sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4','V1_Inh_L4']))
    exc_sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4']))

    TrialAveragedFiringRate(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()

    if bars:
        TrialAveragedFiringRate(param_filter_query(data_store,st_name="FlashedBar"),ParameterSet({})).analyse()

    if gratings:
        TrialAveragedFiringRate(param_filter_query(data_store,st_name="FullfieldDriftingSinusoidalGrating"),ParameterSet({})).analyse()

    Irregularity(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()

    PSTH(param_filter_query(data_store),ParameterSet({'bin_length' : 10.0 })).analyse()

    NeuronToNeuronAnalogSignalCorrelations(param_filter_query(data_store,analysis_algorithm='PSTH',st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({'convert_nan_to_zero' : True})).analyse()
    PopulationMeanAndVar(param_filter_query(data_store,st_direct_stimulation_name="None",st_name='InternalStimulus'),ParameterSet({})).analyse()

    dsv = param_filter_query(data_store,sheet_name=exc_sheets)
    ActionPotentialRemoval(dsv,ParameterSet({'window_length': 5.0})).analyse()

    dsv = param_filter_query(data_store,analysis_algorithm='ActionPotentialRemoval')
    TrialVariability(data_store,ParameterSet({'vm': True,  'cond_exc': True, 'cond_inh': True})).analyse()
    TrialMean(data_store,ParameterSet({'vm': True,  'cond_exc': True, 'cond_inh': True})).analyse()

    dsv = param_filter_query(data_store,st_name='InternalStimulus',st_direct_stimulation_name="None")
    Analog_MeanSTDAndFanoFactor(dsv,ParameterSet({})).analyse()

    if gratings:

        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name=sheets)    
        GaussianTuningCurveFit(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()
        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name=sheets)   
        Analog_F0andF1(dsv,ParameterSet({})).analyse()

        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate',sheet_name=sheets)  
        PeriodicTuningCurvePreferenceAndSelectivity_VectorAverage(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()

        ModulationRatio(param_filter_query(data_store,sheet_name=exc_sheets,st_contrast=[100]),ParameterSet({})).analyse()

        dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm='TrialAveragedFiringRate')    
        CircularVarianceOfTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation'})).analyse()

    data_store.save()



def perform_analysis_and_visualization(data_store,gratings,bars,nat_movies=False):
    
    sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4','V1_Inh_L4']))
    exc_sheets = list(set(data_store.sheets()) & set(['V1_Exc_L4']))
    
    NeuronAnnotationsToPerNeuronValues(data_store,ParameterSet({})).analyse()
    
    analog_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_esyn_ids()
    analog_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_esyn_ids()
    spike_ids = param_filter_query(data_store,sheet_name="V1_Exc_L4").get_segments()[0].get_stored_spike_train_ids()
    spike_ids_inh = param_filter_query(data_store,sheet_name="V1_Inh_L4").get_segments()[0].get_stored_spike_train_ids()
    
    
    l4_exc_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Exc_L4')
    l4_exc_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Exc_L4')
    l4_exc = analog_ids[numpy.argmin([circular_dist(o,numpy.pi/2,numpy.pi)  for (o,p) in zip(l4_exc_or[0].get_value_by_id(analog_ids),l4_exc_phase[0].get_value_by_id(analog_ids))])]
    l4_inh_or = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentOrientation', sheet_name = 'V1_Inh_L4')
    l4_inh_phase = data_store.get_analysis_result(identifier='PerNeuronValue',value_name = 'LGNAfferentPhase', sheet_name = 'V1_Inh_L4')
    l4_inh = analog_ids_inh[numpy.argmin([circular_dist(o,numpy.pi/2,numpy.pi)  for (o,p) in zip(l4_inh_or[0].get_value_by_id(analog_ids_inh),l4_inh_phase[0].get_value_by_id(analog_ids_inh))])]
    l4_exc_or_many = numpy.array(l4_exc_or[0].ids)[numpy.nonzero(numpy.array([circular_dist(o,0,numpy.pi)  for (o,p) in zip(l4_exc_or[0].values,l4_exc_phase[0].values)]) < 0.1)[0]]
    
    l4_exc_or_many = list(set(l4_exc_or_many) &  set(spike_ids))
        
    orr = list(set([MozaikParametrized.idd(s).orientation for s in queries.param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',st_contrast=100).get_stimuli()]))                
        
    l4_exc_or_many_analog = numpy.array(analog_ids)[numpy.nonzero(numpy.array([circular_dist(l4_exc_or[0].get_value_by_id(i),0,numpy.pi)  for i in analog_ids]) < 0.1)[0]]
    l4_inh_or_many_analog = numpy.array(analog_ids_inh)[numpy.nonzero(numpy.array([circular_dist(l4_inh_or[0].get_value_by_id(i),0,numpy.pi)  for i in analog_ids_inh]) < 0.15)[0]]
    
    lgn_on_ids = param_filter_query(data_store,sheet_name="X_ON").get_segments()[0].get_stored_spike_train_ids()
    lgn_off_ids = param_filter_query(data_store,sheet_name="X_OFF").get_segments()[0].get_stored_spike_train_ids()    

    #analysis(data_store,analog_ids,analog_ids_inh,gratings,bars)

    if bars:
        dsv = queries.param_filter_query(data_store,st_name='FlashedBar')
        for ads in dsv.get_analysis_result():
            sid = MozaikParametrized.idd(ads.stimulus_id)
            sid.x=0
            ads.stimulus_id = str(sid)
        for seg in dsv.get_segments():    
            sid = MozaikParametrized.idd(seg.annotations['stimulus'])
            sid.x=0
            seg.annotations['stimulus'] = str(sid)

    def overviews(dsv,name_prefix):

        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L4ExcAnalog1.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L4InhAnalog1.png').plot()    
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L4ExcAnalog2.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L4InhAnalog2.png').plot()    
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L4ExcAnalog3.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'L4InhAnalog3.png').plot()    

        RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'SSExcRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})
        RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'SSInhRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})
    
        RasterPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neurons' : lgn_on_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'XONRaster.png').plot({'SpikeRasterPlot.group_trials':True})
        RasterPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neurons' : lgn_off_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name=name_prefix+'XOFFRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})



    # self sustained plotting
    dsv = param_filter_query(data_store,st_name='InternalStimulus')   
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalog1.png').plot()
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[0], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhAnalog1.png').plot()    
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalog2.png').plot()
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[1], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhAnalog2.png').plot()    
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcAnalog3.png').plot()
    OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : analog_ids_inh[2], 'sheet_activity' : {}, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhAnalog3.png').plot()    
    
    RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : spike_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSExcRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})
    RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neurons' : spike_ids_inh,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='SSInhRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})

    if bars:
        dsv = queries.param_filter_query(data_store,st_name='FlashedBar',sheet_name = 'V1_Exc_L4',st_relative_luminance=1)
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='BarExcAnalog1.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[1], 'sheet_activity' : {}, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='BarExcAnalog2.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[2], 'sheet_activity' : {}, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='BarExcAnalog3.png').plot()


        dsv = queries.param_filter_query(data_store,st_name='FlashedBar',y_axis_name=['Vm (no AP) trial-to-trial mean','exc. conductance trial-to-trial mean','inh. conductance trial-to-trial mean'],sheet_name = 'V1_Exc_L4',analysis_algorithm='TrialMean',st_relative_luminance=1)
        PlotTemporalTuningCurve(dsv, ParameterSet({'parameter_name' : 'y', 'neurons': list(analog_ids[:4]), 'sheet_name' : 'V1_Exc_L4','mean' : False}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='AnalogSignalEvolutionOnBar1.png').plot({'*.x_lim'  :(0,200),'*.y_label'  : 'offset'})
        dsv = queries.param_filter_query(data_store,st_name='FlashedBar',y_axis_name=['Vm (no AP) trial-to-trial mean','exc. conductance trial-to-trial mean','inh. conductance trial-to-trial mean'],sheet_name = 'V1_Exc_L4',analysis_algorithm='TrialMean',st_relative_luminance=0)
        PlotTemporalTuningCurve(dsv, ParameterSet({'parameter_name' : 'y', 'neurons': list(analog_ids[:4]), 'sheet_name' : 'V1_Exc_L4','mean' : False}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='AnalogSignalEvolutionOffBar1.png').plot({'*.x_lim'  :(0,200),'*.y_label'  : 'offset'})

        dsv = queries.param_filter_query(data_store,st_name='FlashedBar',y_axis_name=['Vm (no AP) trial-to-trial mean','exc. conductance trial-to-trial mean','inh. conductance trial-to-trial mean'],sheet_name = 'V1_Exc_L4',analysis_algorithm='TrialMean',st_relative_luminance=1)
        PlotTemporalTuningCurve(dsv, ParameterSet({'parameter_name' : 'y', 'neurons': list(analog_ids[4:8]), 'sheet_name' : 'V1_Exc_L4','mean' : False}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='AnalogSignalEvolutionOnBar2.png').plot({'*.x_lim'  :(0,200),'*.y_label'  : 'offset'})
        dsv = queries.param_filter_query(data_store,st_name='FlashedBar',y_axis_name=['Vm (no AP) trial-to-trial mean','exc. conductance trial-to-trial mean','inh. conductance trial-to-trial mean'],sheet_name = 'V1_Exc_L4',analysis_algorithm='TrialMean',st_relative_luminance=0)
        PlotTemporalTuningCurve(dsv, ParameterSet({'parameter_name' : 'y', 'neurons': list(analog_ids[4:8]), 'sheet_name' : 'V1_Exc_L4','mean' : False}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='AnalogSignalEvolutionOffBar2.png').plot({'*.x_lim'  :(0,200),'*.y_label'  : 'offset'})

        #dsv = queries.param_filter_query(data_store,st_name='FlashedBar',y_axis_name=['Vm (no AP) trial-to-trial mean','exc. conductance trial-to-trial mean','inh. conductance trial-to-trial mean'],sheet_name = 'V1_Exc_L4',analysis_algorithm='TrialMean',st_relative_luminance=1)
        #PlotTemporalTuningCurve(dsv, ParameterSet({'parameter_name' : 'y', 'neurons': list(analog_ids[8:12]), 'sheet_name' : 'V1_Exc_L4','mean' : False}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='AnalogSignalEvolutionOnBar3.png').plot({'*.x_lim'  :(0,200),'*.y_label'  : 'offset'})
        #dsv = queries.param_filter_query(data_store,st_name='FlashedBar',y_axis_name=['Vm (no AP) trial-to-trial mean','exc. conductance trial-to-trial mean','inh. conductance trial-to-trial mean'],sheet_name = 'V1_Exc_L4',analysis_algorithm='TrialMean',st_relative_luminance=0)
        #PlotTemporalTuningCurve(dsv, ParameterSet({'parameter_name' : 'y', 'neurons': list(analog_ids[8:12]), 'sheet_name' : 'V1_Exc_L4','mean' : False}),fig_param={'dpi' : 100,'figsize': (20,10)},plot_file_name='AnalogSignalEvolutionOffBar3.png').plot({'*.x_lim'  :(0,200),'*.y_label'  : 'offset'})
	
	if False:
            dsv = param_filter_query(data_store,st_name='FlashedBar',analysis_algorithm='TrialMean',st_relative_luminance=1.0,y_axis_name=['exc. conductance trial-to-trial mean'],sheet_name='V1_Exc_L4',st_y=-0.5769230769230769)
	    mozaik.visualization.plotting.AnalogSignalListPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4','neurons' : analog_ids.tolist(), 'mean' : True}),fig_param={'dpi' : 100,'figsize': (10,4)},plot_file_name='BarMeanExcCondON.png').plot()
	    dv = param_filter_query(data_store,st_name='FlashedBar',analysis_algorithm='TrialMean',st_relative_luminance=1.0,y_axis_name=['inh. conductance trial-to-trial mean'],sheet_name='V1_Exc_L4',st_y=-0.5769230769230769)
	    mozaik.visualization.plotting.AnalogSignalListPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4','neurons' : analog_ids.tolist(), 'mean' : True}),fig_param={'dpi' : 100,'figsize': (10,4)},plot_file_name='BarMeanInhCondON.png').plot()
            dsv = param_filter_query(data_store,st_name='FlashedBar',analysis_algorithm='TrialMean',st_relative_luminance=1.0,y_axis_name=['vm trial-to-trial mean'],sheet_name='V1_Exc_L4',st_y=-0.5769230769230769)
            mozaik.visualization.plotting.AnalogSignalListPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4','neurons' : analog_ids.tolist(), 'mean' : True}),fig_param={'dpi' : 100,'figsize': (10,4)},plot_file_name='BarMeanVmON.png').plot()

            dsv = param_filter_query(data_store,st_name='FlashedBar',analysis_algorithm='TrialMean',st_relative_luminance=0.0,y_axis_name=['exc. conductance trial-to-trial mean'],sheet_name='V1_Exc_L4',st_y=-0.5769230769230769)
            mozaik.visualization.plotting.AnalogSignalListPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4','neurons' : analog_ids.tolist(), 'mean' : True}),fig_param={'dpi' : 100,'figsize': (10,4)},plot_file_name='BarMeanExcCondOFF.png').plot()
            dsv = param_filter_query(data_store,st_name='FlashedBar',analysis_algorithm='TrialMean',st_relative_luminance=0.0,y_axis_name=['inh. conductance trial-to-trial mean'],sheet_name='V1_Exc_L4',st_y=-0.5769230769230769)
            mozaik.visualization.plotting.AnalogSignalListPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4','neurons' : analog_ids.tolist(), 'mean' : True}),fig_param={'dpi' : 100,'figsize': (10,4)},plot_file_name='BarMeanInhCondOFF.png').plot()
            dsv = param_filter_query(data_store,st_name='FlashedBar',analysis_algorithm='TrialMean',st_relative_luminance=0.0,y_axis_name=['vm trial-to-trial mean'],sheet_name='V1_Exc_L4',st_y=-0.5769230769230769)
            mozaik.visualization.plotting.AnalogSignalListPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4','neurons' : analog_ids.tolist(), 'mean' : True}),fig_param={'dpi' : 100,'figsize': (10,4)},plot_file_name='BarMeanVmOFF.png').plot()

        dsv = queries.param_filter_query(data_store,st_name='FlashedBar',st_relative_luminance=1.0)
        RasterPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neurons' : lgn_on_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Bar_XONRaster.png').plot({'SpikeRasterPlot.group_trials':True})
        RasterPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neurons' : lgn_off_ids,'trial_averaged_histogram': False, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='Bar_XOFFRasterL4.png').plot({'SpikeRasterPlot.group_trials':True})
	
	#dsv = param_filter_query(data_store,st_name=['FlashedBar'],st_relative_luminance=0,st_y=-0.5769230769230769)    
	#RasterPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neurons' : lgn_off_ids, 'trial_averaged_histogram' : True, 'spontaneous' : False}),plot_file_name="Bar_LGN_OFF_Bar_OFF.png",fig_param={'dpi' : 100,'figsize': (14,4)}).plot()
	#RasterPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neurons' : lgn_on_ids, 'trial_averaged_histogram' : True, 'spontaneous' : False}),plot_file_name="Bar_LGN_ON_Bar_OFF.png",fig_param={'dpi' : 100,'figsize': (14,4)}).plot()

	#dsv = param_filter_query(data_store,st_name=['FlashedBar'],st_relative_luminance=1,st_y=-0.5769230769230769)    
	#RasterPlot(dsv,ParameterSet({'sheet_name' : 'X_OFF', 'neurons' : lgn_off_ids, 'trial_averaged_histogram' : True, 'spontaneous' : False}),plot_file_name="Bar_LGN_OFF_Bar_ON.png",fig_param={'dpi' : 100,'figsize': (14,4)}).plot()
	#RasterPlot(dsv,ParameterSet({'sheet_name' : 'X_ON', 'neurons' : lgn_on_ids, 'trial_averaged_histogram' : True, 'spontaneous' : False}),plot_file_name="Bar_LGN_ON_Bar_ON.png",fig_param={'dpi' : 100,'figsize': (14,4)}).plot()

    


    if nat_movies:        
	dsv = param_filter_query(data_store,st_name='NaturalImageWithEyeMovement')   
        overviews(dsv,"NI_")

        TrialToTrialVariabilityComparison(data_store,ParameterSet({'sheet_name1' : 'V1_Exc_L4','sheet_name2' : 'V1_Exc_L4','data_dg' : 0.93 , 'data_ni' : 1.19}),fig_param={'dpi' : 200,'figsize': (15,7.5)},plot_file_name='TrialToTrialVariabilityComparison.png').plot()

        dsv = param_filter_query(data_store,st_name='NaturalImageWithEyeMovement')    
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : l4_exc, 'sheet_activity' : {}, 'spontaneous' : True}),plot_file_name='NMExc.png',fig_param={'dpi' : 100,'figsize': (28,12)}).plot({'Vm_plot.y_lim' : (-70,-50),'Conductance_plot.y_lim' : (0,50.0)})
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Inh_L4', 'neuron' : l4_inh, 'sheet_activity' : {}, 'spontaneous' : True}),plot_file_name='NMInh.png',fig_param={'dpi' : 100,'figsize': (28,12)}).plot({'Vm_plot.y_lim' : (-70,-50),'Conductance_plot.y_lim' : (0,50.0)})
                    
        #TrialCrossCorrelationAnalysis(data_store,ParameterSet({'neurons1' : list(analog_ids),'sheet_name1' : 'V1_Exc_L4','neurons2' : list([]),'sheet_name2' : 'V1_Exc_L2/3', 'window_length' : 250}),fig_param={"dpi" : 100,"figsize": (15,6.5)},plot_file_name="trial-to-trial-cross-correlation.png").plot({'*.Vm.title' : None, '*.fontsize' : 19})


    if gratings:    
	dsv = queries.param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name = 'V1_Exc_L4',st_orientation=0,st_trial=1)
        RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : list(spike_ids),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedExcRasterTrial1.png').plot({'SpikeRasterPlot.group_trials':True})
	dsv = queries.param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name = 'V1_Exc_L4',st_orientation=0,st_trial=2)
        RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : list(spike_ids),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedExcRasterTrial2.png').plot({'SpikeRasterPlot.group_trials':True})
	dsv = queries.param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name = 'V1_Exc_L4',st_orientation=0,st_trial=3)
        RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : list(spike_ids),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedExcRasterTrial3.png').plot({'SpikeRasterPlot.group_trials':True})
	dsv = queries.param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name = 'V1_Exc_L4',st_orientation=0,st_trial=4)
        RasterPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neurons' : list(spike_ids),'trial_averaged_histogram': False, 'spontaneous' : True}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='EvokedExcRasterTrial4.png').plot({'SpikeRasterPlot.group_trials':True})

        dsv = queries.param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',sheet_name = 'V1_Exc_L4',st_orientation=0)
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[0], 'sheet_activity' : {}, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='GratingExcAnalog1.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[1], 'sheet_activity' : {}, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='GratingBarExcAnalog2.png').plot()
        OverviewPlot(dsv,ParameterSet({'sheet_name' : 'V1_Exc_L4', 'neuron' : analog_ids[2], 'sheet_activity' : {}, 'spontaneous' : False}),fig_param={'dpi' : 100,'figsize': (28,12)},plot_file_name='GratingBarExcAnalog3.png').plot()

        # tuninc curves
        #dsv = param_filter_query(data_store,st_name='FullfieldDriftingSinusoidalGrating',analysis_algorithm=['TrialAveragedFiringRate','Analog_F0andF1'])    
        #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids[:6]), 'sheet_name' : 'V1_Exc_L4','centered'  : True,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='TCsExcL4.png',fig_param={'dpi' : 100,'figsize': (15,7.5)}).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})
        #PlotTuningCurve(dsv,ParameterSet({'parameter_name' : 'orientation', 'neurons': list(analog_ids_inh[:6]), 'sheet_name' : 'V1_Inh_L4','centered'  : True,'mean' : False, 'polar' : False, 'pool'  : False}),plot_file_name='TCsInhL4.png',fig_param={'dpi' : 100,'figsize': (15,7.5)}).plot({'TuningCurve F0_Inh_Cond.y_lim' : (0,180) , 'TuningCurve F0_Exc_Cond.y_lim' : (0,80)})


        Kremkow_plots.OrientationTuningSummary(data_store,ParameterSet({'exc_sheet_name': 'V1_Exc_L4','inh_sheet_name': 'V1_Inh_L4'}),fig_param={'dpi' : 100,'figsize': (15,9)},plot_file_name='OrientationTuningSummaryL4.png').plot()            

        OrientationTuningSummaryAnalogSignals(data_store,ParameterSet({'exc_sheet_name1': 'V1_Exc_L4','inh_sheet_name1': 'V1_Inh_L4','exc_sheet_name2': 'None','inh_sheet_name2': 'None'}),fig_param={'dpi' : 200,'figsize': (18,12)},plot_file_name='OrientationTuningSummaryAnalogSignals.png').plot({'*.fontsize' : 19,'*.y_lim' : (0,None)})            
        
	
        MRfig(param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'None'}),plot_file_name='MR.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()
        #MRfigReal(param_filter_query(data_store,sheet_name=['V1_Exc_L4'],st_contrast=[100],st_name='FullfieldDriftingSinusoidalGrating'),ParameterSet({'SimpleSheetName' : 'V1_Exc_L4','ComplexSheetName' : 'None'}),plot_file_name='MR.png',fig_param={'dpi' : 100,'figsize': (19,12)}).plot()
        
	TrialToTrialVariabilityComparison(data_store,ParameterSet({'sheet_name1' : 'V1_Exc_L4','sheet_name2' : 'V1_Exc_L4','data_dg' : 0.93 , 'data_ni' : 1.19}),fig_param={'dpi' : 200,'figsize': (15,7.5)},plot_file_name='TrialToTrialVariabilityComparison.png').plot()
        
        TrialCrossCorrelationAnalysis(data_store,ParameterSet({'neurons1' : list(analog_ids),'sheet_name1' : 'V1_Exc_L4','neurons2' : [],'sheet_name2' : 'V1_Exc_L2/3', 'window_length' : 250}),fig_param={"dpi" : 100,"figsize": (15,6.5)},plot_file_name="trial-to-trial-cross-correlation.png").plot({'*.Vm.title' : None, '*.fontsize' : 19})

    SpontActOverview(data_store,ParameterSet({'l4_exc_neuron' : analog_ids[0], 'l4_inh_neuron' : analog_ids_inh[0],'l23_exc_neuron' : -1,'l23_inh_neuron' : -1}),plot_file_name='SpontActOverview.png', fig_param={'dpi' : 200,'figsize': (18,14.5)}).plot()
    SpontStatisticsOverview(data_store,ParameterSet({}), fig_param={'dpi' : 200,'figsize': (18,12)},plot_file_name='SpontStatisticsOverview.png').plot()

