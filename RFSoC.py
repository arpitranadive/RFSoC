# rfSoC driver with inbuilt pulse-generation
# Recommended practice is to add pulses to parameter snaps.
#                                                    -- Arpit


# Updated version to work with QCodes 0.31.0, coded on 28 Jan. 2022
# Long DC pulses were also implemented but can show some data corruption bugs in extreme scenarios.





import time
import datetime
import numpy as np
import sys
import struct
import ctypes  # only for DLL-based instrument
import pickle as pk

import qcodes as qc
from qcodes import (Instrument, VisaInstrument,
					ManualParameter, MultiParameter,
					validators as vals)
from qcodes.instrument.channel import InstrumentChannel
from qcodes.instrument.parameter import ParameterWithSetpoints, Parameter

#import SequenceGeneration_v2 as sqg
from qcodes.utils.delaykeyboardinterrupt import DelayedKeyboardInterrupt
from qcodes.utils.validators import Numbers, Arrays

from pprint import pprint
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd
from IPython.display import display, HTML
from ipywidgets import IntProgress
import matplotlib.pyplot as plt

import functools
import operator
from itertools import chain

sys.path.append(r"C:\Users\nicolas.roch\measurement\Experiment\Drivers")
#from progress_barV2 import bar

import logging
log = logging.getLogger(__name__)





class GeneratedSetPoints(Parameter):
	"""
	A parameter that generates a setpoint array from start, stop and num points
	parameters.
	"""
	def __init__(self, startparam, stopparam, numpointsparam, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self._startparam = startparam
		self._stopparam = stopparam
		self._numpointsparam = numpointsparam

	def get_raw(self):
		return np.linspace(self._startparam(), self._stopparam() -1,
							  self._numpointsparam())


class RAW(Parameter):

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self._channel = self._instrument._adc_channel

	def get_raw(self):
		time.sleep(0.2)

		dataI, dataQ = self._instrument._parent.get_single_readout_pulse()

		# print(self._channel)
		if self._channel in np.arange(1,9):
			data_ret = dataI[self._channel-1]
		else:
			log.warning('Wrong parameter.')

		return data_ret


class RAW_ALL(Parameter):

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)

	def get_raw(self):
		time.sleep(0.2)

		dataI, dataQ = self._instrument.get_readout_pulse()

		data_ret = dataI

		return data_ret


class IQINT(Parameter):

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self._channel = self._instrument._adc_channel

	def get_raw(self):
		time.sleep(0.2)

		dataI, dataQ = self._instrument._parent.get_single_readout_pulse()

		# print(self._channel)
		if self._channel in np.arange(1,9):
			data_retI = dataI[self._channel-1]
			data_retQ = dataQ[self._channel-1]
		else:
			log.warning('Wrong parameter.')

		return data_retI, data_retQ


class IQINT_ALL(Parameter):

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)


	def get_raw(self):
		time.sleep(0.2)

		data_retI, data_retQ = self._instrument.get_readout_pulse()

		return data_retI, data_retQ


class IQINT_ALL_read_header(Parameter):

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)


	def get_raw(self):
		time.sleep(0.2)

		data_retI, data_retQ = self._instrument.get_readout_pulse()

		return data_retI, data_retQ


class IQINT_AVG(Parameter):

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)

	def get_raw(self):

		time.sleep(0.2)

		data_retI, data_retQ = self._instrument.get_readout_pulse()#.get_readout_pulse() ###### Martina 08/11/2020

		Sq_I = [[],[],[],[],[],[],[],[]]
		Sq_Q = [[],[],[],[],[],[],[],[]]

		Sq_I_list = [[],[],[],[],[],[],[],[]]
		Sq_Q_list = [[],[],[],[],[],[],[],[]]

		for i in range(8):

			if len(data_retI[i])>0:

				for j in range(len(data_retI[i])):

					if len(data_retI[i][j])>0:

						Sq_I[i] = np.append(Sq_I[i],np.mean(data_retI[i][j]))
						Sq_Q[i] = np.append(Sq_Q[i],np.mean(data_retQ[i][j]))
						Sq_I_list[i].append(np.mean(data_retI[i][j]))
						Sq_Q_list[i].append(np.mean(data_retQ[i][j]))

		for i in range(8):

			Sq_I[i] = np.array(Sq_I_list[i])
			Sq_Q[i] = np.array(Sq_Q_list[i])

		return Sq_I,Sq_Q


class ADC_power(Parameter):

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)

	def get_raw(self):

		time.sleep(0.2)

		data_retI, data_retQ = self._instrument.get_readout_pulse()#.get_readout_pulse() ###### Martina 08/11/2020

		Sq_I_list = [[],[],[],[],[],[],[],[]]
		Sq_Q_list = [[],[],[],[],[],[],[],[]]
		Pow = [[],[],[],[],[],[],[],[]]

		for i in range(8):

			if len(data_retI[i])>0:

				for j in range(len(data_retI[i])):

					if len(data_retI[i][j])>0:

						Sq_I_list[i].append(np.mean(data_retI[i][j]**2))
						Sq_Q_list[i].append(np.mean(data_retQ[i][j]**2))

		for i in range(8):

			Pow[i] = (np.array(Sq_I_list[i]) + np.array(Sq_Q_list[i]))/(50*2)    # RMS power

		return Pow


class ADC_power_dBm(Parameter):

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)

	def get_raw(self):

		time.sleep(0.2)

		data_retI, data_retQ = self._instrument.get_readout_pulse()#.get_readout_pulse() ###### Martina 08/11/2020

		Sq_I_list = [[],[],[],[],[],[],[],[]]
		Sq_Q_list = [[],[],[],[],[],[],[],[]]
		Pow = [[],[],[],[],[],[],[],[]]

		for i in range(8):

			if len(data_retI[i])>0:

				for j in range(len(data_retI[i])):

					if len(data_retI[i][j])>0:

						Sq_I_list[i].append(np.mean(data_retI[i][j]**2))
						Sq_Q_list[i].append(np.mean(data_retQ[i][j]**2))

		for i in range(8):

			Pow[i] = 10*np.log10(1e3*(np.array(Sq_I_list[i]) + np.array(Sq_Q_list[i]))/(50*2))        # RMS power

		return Pow


#TODO : add the different results we would return
class AcqChannel(InstrumentChannel):

	#initializing the array storing channel instances
	objs = []

	def __init__(self, parent: 'RFSoC', name: str, channel: int):

		"""
		Args:
			parent: Instrument that this channel is bound to.
			name: Name to use for this channel.
			channel: channel on the card to use
		"""
		if channel not in np.arange(1,9):
			raise ValueError('channel number must be in between 1 and 8')

		self._adc_channel = channel

		super().__init__(parent, name)

		AcqChannel.objs.append(self)

		self.add_parameter(name='status',
						   label = 'ADC{} status'.format(self._adc_channel),
						   set_cmd='ADC:ADC{} {}'.format(self._adc_channel,'{:d}'),
						   val_mapping={'ON': 1, 'OFF': 0},
						   initial_value='OFF',
						   snapshot_value = False
						   )

		#TODO : add allowed values of decimation and mixer frea
		self.add_parameter(name='decfact',
						   label='ADC{} decimation factor'.format(self._adc_channel),
						   #the decimation works by tiles of two adcs
						   set_cmd='ADC:TILE{}:DECFACTOR {}'.format((self._adc_channel-1)//2,'{:d}'),
						   snapshot_value = False
						   )

		self.add_parameter(name='fmixer',
						   label = 'ADC{} mixer frequency'.format(self._adc_channel),
						   set_cmd='ADC:ADC{}:MIXER {}'.format(self._adc_channel,'{:.4f}'),
						   # get_parser = self.MHz_to_Hz,
						   # set_parser = self.Hz_to_MHz,
						   snapshot_value = False
						   )

		self.status('OFF')

	def MHz_to_Hz(self,value):
		return value*1e6

	def Hz_to_MHz(self,value):
		return value*1e-6

class RFSoC(VisaInstrument):

	# all instrument constructors should accept **kwargs and pass them on to
	# super().__init__
	def __init__(self, name, address, **kwargs):
		# supplying the terminator means you don't need to remove it from every
		# response
		# super().__init__(name, address, terminator='\r\n', **kwargs)
		super().__init__(name, address, terminator='\n', **kwargs)


		self.dummy_array_size_8 = 8
		self.sampling_rate = 2e9
		self.FPGA_clock = 250e6

		self.DAC_amplitude_calib = [1, 1, 1, 1, 1, 1, 1, 1]
		self.ADC_amplitude_calib = [1, 1, 1, 1, 1, 1, 1, 1]

		self.pulses = pd.DataFrame()
		self.seq_looping = pd.DataFrame()
		self.ADC_ch_active = np.zeros(8)
		self.length_vec = [[],[],[],[],[],[],[],[]]
		self.ADC_events = np.array([])
		self.global_sequence_str = ''
		self.global_sequence = []
		self.global_sequence_info = []

		self.display_sequence = True
		self.display_IQ_progress = True
		self.debug_mode = False
		self.debug_bypass_count_check = False
		self.debug_mode_plot_waveforms = False
		self.debug_mode_waveform_string = False
		self.loop_time = False

		self.global_loop = (False,0) 
		self.n_points_total = 0

		self._mux_config_matrix = np.array([[0,1,2,2],
											[1,0,3,3],
										['x','x',0,1],
										['x','x',1,0]])

		self.raw_dump_location = "C:/Data_tmp"


		#Add the channels to the instrument
		for adc_num in np.arange(1,9):
			adc_name='ADC{}'.format(adc_num)
			adc=AcqChannel(self,adc_name,adc_num)
			self.add_submodule(adc_name, adc)


		self.add_parameter('sequence_str',
							get_parser=str,
							initial_value='',
							parameter_class=ManualParameter)

		self.add_parameter('n_points',
							get_parser=int,
							initial_value = int(1),
							parameter_class=ManualParameter)

		self.add_parameter('acquisition_mode',
							label='ADCs acquisition mode',
							get_parser=str,
							vals = vals.Enum('RAW','INT'),
							parameter_class=ManualParameter
							)

		self.add_parameter( name = 'output_format',
							#Format(string) : 'BIN' or 'ASCII'
							label='Output format',
							vals = vals.Enum('ASCII','BIN'),
							set_cmd='OUTPUT:FORMAT ' + '{}',
							initial_value='BIN',
							# get_cmd='OUTPUT:FORMAT?',
							#snapshot_get  = False,
							get_parser=str,
							snapshot_value = False)

		self.add_parameter(name='RAW_ALL',
						   unit='V',
						   label='Raw adc for all channel',
						   parameter_class=RAW_ALL,
						   snapshot_value = False)

		self.add_parameter(name='IQINT_ALL',
						   unit='V',
						   label='Integrated I Q for all channels',
						   parameter_class=IQINT_ALL,
						   snapshot_value = False)

		self.add_parameter(name='IQINT_ALL_read_header',
						   unit='V',
						   label='Integrated I Q for all channels with header check',
						   parameter_class=IQINT_ALL_read_header,
						   snapshot_value = False)

		self.add_parameter(name='IQINT_two_mode_squeezing',
						   unit='V',
						   label='Integrated I Q for 2 channels with header check',
						   parameter_class = ManualParameter,
						   snapshot_value = False)

		self.add_parameter(name='IQINT_AVG',
						   unit='V',
						   label='Integrated averaged I Q for all channels',
						   parameter_class=IQINT_AVG,
						   vals=Arrays(shape=(2, 8, 2)), ### ME
						   snapshot_value = False)

		self.add_parameter('channel_axis',
							unit = 'channel index',
							label = 'Channel index axis for ADC power mode',
							parameter_class = GeneratedSetPoints,
							startparam = self.int_0,
							stopparam = self.int_8,
							numpointsparam = self.int_8,
							snapshot_value = False,
							vals=Arrays(shape=(self.dummy_array_size_8,))
							)

		self.add_parameter(name='ADC_power',
						   unit='W',
						   label='Array of incident power on ADC channels',
						   parameter_class=ADC_power,
						   vals=Arrays(shape=(self.dummy_array_size_8,)),
						   # vals=Arrays(shape=(self.dummy_array_size_8,)),
						   snapshot_value = False)

		self.add_parameter(name='ADC_power_dBm',
						   unit='dBm',
						   label='Array of incident power on ADC channels',
						   parameter_class=ADC_power_dBm,
						   vals=Arrays(shape=(self.dummy_array_size_8,)),
						   snapshot_value = False)

		# legacy used for process_sequencing
		self.add_parameter(name='freq_sync',
						   unit='Hz',
						   initial_value=10e6,
						   label='Reference frequency for synchronizing pulse repetition',
						   get_parser=float,
						   parameter_class=ManualParameter,
						   snapshot_value = False)

		self.add_parameter(name='time_vec',
						   unit='s',
						   label='Time vector used for time sequence measurement'
						   )



	def pulses_sequence(self):

		"""
		This function takes the pulses sequence stored in self.pulses and convert it in a Panda Dataframe 
		that can be used by 'LUT_and_address_filling' and 
		process_sequencing_IQ_table() (or process_sequencing() in the legacy mode)
		If self.display_sequence is True it also provide à graphical representation of the pulse sequence. 
		Note that the displayed sequence just represent the first instance of the loop set by self.n_points.
		Also only the first iteration of the loop set by a n_rep>1 is shown (in dark color).
		"""

		self.global_loop = (False,0)

		log.info('Started sequence processing'+'  \n')

		pulses_raw_df = self.pulses
		seq_looping = self.seq_looping

		pulses_raw_df['index'] = range(1, len(pulses_raw_df) + 1)
		pulses_raw_df = pulses_raw_df.set_index('index')

		# --- Check pulse labeling 

		if len(set(pulses_raw_df['label'])) < len(pulses_raw_df['label']):

			log.error('Duplicate Labels: Labels need to be unique for consistent identification of pulse hierarchy.')

		try:

			if len(set(seq_looping['label'])) < len(seq_looping['label']):

				log.error('Duplicate Labels: Labels need to be unique for consistent identification of looping hierarchy.')

		except:

			log.info('No looping programmed.')

		# --- check pulse timing compatibility with RFSoC timing constrains 

		for index, row in pulses_raw_df.iterrows():

			if row['length'] == self.time_conversion(row['length']) and row['start'] == self.time_conversion(row['start']):

				pass

			else:

				pulses_raw_df.at[index,'length'] = self.time_conversion(row['length'])
				pulses_raw_df.at[index,'start'] = self.time_conversion(row['start'])
				
				log.error('Incompatible timing: RFSoC supports events only in factors of 4 ns.')

		# --- fixed hierarchy resolution, please let me know if issues are found -Arpit (202301)

		pulses_raw_df.set_index('label', inplace = True)

		resolve_hierarchy = True
		while resolve_hierarchy:

			for index, row in pulses_raw_df.iterrows():

				if row['parent'] != None:

					if pulses_raw_df.loc[row['parent']]['parent'] == None and pulses_raw_df.loc[row['parent']]['partner'] == None:

						pulses_raw_df.loc[index,'start'] = pulses_raw_df.loc[index,'start'] +  pulses_raw_df.loc[row['parent']]['start'] + (pulses_raw_df.loc[row['parent']]['length'])
						pulses_raw_df.loc[index,'parent'] = None

				if row['partner'] != None:

					if pulses_raw_df.loc[row['partner']]['parent'] == None and pulses_raw_df.loc[row['partner']]['partner'] == None  and row['module'] == 'DAC':

						pulses_raw_df.loc[index,'start'] = pulses_raw_df.loc[index,'start'] +  pulses_raw_df.loc[row['partner']]['start']
						pulses_raw_df.loc[index,'param']['partnet_length'] = pulses_raw_df.loc[row['partner']]['length'] 
						pulses_raw_df.loc[index,'param']['partnet_start'] = pulses_raw_df.loc[row['partner']]['start'] 
						pulses_raw_df.loc[index,'param']['partnet_label'] = row['partner']
						pulses_raw_df.loc[index,'partner'] = None
						
					elif pulses_raw_df.loc[row['partner']]['parent'] == None and pulses_raw_df.loc[row['partner']]['partner'] == None  and row['module'] == 'ADC':

						pulses_raw_df.loc[index,'start'] = pulses_raw_df.loc[index,'start'] +  pulses_raw_df.loc[row['partner']]['start']
						pulses_raw_df.loc[index,'param'][0]['partnet_length'] = pulses_raw_df.loc[row['partner']]['length'] 
						pulses_raw_df.loc[index,'param'][0]['partnet_start'] = pulses_raw_df.loc[row['partner']]['start'] 
						pulses_raw_df.loc[index,'param'][0]['partnet_label'] = row['partner']
						pulses_raw_df.loc[index,'partner'] = None

			resolve_hierarchy = False
			for val in pulses_raw_df['parent']:
				if val != None:
					resolve_hierarchy = True
			for val in pulses_raw_df['partner']:
				if val != None:
					resolve_hierarchy = True

		if self.debug_mode:

			print('Hierarchy resolution...')
			display(pulses_raw_df)

		# --- Initialization 

		pulses_df = pd.DataFrame()
		loops_df = pd.DataFrame()

		time_ADC = [0,0,0,0,0,0,0,0]
		time_DAC = [0,0,0,0,0,0,0,0]
		length_vec = [[],[],[],[],[],[],[],[]]
		# color of the displayed pulses with n_rep=1
		color_dict = {}
		wait_color = int("D3D3D3", 16)
		DAC_color = int("c1666b", 16)		#306cc7
		ADC_color = int("4281a4", 16)		#db500b
		loop_color = int("48a9a6", 16)
		# # color of the displayed pulses with n_rep>1
		# wait_color_rep = int("808080", 16)
		# DAC_color_rep = int("234e90", 16)
		# ADC_color_rep = int("913608", 16)

		wait_count = 0 # counter used to label the wait 
		termination_time = 0 # keep track of the latest event
		ch_demod = None # either or not the ADC LUT will be used during the sequence

		self.ADC_events = np.array([])


		# --- Beginning of the treatment (loops)

		for index, row in seq_looping.iterrows():

			if self.debug_mode:

				print('\n------------- Loop processing -------------\n')

			if row['mode'] == 'global':

				self.global_loop = (True,int(row['rep']))

				label = 'looping ' + str(int(row['rep'])) + ' times'
				channel = 'loop global'
				start = 0
				stop = -1
				time = 0
				mode = row['mode']
				color = color = '#{0:06X}'.format(loop_color)

				loops_df = pd.concat([loops_df, pd.DataFrame.from_records([dict(label=label, start=start, stop=stop, time=time, mode=mode, Channel=channel, color=str(color))])], ignore_index=True)

				if self.debug_mode:

					print('global loop detected for n = '+str(self.global_loop[1]))


			if self.debug_mode:

				display(loops_df)
				print('\n------------- Loop processing finished -------------\n')

		# --- Beginning of the treatment (DAC)

		tmp_df = pulses_raw_df.loc[pulses_raw_df['module'] == 'DAC']

		if self.debug_mode:
			
			print('\n------------- DAC processing input -------------\n')
			display(tmp_df)

		for index, row in tmp_df.iterrows():

			if row['start'] > time_DAC[int(row['channel'])-1]:
				
				label = 'wait' + str(wait_color-int("D3D3D3", 16)+1)
				start = time_DAC[int(row['channel'])-1]
				stop = row['start']
				time = row['start'] - time_DAC[int(row['channel'])-1]
				module = row['module']
				Channel = 'DAC ch' + str(int(row['channel']))
				mode = 'wait'
				color_dict[str(wait_color)] = '#{0:06X}'.format(wait_color)
				color = '#{0:06X}'.format(wait_color)
				wait_color += 1
				param = row['param']
				ch_num = row['channel']
				
				pulses_df = pd.concat([pulses_df, pd.DataFrame.from_records([dict(label=label, start=start, stop=stop, time=time, module=module , Channel=Channel, mode=mode, color=str(color), param=param, ch_num=ch_num)])], ignore_index=True)

			start = row['start']
			length = row['length']
			stop = start + length


			label = index
			time = length
			module = row['module']
			Channel = 'DAC ch' + str(int(row['channel']))
			mode = row['mode']
			param = row['param']
			ch_num = row['channel']
			LUT = not(np.isnan(row['LUT']))
			start_pointer = row['starting_pointer']


			color = '#{0:06X}'.format(DAC_color)

			pulses_df = pd.concat([pulses_df, pd.DataFrame.from_records([dict(label=label, start=start, stop=stop, time=time, module=module, Channel=Channel, mode=mode, color=str(color), param=param, ch_num=ch_num, LUT=LUT, start_pointer=start_pointer)])], ignore_index=True)



			time_DAC[int(row['channel'])-1] = stop

			if stop>termination_time:

				termination_time = stop

		if self.debug_mode:

			print('\n------------- xxxxxxxxxxxxxxxxxxxxxxx -------------\n')
			display(pulses_df)
			print('\n------------- DAC processing finished -------------\n')


		# --- Beginning of the treatment (ADC)

		tmp_df = pulses_raw_df.loc[pulses_raw_df['module'] == 'ADC']

		if self.debug_mode:
			
			print('\n------------- ADC processing input -------------\n')
			display(tmp_df)

		for index, row in tmp_df.iterrows():

			if row['start'] > time_ADC[int(row['channel'])-1]:
			
				label = 'wait' + str(wait_color-int("D3D3D3", 16)+1)
				start = time_ADC[int(row['channel'])-1]
				stop = row['start']
				time = row['start'] - time_ADC[int(row['channel'])-1]
				module = row['module']
				Channel = 'ADC ch' + str(int(row['channel']))
				mode = 'wait'
				color_dict[str(wait_color)] = '#{0:06X}'.format(wait_color)
				color = '#{0:06X}'.format(wait_color)
				wait_color += 1
				param = row['param']
				ch_num = row['channel']
				
				pulses_df = pd.concat([pulses_df, pd.DataFrame.from_records([dict(label=label, start=start, stop=stop, time=time, module=module , Channel=Channel, mode=mode, color=str(color), param=param, ch_num=ch_num)])], ignore_index=True)

			start = row['start']
			length = row['length']
			stop = start + length
			label = index
			time = length
			module = row['module']
			Channel = 'ADC ch' + str(int(row['channel']))
			ch_num = row['channel']
			color = '#{0:06X}'.format(ADC_color)

			# fill the parameters only used for the 'IQ_table' mode 
			LUT = not(np.isnan(row['LUT']))
			ch_demod = row['demodulation_channels']
			param = row['param']
			mode = row['mode']
			start_pointer = 0 # shift in the starting pointer not implemented so far


			pulses_df = pd.concat([pulses_df, pd.DataFrame.from_records([dict(label=label, start=start, stop=stop, time=time, module=module, Channel=Channel, mode=mode, color=str(color), param=param, ch_num=ch_num,start_pointer=start_pointer, LUT=LUT, ch_demod=ch_demod)])], ignore_index=True)

			self.ADC_events = np.append(self.ADC_events,ch_demod)
			length_vec[int(row['channel'])-1].append(float(length)*1e-6*self.sampling_rate)

			time_ADC[int(row['channel'])-1] = stop

			if stop>termination_time:
			# update of the termination time used to compute if a wait will be added before the next pulse
				termination_time = stop

		if self.debug_mode:

			print('\n------------- xxxxxxxxxxxxxxxxxxxxxxx -------------\n')
			display(pulses_df)
			print('\n------------- ADC processing finished -------------\n')


		# --- Pulse sequences treatment solving inconsistency regarding stop time when n_rep>1 and more than one pulses is used. 

		# pulses_fin = pd.DataFrame()
		# for idx, ch in enumerate(channel_list):
		# 	pulses_ch =  pulses_df.loc[pulses_df['Channel']==ch]
		# 	pulses_ch = pulses_ch.sort_values('start')
		# 	for i in range(len(pulses_ch) - 1):
		# 		pulses_ch.iloc[i, pulses_ch.columns.get_loc('stop')] = pulses_ch.iloc[i + 1]['start']
		# 		pulses_ch.iloc[i, pulses_ch.columns.get_loc('time')] = pulses_ch.iloc[i]['stop'] - pulses_ch.iloc[i]['start']
		# 	pulses_fin = pd.concat([pulses_fin, pulses_ch], ignore_index=True)
		pulses_fin = pulses_df

		# --- Sequence plotting 

		if self.display_sequence:

			# --- Prepare channel list for display

			channel_set = list(set(list(pulses_df['Channel'])))
			channel_list = []

			for i in range(1,9):

				if 'DAC ch'+str(i) in channel_set:

					channel_list.append('DAC ch'+str(i))

				if 'ADC ch'+str(i) in channel_set:

					channel_list.append('ADC ch'+str(i))

			if self.global_loop[0]:

				channel_list.append('loop global')

			# --- Prepare the figure

			pulses_loop_list = []
			for idx, ch in enumerate(channel_list):

				if ch[0:4] != 'loop':

					pulses_loop =  pulses_fin.loc[pulses_fin['Channel']==ch]
					pulses_loop = pulses_loop.sort_values('start')

					pulses_loop_list.append(pulses_loop)

				elif ch == 'loop global':

					pulses_loop = loops_df.loc[loops_df['mode']=='global']
					pulses_loop.at[0,'stop'] = termination_time
					pulses_loop.at[0,'time'] = termination_time
					pulses_loop_list.append(pulses_loop)

			fig = make_subplots(rows=len(channel_list), cols=1, shared_xaxes=True)   

			if self.debug_mode:

				print('\n------------- Sequence display table -------------\n')         

			for idx, pulses_loop in enumerate(pulses_loop_list):

				if self.debug_mode:

					display(pulses_loop)

				fig.add_trace(go.Bar(x=pulses_loop.time, y=pulses_loop.Channel, orientation='h', text=pulses_loop.label, marker=dict(color=pulses_loop.color), 
										customdata=np.dstack((pulses_loop.label, pulses_loop.start, pulses_loop.stop))[0],
										hovertemplate='label: %{customdata[0]: s} <br>start: %{customdata[1]:.3f} <br>stop: %{customdata[2]:.3f}'), idx + 1, 1)

			if self.debug_mode:

				print('\n------------- Sequence display table -------------\n')

			fig.update_layout(showlegend=False,width=1600,height=200+80*idx,)
			fig.show()


		# store the length of the pulses and the channel used, it is use for data shaping in get_readout_pulse()
		self.length_vec = length_vec


		return pulses_fin



	def LUT_and_adress_filling(self, pulses_df, ADC=True, DAC=True):

		"""
		This function take the pulses sequence provided and store the pulses in the LUT of the DAC/ADC. 
		It also define the pointers that can be used. 

		For now three pointers can be defined in the ADC: start/loop/stop. 
		However, the start pointer cannot be shifted for now (it is set at the beginning of the pulses).
		Moreover, the loop pointer is set equal to the start pointer.


		Parameters:
		pulses_df -- panda DataFrame containing the pulse sequence. 
		ADC -- LUT and pointer address for the ADC (default True).
		DAC -- LUT and pointer address for the DAC (default True).

		Return: 
		DAC_pulses_pointer, ADC_pulses_pointer -- two lists containing the DAC and ADC pointers.

		"""

		DAC_pulses_array = [np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])]
		DAC_pulses_pointer = [[],[],[],[],[],[],[],[]]

		ADC_I_pulses_array = [np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])]
		ADC_Q_pulses_array = [np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])]
		ADC_pulses_pointer = [[[], [], []], [[], [], []], [[], [], []], [[], [], []], [[], [], []] , [[], [], []] , [[], [], []], [[], [], []]]


		if DAC:

			LUT_df = pulses_df.loc[(pulses_df['module'] == 'DAC') & (pulses_df['mode'] != 'wait')]

			if self.debug_mode:

				print('\n------------------------------- DAC LUT generation and upload ----------------------------------\n')
				display(LUT_df)

			event_time_list = list(dict.fromkeys(LUT_df['start']))
			event_time_list.sort()


			for event_time in event_time_list:

				event_df = LUT_df.loc[LUT_df['start'] == event_time]

				if self.debug_mode:

					print('\nCreating LUT map for events at ' + str(event_time))
					display(event_df)

				for index, row in event_df.iterrows():

					ch_num = int(row['ch_num'])

					if not(np.isnan(row['start_pointer'])):

						if row['LUT']:

							pulse_addr = round((len(DAC_pulses_array[ch_num-1]) + row['start_pointer']*1e-6*self.sampling_rate)/11)
							last_pointer = pulse_addr

							DAC_waveform = self.pulse_gen_SCPI(row['mode'],row['param'],row['time'],ch_num, module='DAC', pulse_label=row['label'])
							# adding pulse to waveform
							DAC_pulses_array[ch_num-1] = np.append(DAC_pulses_array[ch_num-1],DAC_waveform)

						else:
							
							pulse_addr =  last_pointer + round((row['start_pointer']*1e-6*self.sampling_rate)/8)
							
						DAC_pulses_pointer[ch_num-1].append(pulse_addr)

						if self.debug_mode:

							print('For pulse \"' + row['label'] + '\" LUT pointer at: ' + str(pulse_addr))
							print()

					else:

						log.error('Start pointer not found for ' + row['label'])


			for i in range(8):

				if len(DAC_pulses_array[i])>0:

					self.write('DAC:DATA:CH{}:CLEAR'.format(str(i+1)))

					if self.debug_mode:

						LUT_points = DAC_pulses_array[i].copy()
						LUT_points = self.DAC_amplitude_calib[i]*(LUT_points.reshape((-1, 11)).T[:-3].T.reshape((-1)))/(2**13)

						fig = plt.figure(figsize=(18,5))
						plt.plot(np.array(range(len(LUT_points)))/2000,LUT_points)
						plt.grid()
						plt.title('LUT for DAC ch '+str(i+1),fontsize = 18)
						# plt.legend(fontsize = 14)
						plt.show()

					DAC_SCPI_cmd = 'DAC:DATA:CH' + str(i+1) + ' 0,' + ','.join((DAC_pulses_array[i].astype(int)).astype(str)) 

					if self.debug_mode:

						print('\nDAC sequence for CH '+str(i+1)+':\n',DAC_SCPI_cmd[:40]+'...')

					log.info('Writing waveform for CH'+str(i+1)+'  \n')
					self.write(DAC_SCPI_cmd)

			log.info('Waveform processing complete' + '\n')

			if self.debug_mode:

				print('\n------------------------------- DAC LUT processing finished ----------------------------------\n')

		if ADC and self.acquisition_mode()=='INT':

			LUT_df = pulses_df.loc[(pulses_df['module'] == 'ADC') & (pulses_df['mode'] != 'wait')]

			if self.debug_mode:

				print('\n------------------------------- ADC LUT generation and upload ----------------------------------\n')
				display(LUT_df)

			ch_demod_list = LUT_df['ch_demod'].tolist()
			ch_demod_list = [item for sublist in ch_demod_list for item in sublist]
			ch_demod_list = set(ch_demod_list)

			for i in range(1,9):

				if i in ch_demod_list:

					if i==1:
						self.ADC1.status('ON')
						self.ADC1.fmixer(0)
					if i==2:
						self.ADC2.status('ON')
						self.ADC2.fmixer(0)
					if i==3:
						self.ADC3.status('ON')
						self.ADC3.fmixer(0)
					if i==4:
						self.ADC4.status('ON')
						self.ADC4.fmixer(0)
					if i==5: 
						self.ADC5.status('ON')
						self.ADC7.fmixer(0)
					if i==6:
						self.ADC6.status('ON')
						self.ADC6.fmixer(0)
					if i==7:
						self.ADC7.status('ON')
						self.ADC7.fmixer(0)
					if i==8:
						self.ADC8.status('ON')
						self.ADC8.fmixer(0)

					if self.debug_mode:

						print('Activating ADC channel ' + str(i) + ' for data acquisition.\n')


			event_time_list = list(dict.fromkeys(LUT_df['start']))
			event_time_list.sort()

			for event_time in event_time_list:

				event_df = LUT_df.loc[LUT_df['start'] == event_time]

				if self.debug_mode:

					print('\nCreating LUT map for events at ' + str(event_time))
					display(event_df)

				for index, row in event_df.iterrows():

					if row['mode']=='sin':

						if not(np.isnan(row['start_pointer'])):

							ch_num = int(row['ch_num'])
							ch_demod = row['ch_demod']

							for idx, chd in enumerate(ch_demod):

								if row['LUT']:

									if self.debug_mode:

										print('Uploading new LUT for ADC acquisition \"' + row['label'] + '\" with start pointer at: %i'%round(len(ADC_I_pulses_array[chd-1])/8))
										# print('Starting pointer set at: %i'%round(len(ADC_I_pulses_array[chd-1])/8))

									pulse_addr_start = len(ADC_I_pulses_array[chd-1])//8
									pulse_addr_loop = pulse_addr_start

									param_I = row['param'][idx]
									freq_demod = param_I['freq']

									# -- Error statement needs more explanation : marked for removal
									#
									# if 1e3%freq_demod!=0:
									# 	log.error('Demodulation frequency is not a multiple of the sampling_rate. \
									# 			   As a result, you are not using the loop pointer and the memory \
									# 			   usage of the ADC LUT is sub-optimal.')
									# 	time_vec = row['time']
									# else:
									# 	period_demod = int(1e3/freq_demod) # in ns
									# 	period_loop = int(1e9/self.sampling_rate*8) # in ns
									# 	time_vec = np.lcm(period_demod, period_loop)*1e-3 # in µs

									param_Q = param_I.copy()
									param_I['phase_offset'] = 0
									param_Q['phase_offset'] = np.pi/2

									waveform_I = self.pulse_gen_ADC_LUT(param_I,row['time'],chd)
									waveform_Q = self.pulse_gen_ADC_LUT(param_Q,row['time'],chd)

									if self.debug_mode:

										fig = plt.figure(figsize=(18,5))
										plt.plot(np.array(range(len(waveform_I)))/2000,waveform_I/(2**13),label='I')
										plt.plot(np.array(range(len(waveform_Q)))/2000,waveform_Q/(2**13),label='Q')
										plt.grid()
										plt.title('LUT for ADC pulse \"'+row['label']+'\"',fontsize = 18)
										plt.legend(fontsize = 14)
										plt.show()

									ADC_I_pulses_array[chd-1] = np.append(ADC_I_pulses_array[chd-1], waveform_I)
									ADC_Q_pulses_array[chd-1] = np.append(ADC_Q_pulses_array[chd-1], waveform_Q)

								else:

									log.error('ADC LUT sharing may have bugs and not supported yet.')

									# pulse_addr_start = ADC_pulses_pointer[chd-1][0][-1]
									# pulse_addr_stop = ADC_pulses_pointer[chd-1][2][-1]
									# pulse_addr_loop = ADC_pulses_pointer[chd-1][1][-1]



								pulse_addr_stop = len(ADC_I_pulses_array[chd-1])//8 - 1


								ADC_pulses_pointer[chd-1][0].append(pulse_addr_start)
								ADC_pulses_pointer[chd-1][1].append(pulse_addr_loop)
								ADC_pulses_pointer[chd-1][2].append(pulse_addr_stop)

								if self.debug_mode:
									print('For pulse ' + row['label'] + 'Saved start/loop/stop pointer: ', str(pulse_addr_start), str(pulse_addr_loop), str(pulse_addr_stop))
									print()

						else:

							log.error('Start pointer not found for ' + row['label'])

					else:

						log.error('Only sinusoidal demodulation is supported in current version.')

			if self.debug_mode:
				print('\nADC pointer array:')
				print(ADC_pulses_pointer)
				print()

			# -- Unclear reset, ADC memory is 128,000 points - marked for removal
			#
			# reset memory :
			#
			# for i in range(4):
			# 	empty_list = str([0]*16384)[1:-1]
			# 	empty_list = empty_list.replace(' ', '')
			# 	ADC_SCPI_cmd_I = 'ADC:I' + str(i+1) +' 0,' + empty_list
			# 	ADC_SCPI_cmd_Q = 'ADC:Q' + str(i+1) +' 0,' + empty_list
			# 	self.write(ADC_SCPI_cmd_I)
			# 	self.write(ADC_SCPI_cmd_Q)

			for i in range(4):

				if len(ADC_I_pulses_array[i])>0:

					if self.debug_mode:

						fig = plt.figure(figsize=(18,5))
						plt.plot(np.array(range(len(ADC_I_pulses_array[i])))/2000,ADC_I_pulses_array[i]/(2**13),label='I')
						plt.plot(np.array(range(len(ADC_Q_pulses_array[i])))/2000,ADC_Q_pulses_array[i]/(2**13),label='Q')
						plt.grid()
						plt.title('LUT for ADC channel(mixer) \"'+str(i+1)+'\"',fontsize = 18)
						plt.legend(fontsize = 14)
						plt.show()

					ADC_SCPI_cmd_I = 'ADC:I' + str(i+1) +' 0,' + ','.join((ADC_I_pulses_array[i].astype(int)).astype(str))
					ADC_SCPI_cmd_Q = 'ADC:Q' + str(i+1) +' 0,' + ','.join((ADC_Q_pulses_array[i].astype(int)).astype(str))

					if self.debug_mode:

						print('ADC sequence for I'+str(i+1)+':', ADC_SCPI_cmd_I[0:50])
						print('ADC sequence for Q'+str(i+1)+':', ADC_SCPI_cmd_Q[0:50])

					log.info('Writing waveform for I'+str(i+1)+'  \n')
					self.write(ADC_SCPI_cmd_I)

					log.info('Writing waveform for Q'+str(i+1)+'  \n')
					self.write(ADC_SCPI_cmd_Q)

			log.info('Waveform processing complete' + '\n')

		elif ADC and self.acquisition_mode()=='RAW':

			LUT_df = pulses_df.loc[(pulses_df['module'] == 'ADC') & (pulses_df['mode'] != 'wait')]

			ch_on_list = LUT_df['Channel'].tolist()

			for i in ch_on_list:

				if int(i[-1])==1:
					self.ADC1.status('ON')
					self.ADC1.fmixer(0)
				if int(i[-1])==2:
					self.ADC2.status('ON')
					self.ADC2.fmixer(0)
				if int(i[-1])==3:
					self.ADC3.status('ON')
					self.ADC3.fmixer(0)
				if int(i[-1])==4:
					self.ADC4.status('ON')
					self.ADC4.fmixer(0)
				if int(i[-1])==5: 
					self.ADC5.status('ON')
					self.ADC7.fmixer(0)
				if int(i[-1])==6:
					self.ADC6.status('ON')
					self.ADC6.fmixer(0)
				if int(i[-1])==7:
					self.ADC7.status('ON')
					self.ADC7.fmixer(0)
				if int(i[-1])==8:
					self.ADC8.status('ON')
					self.ADC8.fmixer(0)

				if self.debug_mode:

					print('Activating ADC channel ' + str(i) + ' for data acquisition.\n')

			ADC_pulses_pointer = None

		return DAC_pulses_pointer, ADC_pulses_pointer


	def process_sequencing_IQ_table(self):

		"""
		This function take the stored pulses sequence and create the SCPI command 
		that is sent to the RFSoC sequencer. 
		It also fills the ADC/DAC LUT via LUT_and_adress_filling()
		"""

		# --- Charging the pulse sequence and the event list

		pulses_df = self.pulses_sequence()
		seq_looping = self.seq_looping
		event_time_list = list(dict.fromkeys(pulses_df['start']))
		event_time_list.sort()
		# termination time used to close the sequence 
		termination_time = np.max((np.array(list(pulses_df['stop']))))

		# --- Initialization of the parameters
		
		# array containing the instructions played by the sequencer
		global_sequence = np.array([])
		global_sequence_info = []
		# load the pointer
		DAC_pulses_pointer, ADC_pulses_pointer = self.LUT_and_adress_filling(pulses_df)


		# lists of boolean storing the ADC/DAC state
		ADC_state = [0,0,0,0,0,0,0,0]
		DAC_state = [0,0,0,0,0,0,0,0]
		ADC_state_prev = [0,0,0,0,0,0,0,0]
		DAC_state_prev = [0,0,0,0,0,0,0,0]
		self.ADC_ch_active = np.array([0,0,0,0,0,0,0,0])

		# map for LUT ADC
		LUT_ADC_map = np.array(np.identity(4))*0

		# lists storing the ADC/DAC pointer index 
		pointer_dac = [0,0,0,0,0,0,0,0]
		pointer_adc = [0,0,0,0,0,0,0,0]
		# boolean checking if the repetition started
		rep_started = False
		# array used to store with physical ADC is send to which ADC LUT
		mux_state = np.zeros((4, 4), dtype='int')



		# # ---  If a repetitions is used it check the number of pulses per repetitions
		# nb_pulses_dac = np.zeros(8)
		# nb_pulses_adc = np.zeros(8)
		# pulses_counter_dac = np.ones(8)
		# pulses_counter_adc = np.ones(8)

		# pulses_rep_all = pulses_df.loc[(pulses_df['rep_nb']>1) & (pulses_df['mode']!='wait')]
		# pulses_rep_dac = pulses_rep_all.loc[(pulses_rep_all['module']=='DAC')]
		# pulses_rep_adc = pulses_rep_all.loc[(pulses_rep_all['module']=='ADC')]

		# # count the pulses number per dac 

		# if len(pulses_rep_dac)>0:

		# 	for index, row in pulses_rep_dac.iterrows():

		# 		ch_num = int(row['ch_num'])
		# 		nb_pulses_dac[ch_num - 1] +=1 

		# # same for the adc 

		# if len(pulses_rep_adc)>0:
		# 	for index, row in pulses_rep_adc.iterrows():
		# 		ch_demod = row['ch_demod']
		# 		for chd in ch_demod: 
		# 			nb_pulses_adc[chd - 1] +=1
		
		# # --- Start of the sequence filling 

		if self.debug_mode:

			print('\n\n\n------------------------------- Generating the sequence commands ----------------------------------\n')
			display(pulses_df.sort_values('start'))

		for event_time in event_time_list:

			if event_time>0:

				# add a waiting time any time it is needed 
				global_sequence = np.append(global_sequence,1)
				wait_time = int(round((event_time-event_time_prev)*250) - 1)
				global_sequence = np.append(global_sequence, wait_time)

				global_sequence_info.append('adding wait till the event at ' + str(event_time))

				if self.debug_mode:

					print(global_sequence_info[-1])
					print('Seq instruction : ',global_sequence[-2],global_sequence[-1])
					print()

			event_time_prev = event_time

			# take all events occurring at a given time
			tmp_df = pulses_df.loc[pulses_df['start'] == event_time]
			tmp_df = tmp_df.sort_values(by='module', ascending=False)

			if self.debug_mode:

					print('Creating sequence entry for events at ' + str(event_time))
					display(tmp_df)

			for index, row in tmp_df.iterrows():

				ch_num = int(row['ch_num'])

				if row['module'] == 'DAC':

					if row['mode'] != 'wait':

						# set the starting pointer to the correct address 
						global_sequence = np.append(global_sequence,4096+ch_num)
						global_sequence = np.append(global_sequence,DAC_pulses_pointer[ch_num-1][pointer_dac[ch_num-1]])

						global_sequence_info.append('Adding sequencer command to point to address of \"' + row['label'] + '\"')

						pointer_dac[ch_num-1] += 1

						if self.debug_mode:

							print(global_sequence_info[-1])
							print('Seq instruction : ',global_sequence[-2],global_sequence[-1])
							print()

						# the corresponding DAC is set to ON
						DAC_state[(ch_num-1)] = 1

						if self.debug_mode:

							print('DAC state is :', DAC_state)
							print()

					if row['mode'] == 'wait':

						# we switch off the corresponding DAC, the waiting command is set at the next iteration
						DAC_state[(ch_num-1)] = 0

						if self.debug_mode:

							print('DAC state is :', DAC_state)


				if row['module'] == 'ADC':

					ch_demod = row['ch_demod']
					ADC_ch = int(row['ch_num'])
					LUT_ADC_map_prev = LUT_ADC_map.copy()

					if row['mode'] != 'wait':

						# change the ADC to ON 

						for ch_demod_i in ch_demod:
							
							ADC_state[ch_demod_i-1] = 1

							if self.debug_mode:

								print('ADC state is :', ADC_state)
								print()

							# add the number of points the ADC LUT should take 
							global_sequence = np.append(global_sequence,4106 + ch_demod_i)
							global_sequence = np.append(global_sequence,int(row['time']*1e-6*self.sampling_rate))

						global_sequence_info.append('Adding sequencer command to set number of acquisition points to ' + str(int(row['time']*1e-6*self.sampling_rate)))

						if self.debug_mode:

							print(global_sequence_info[-1])
							print('Seq instruction : ',global_sequence[-2],global_sequence[-1])
							print()

						if self.acquisition_mode()=='INT':

							# --- Update the ADC to ADC LUT connection

							if self.debug_mode:

								print('ADC 1-4 LUT Map : ')
								pprint(LUT_ADC_map)

							for k in list(map(int,ch_demod)):

								LUT_ADC_map[ADC_ch-1][k-1] = 1

							# --- test validity of ADC map

							mux_state = np.sum(LUT_ADC_map, axis=0)
							mux_state_valid = True

							for k in range(len(mux_state)):

								if mux_state[k]>1:

									mux_state_valid = False
									log.error('mux state validity issue, are you trying to use one mixer for multiple demods?')

							for k in range(2,4):

								for l in range(0,2):

									if LUT_ADC_map[k,l] != 0:

										mux_state_valid = False
										log.error('mux state validity issue, invalid channel mixing!')

							if self.debug_mode:

								print('Updated ADC 1-4 LUT Map : ')
								pprint(LUT_ADC_map)
								print('mux state validity registered as: ' + str(mux_state_valid))


							# --- Sequence filling for ADC type command
							if mux_state_valid and not(np.array_equal(LUT_ADC_map_prev,LUT_ADC_map)):

								mux_config_for_ch = self._mux_config_matrix[ADC_ch-1]

								for k in range(4):

									chd = k+1

									if LUT_ADC_map[ADC_ch-1][k] == 1:

										param_id_1 = 4128 + 2*k
										param_id_2 = 4129 + 2*k

										if k<2:

											mux_b31 = '0'
											mux_b30 = '0'
											mux_b30_b29 = mux_b30 + str(mux_config_for_ch[k])
											mux_b28 = '1'

										else:

											mux_b31 = '0'
											mux_b30_b29 = bin(int(mux_config_for_ch[k]))[2:].zfill(2)
											mux_b28 = '1'

										mux_b27_to_b14 = '0'*14

										iqram_addr_start = bin(ADC_pulses_pointer[chd - 1][0][pointer_adc[chd - 1]])[2:].zfill(14)
										iqram_addr_loop = bin(ADC_pulses_pointer[chd - 1][1][pointer_adc[chd - 1]])[2:].zfill(14)
										iqram_addr_stop = bin(ADC_pulses_pointer[chd - 1][2][pointer_adc[chd - 1]])[2:].zfill(14)

										param_val_1 = int(mux_b31 + mux_b30_b29 + mux_b28 + mux_b27_to_b14 + iqram_addr_start, 2)
										param_val_2 = int('0'*2 + iqram_addr_loop + '0'*2 + iqram_addr_stop, 2)
										
										global_sequence = np.append(global_sequence, param_id_1)
										global_sequence = np.append(global_sequence, param_val_1)
										global_sequence = np.append(global_sequence, param_id_2)
										global_sequence = np.append(global_sequence, param_val_2)

										global_sequence_info.append('Adding sequencer command to configure mux for ADC ch ' + str(ADC_ch) + ' and demod_ch ' + str(k))
										global_sequence_info.append('Adding sequencer command to configure mux for ADC ch ' + str(ADC_ch) + ' and demod_ch ' + str(k))

										if self.debug_mode:

											print(global_sequence_info[-1])
											print('Seq instruction : ',global_sequence[-4],global_sequence[-3])
											print('Seq instruction : ',global_sequence[-2],global_sequence[-1])
											print()

										# the pointer index of the given channel is updated or not depending on the situation
										pointer_adc[k] += 1

						elif self.acquisition_mode()=='RAW':

							param_id_1 = 4128 + 2*(ADC_ch-1)

							mux_b31 = '0'
							mux_b30_b29 = '00'
							mux_b28 = '0'

							mux_b27_to_b14 = '0'*14
							mux_b13_to_b0 = '0'*14

							param_val_1 = int(mux_b31 + mux_b30_b29 + mux_b28 + mux_b27_to_b14 + mux_b13_to_b0, 2)

							global_sequence = np.append(global_sequence, param_id_1)
							global_sequence = np.append(global_sequence, param_val_1)

							global_sequence_info.append('Adding sequencer command to turn off mux for ADC ch ' + str(ADC_ch))

							if self.debug_mode:

								print(global_sequence_info[-1])
								print('Seq instruction : ',global_sequence[-4],global_sequence[-3])
								print('Seq instruction : ',global_sequence[-2],global_sequence[-1])
								print()

					elif row['mode'] == 'wait':

						ADC_state[ADC_ch-1] = 0

						if self.debug_mode:

							print('Wait pulse - ADC state is :', ADC_state)
							print()

			# ---  Update of the ADC and DAC state at every time step 

			if DAC_state != DAC_state_prev or ADC_state != ADC_state_prev:

				if self.debug_mode:

					print('ADC state updated from %s to %s'%(ADC_state_prev, ADC_state))
					print('DAC state updated from %s to %s'%(DAC_state_prev, DAC_state))

				# update the ADC state 
				bin_adc_cmd = ''.join(map(str,reversed(ADC_state)))

				# update the DAC state
				bin_dac_cmd = ''

				for i in reversed(range(8)):

					if DAC_state[i] != DAC_state_prev[i]:

						if DAC_state[i] == 1:

							bin_dac_cmd += '011'
							DAC_state_prev[i] = 1

						else:

							bin_dac_cmd += '001'
							DAC_state_prev[i] = 0

					else:

						bin_dac_cmd +='000'

				# add the two bit strings and add the command to update states 
				bin_ADC_DAC_state = bin_adc_cmd + bin_dac_cmd

				global_sequence = np.append(global_sequence,4096)
				global_sequence = np.append(global_sequence,int(bin_ADC_DAC_state,2))

				global_sequence_info.append('Updating ADC DAC state')

				ADC_state_prev = ADC_state.copy()

				if self.debug_mode:

					print('Bit string of the DAC state: ',bin_dac_cmd)
					print('Bit string of the ADC state: ',bin_adc_cmd)
					print()
					print('Seq instruction : ',global_sequence[-2],global_sequence[-1])
					print()

			else :

				if self.debug_mode:

					print('No ADC or DAC update at this state')
					print('ADC previous: %s  / ADC step : %s'%(ADC_state_prev, ADC_state))
					print('DAC previous: %s  / DAC step : %s'%(DAC_state_prev, DAC_state))
					print()



		# --- Add a last waiting time if needed

		wait_term = int(round((termination_time - event_time)*250))-1
		global_sequence = np.append(global_sequence,1)
		global_sequence = np.append(global_sequence,wait_term)

		global_sequence_info.append('Terminate by wait of : %f' %(wait_term + 1))

		if self.debug_mode:

			print(global_sequence_info[-1])
			print('Seq instruction : ',global_sequence[-2],global_sequence[-1])

		# Switch off all the DAC/ADC 
		global_sequence = np.append(global_sequence,4096)
		global_sequence = np.append(global_sequence,int('00000000001001001001001001001001',2))

		global_sequence_info.append('Reached end - switch off all the DAC/ADC')

		# --- Set the Acquisition mode 
		if self.acquisition_mode() == 'RAW':
			acq_mode = 0
		elif self.acquisition_mode() == 'INT':
			acq_mode = 286331153
		else:
			log.error('Invalid acquisition mode\n')

		# if global loop is programmed we add it to the sequence 
		try:
			global_rep = int(seq_looping['rep'].where(seq_looping['mode'] == 'global'))
		except:
			global_rep = 1

		if global_rep > 1:

			global_sequence_str = 'SEQ 0,1,9,4106,' + str(acq_mode) + ',258,' + str(int(global_rep-1)) + ',' + ','.join((global_sequence.astype(int)).astype(str))  + ',514,0,0,0'

		else :

			global_sequence_str = 'SEQ 0,1,9,4106,' + str(acq_mode) + ',' + ','.join((global_sequence.astype(int)).astype(str)) + ',0,0'

		if self.debug_mode:

			print('Sequence programmer command: ',global_sequence_str[:50])

		# --- Send the SCPI command 
		log.info('Writing global sequence' + '\n')
		self.write(global_sequence_str)

		# Update class variables 
		self.n_points_total = global_rep
		self.ADC_ch_active = ADC_state
		self.global_sequence_str  = global_sequence_str
		self.global_sequence  = global_sequence
		self.global_sequence_info = global_sequence_info



	def pulse_gen_ADC_LUT(self,param,duration,ch):

		"""
		This functions separates the pulse generation for ADC LUTs, 
		disentangling the calibration of channel and module specific 
		calibration.
		"""

		period = 1./self.sampling_rate
		time_vec = np.arange(period,duration*1e-6+period/2,period)

		wavepoints = ((2**13)*self.ADC_amplitude_calib[ch-1]*param['amp'] - 1)*np.sin(-param['phase_offset'] + 2*np.pi*param['freq']*1e6*time_vec)

		idx_of = np.nonzero(((wavepoints > 8191) | (wavepoints < -8192)))[0]
		if len(np.nonzero(wavepoints[idx_of] != 16381)[0])>0:
			if module =='ADC':
				log.error('Error when filling the ADC memory: maximal amplitude is over the resolution') 
			elif module == 'ADC':
				log.error('Error when filling the ADC memory: maximal amplitude is over the resolution') 

		return wavepoints



	def pulse_gen_SCPI(self,mode,param,duration,ch, module, pulse_label='unknown'):

		"""
		This function convert the pulses parameters into a 1D array of Volt amplitude in bit scale.
		This is used in LUT_and_adress_filling() to create the SCPI command. 

		Parameters: 
		function -- (str) function to generate, can be 'sin+sin' for a sum of sin
		'sin' for a sin, 'trigger' for a trigger and 'DC' for a DC channel 
		param --  list containing the function parameters
		duration -- (float) pulse duration 
		ch --  (int) LUT channel
		module -- (str) 'DAC' or 'ADC'

		return: 
		wavepoints -- 1D array of points that will be stored in the LUT 
		"""

		period = 1./self.sampling_rate
		time_vec = np.arange(period,duration*1e-6+period/2,period)
		# stop vector used to separate different pulses in a channel 
		stop_vector = np.array([0,0,0,0,0,0,0,0,0,0,16383])

		if mode == 'sin+sin':

			wavepoints1 = ((2**13)*self.DAC_amplitude_calib[ch-1]*param['amp1']*np.sin(-param['phase_offset1'] - 1) + 2*np.pi*param['freq1']*1e6*time_vec)
			wavepoints2 = ((2**13)*self.DAC_amplitude_calib[ch-1]*param['amp2']*np.sin(-param['phase_offset2'] - 1) + 2*np.pi*param['freq2']*1e6*time_vec)

			wavepoints = ((2**13)*self.DAC_amplitude_calib[ch-1]*param['dc_offset'] - 1)+ wavepoints1 + wavepoints2

			# check if some points are above the limit imposed by the 13 bits
			idx_of = np.nonzero(((wavepoints > 8191) | (wavepoints < -8192)))[0]
			if len(np.nonzero(wavepoints[idx_of] != 16381)[0])>0:

				if module =='DAC':
					log.error('Error when filling the DAC memory: maximal amplitude is over the resolution') 
				elif module == 'ADC':
					log.error('Error when filling the ADC memory: maximal amplitude is over the resolution') 


			if self.debug_mode and self.debug_mode_plot_waveforms:

				fig = plt.figure(figsize=(18,5))
				plt.plot(np.array(range(len(wavepoints1)))/2000,wavepoints1/(2**13))
				plt.grid()
				plt.title('LUT for DAC pulse \"'+pulse_label+'\" sin1',fontsize = 18)
				# plt.legend(fontsize = 14)
				plt.show()

				fig = plt.figure(figsize=(18,5))
				plt.plot(np.array(range(len(wavepoints2)))/2000,wavepoints2/(2**13))
				plt.grid()
				plt.title('LUT for DAC pulse \"'+pulse_label+' sin2\"',fontsize = 18)
				# plt.legend(fontsize = 14)
				plt.show()

				fig = plt.figure(figsize=(18,5))
				plt.plot(np.array(range(len(wavepoints)))/2000,wavepoints/(2**13))
				plt.grid()
				plt.title('LUT for DAC pulse \"'+pulse_label+' total\"',fontsize = 18)
				# plt.legend(fontsize = 14)
				plt.show()


			if module == 'DAC':

				# --- When using 'DAC' every 8 points you insert trigger bits that define the trigger states 

				# adding zeros to make length multiple of 8
				wavepoints = np.append(wavepoints,np.zeros(len(time_vec)%8))
				trig_rep_len = int(len(wavepoints)/8)
				# adding trigger vectors
				wavepoints = np.concatenate((wavepoints.reshape(trig_rep_len,8), np.array(np.zeros(trig_rep_len))[:,None]),axis=1)
				wavepoints = np.concatenate((wavepoints, np.array(np.zeros(trig_rep_len))[:,None]),axis=1)

				# adding repetation (0 for once)
				wavepoints = np.concatenate((wavepoints, np.array(np.zeros(trig_rep_len))[:,None]),axis=1)
				# convert to 1D array
				wavepoints = wavepoints.reshape(1,len(wavepoints)*len(wavepoints[0]))

				wavepoints =  wavepoints[0]
				# wavepoints = np.append(wavepoints,stop_vector)

				



		elif mode == 'sin':

			wavepoints = ((2**13)*self.DAC_amplitude_calib[ch-1]*param['dc_offset'] - 1) + ((2**13)*self.DAC_amplitude_calib[ch-1]*param['amp'] - 1)*np.sin(-param['phase_offset'] + 2*np.pi*param['freq']*1e6*time_vec)

			idx_of = np.nonzero(((wavepoints > 8191) | (wavepoints < -8192)))[0]
			if len(np.nonzero(wavepoints[idx_of] != 16381)[0])>0:
				if module =='DAC':
					log.error('Error when filling the DAC memory: maximal amplitude is over the resolution') 
				elif module == 'ADC':
					log.error('Error when filling the ADC memory: maximal amplitude is over the resolution') 

			if self.debug_mode:
				
				fig = plt.figure(figsize=(18,5))
				plt.plot(np.array(range(len(wavepoints)))/2000,wavepoints/(2**13))
				plt.grid()
				plt.title('LUT for DAC pulse \"'+pulse_label+'\"',fontsize = 18)
				# plt.legend(fontsize = 14)
				plt.show()

			if module == 'DAC':

				# adding zeros to make length multiple of 8
				wavepoints = np.append(wavepoints,np.zeros(len(time_vec)%8))
				trig_rep_len = int(len(wavepoints)/8)
				# adding trigger vectors
				wavepoints = np.concatenate((wavepoints.reshape(trig_rep_len,8), np.array(np.zeros(trig_rep_len))[:,None]),axis=1)
				wavepoints = np.concatenate((wavepoints, np.array(np.zeros(trig_rep_len))[:,None]),axis=1)

				# adding repetation (0 for once)
				wavepoints = np.concatenate((wavepoints, np.array(np.zeros(trig_rep_len))[:,None]),axis=1)
				# convert to 1D array
				wavepoints = wavepoints.reshape(1,len(wavepoints)*len(wavepoints[0]))

				wavepoints =  wavepoints[0]
				# wavepoints = np.append(wavepoints,stop_vector)




		elif mode == 'trigger':

			wavepoints = np.zeros_like(time_vec) # only zeros

			# adding zeros to make length multiple of 8
			wavepoints = np.append(wavepoints,np.zeros(len(time_vec)%8))

			trig_rep_len = int(len(wavepoints)/8)
			# adding trigger vectors
			wavepoints = np.concatenate((wavepoints.reshape(trig_rep_len,8), np.array(np.ones(trig_rep_len))[:,None]),axis=1) # marker 1 in UP position
			wavepoints = np.concatenate((wavepoints, np.array(np.zeros(trig_rep_len))[:,None]),axis=1)

			# adding repetation (0 for once)
			wavepoints = np.concatenate((wavepoints, np.array(np.zeros(trig_rep_len))[:,None]),axis=1)

			# convert to 1D array
			wavepoints = wavepoints.reshape(1,len(wavepoints)*len(wavepoints[0]))

			wavepoints =  wavepoints[0]
			wavepoints = np.append(wavepoints,stop_vector)


		elif mode == 'DC':

			amplitude = (2**13)*self.DAC_amplitude_calib[ch-1]*param['amp']

			total_nop = int(duration * 1e-6 / period)
			# we want the decomposition total_nop = 8 * ( (N_segments-1)*16381 + remainder_loops ) + remainder_nop
			remainder_nop = total_nop % 8
			N_loops = total_nop // 8

			N_segments = N_loops // 16381    # we need to loop N_segments full segments (8 points)
			remainder_loops = N_loops % 16381    # (N_segments - 1) will be looped 16381 times and one will be looped the remainder_loops times
			if remainder_loops:
				N_segments += 1

			wavepoints = amplitude * np.ones(8*N_segments + remainder_nop) # reduced number of points used thanks to loops
			# adding zeros to make length multiple of 8
			if remainder_nop:
				wavepoints = np.append(wavepoints,np.zeros(8 - remainder_nop))

			trig_rep_len = int(len(wavepoints)/8)
			# adding trigger vectors
			wavepoints = np.concatenate((wavepoints.reshape(trig_rep_len,8), np.array(np.zeros(trig_rep_len))[:,None]),axis=1)
			wavepoints = np.concatenate((wavepoints, np.array(np.zeros(trig_rep_len))[:,None]),axis=1)

			# add repetitions (remainder_loops for first segment, 0 for last segment, 16381 for the rest of them)
			if remainder_nop:
				rep_arr = 16381 * np.ones(N_segments + 1)
				rep_arr[-1] = 0
			else:
				rep_arr = 16381 * np.ones(N_segments)
			if remainder_loops:
				rep_arr[0] = remainder_loops
			wavepoints = np.concatenate((wavepoints, rep_arr[:,None]),axis=1)

			# convert to 1D array
			wavepoints = wavepoints.reshape(1,len(wavepoints)*len(wavepoints[0]))

			wavepoints = wavepoints[0]
			wavepoints = np.append(wavepoints,stop_vector)


		else:

			log.error('Wrong waveform mode: ',mode)


		if len(wavepoints) > 128000: 
			if module == 'DAC':
				log.error('Error when filling the DAC memory : to many points, maximal number is 128000 while you are asking for %i'%len(wavepoints))
			if module == 'ADC':
				log.error('Error when filling the ADC memory : to many points, maximal number is 128000 while you are asking for %i'%len(wavepoints))

		return wavepoints


	def start_play(self):
		"""
		This function can be used while debugging the driver. 
		It asks three times for the data and print the output
		list and it length
		"""

		self.write("SEQ:STOP")
		time.sleep(0.1)
		self.write("SEQ:START")
		time.sleep(0.1)
		for k in range(3):
			try:
				 r = self.visa_handle.query_binary_values('OUTPUT:DATA?', datatype="h",
							is_big_endian=False)
			except:
				 r=[]
			if self.debug_mode:
				print(len(r), r)


	def stop_play(self):

		self.write("SEQ:STOP")


	def reset_PLL(self, force_reset=False):

		perform_PLL_reset = False

		if force_reset:

			perform_PLL_reset = True

		else:

			pll_status = self.ask('PLLINIT?')

			if pll_status == 0:

				perform_PLL_reset = True

		if perform_PLL_reset:

			self.write("DAC:RELAY:ALL 0")
			self.write("PLLINIT")
			time.sleep(5)
			self.write("DAC:RELAY:ALL 1")
			print('PLL reset complete.')

		else:

			print('PLL already initialized.')


	def reset_output_data(self):
		try:
			r = self.visa_handle.query_binary_values('OUTPUT:DATA?', datatype="h",
				is_big_endian=False)
		except:
			r=[]


	def get_readout_pulse(self):

		'''
		 This function reformat the data reading the header contents.
		'''

		self.reset_output_data()
		mode = self.acquisition_mode()
		n_rep = int(self.n_points_total)
		length_vec = self.length_vec
		ch_vec = self.ADC_events
		N_adc_events = len(ch_vec) # to be discussed with Arpit
		# N_adc_events = len(np.unique(ch_vec))
		if self.debug_mode:
			len_data_all = 0
		#print(length_vec)
		n_pulses = int(len(ch_vec)/len(np.unique(ch_vec))) # will need updating for general case
		# n_pulses = len(length_vec[0]) # to be discussed with Arpit


		ch_active = self.ADC_ch_active

		if mode == 'INT':

			data_unsorted = {}
			count_meas = 0
			empty_packet_count = 0
			run_num = 0

			getting_valid_dataset = True

			if self.display_IQ_progress:

				self.display_IQ_progress_bar = IntProgress(min=0, max=n_rep) # instantiate the bar
				display(self.display_IQ_progress_bar) # display the bar

			self.write("SEQ:STOP")
			time.sleep(0.1)
			self.write("SEQ:START")
			time.sleep(0.1)

			while getting_valid_dataset:

				while (count_meas//(16*N_adc_events))<n_rep:

					if self.loop_time: 
						a = datetime.datetime.now()

					try:
						r = self.visa_handle.query_binary_values('OUTPUT:DATA?', datatype="h",
							is_big_endian=False)
					except:
						r=[]

					if self.debug_mode: 
						if len(r)>1:
							len_data_all +=len(r)
						# print(len_data_all)

					if self.loop_time: 
						b = datetime.datetime.now()
						print('\nget data: ',b-a)

					if r == 'ERR':

						log.error('rfSoC: Instrument returned ERR!')

						# reset measurement
						data_unsorted = {}
						count_meas = 0
						run_num = 0
						self.write("SEQ:STOP")
						time.sleep(2)
						while True:
							try:
								junk = self.visa_handle.query_binary_values('OUTPUT:DATA?', datatype="h",
							is_big_endian=False)
							except:
								junk = []
							time.sleep(0.1)
							if junk == []:
								break
						junk = []
						self.write("SEQ:START")
						time.sleep(0.1)

						continue

					elif len(r)>1:

						data_unsorted['{}'.format(run_num)] = r
						r_size = len(r)
						count_meas += r_size
						if self.display_IQ_progress:
							self.display_IQ_progress_bar.value = count_meas//(16*N_adc_events)
						run_num += 1
						time.sleep(0.01)


					elif r==[]: # new empty packet?

						time.sleep(0.1)
						self.write("SEQ:STOP")
						try:
							junk = self.visa_handle.query_binary_values('OUTPUT:DATA?', datatype="h",
							is_big_endian=False)
						except:
							junk = []
						time.sleep(0.1)
						if junk ==[]:
							break
						junk = []

						self.write("SEQ:START")
						time.sleep(0.1)
						continue

				if count_meas//(16*N_adc_events) == n_rep:

					getting_valid_dataset = False

				elif self.debug_bypass_count_check:

					print('Count bypass triggered : ',count_meas,N_adc_events,n_rep)
					getting_valid_dataset = False

				else:

					log.error('Data corruption: rfSoC did not send all data points({}/'.format(count_meas//(16*N_adc_events))+str(n_rep)+').')

					print(count_meas)
					print(N_adc_events)

					# reset measurement
					data_unsorted = {}
					count_meas = 0
					empty_packet_count = 0
					self.write("SEQ:STOP")
					time.sleep(2)
					while True:
						try:
							junk = self.visa_handle.query_binary_values('OUTPUT:DATA?', datatype="h",
							is_big_endian=False)
						except:
								junk = []

						time.sleep(0.1)
						if junk == []:
							break
					junk = []
					self.write("SEQ:START")
					time.sleep(0.1)

			self.write("SEQ:STOP")

			'''
				Process data
			'''

			# data_unsorted = list(chain.from_iterable(data_unsorted.values()))
			data_unsorted = functools.reduce(operator.iconcat, list(data_unsorted.values()), [])
			data_unsorted = np.array(data_unsorted,dtype=int)

			# separate header from IQ data
			raw_IQ_data_dump_header = data_unsorted.reshape(int(len(data_unsorted)/8),8)[0::2]
			raw_I_data_dump_data = data_unsorted.reshape(int(len(data_unsorted)/8),8)[1::2][:,0:4]
			raw_Q_data_dump_data = data_unsorted.reshape(int(len(data_unsorted)/8),8)[1::2][:,4:8]

			# extract channel number from first byte of header
			ch_num = (raw_IQ_data_dump_header%256).T[0]

			# extract number of accumulated points for normalization from 3rd to 6th byte of header
			num_points = np.frombuffer(np.stack((raw_IQ_data_dump_header.T[1], raw_IQ_data_dump_header.T[2]), axis=1).astype('int16').tobytes(), dtype=np.int_)

			# vectors indicating channel that the data originated from
			ch_1 = ch_num*(ch_num == np.ones(len(ch_num)))
			ch_2 = ch_num*(ch_num == 2*np.ones(len(ch_num)))/2
			ch_3 = ch_num*(ch_num == 3*np.ones(len(ch_num)))/3
			ch_4 = ch_num*(ch_num == 4*np.ones(len(ch_num)))/4
			ch_5 = ch_num*(ch_num == 5*np.ones(len(ch_num)))/5
			ch_6 = ch_num*(ch_num == 6*np.ones(len(ch_num)))/6
			ch_7 = ch_num*(ch_num == 7*np.ones(len(ch_num)))/7
			ch_8 = ch_num*(ch_num == 8*np.ones(len(ch_num)))/8

			# data conversion form four 16 bit integers to one 64 bit longlong (*(0.3838e-3/16))
			I_all_data = 2 + np.frombuffer(raw_I_data_dump_data.astype('int16').tobytes(), dtype=np.longlong)*0.3838e-3/(16*num_points)
			Q_all_data = 2 + np.frombuffer(raw_Q_data_dump_data.astype('int16').tobytes(), dtype=np.longlong)*0.3838e-3/(16*num_points)


			# --- may be adapted for more advanced data shaping 

			I = [((I_all_data*ch_1)[I_all_data*ch_1!=0]-2).reshape(n_rep*ch_active[0],n_pulses).T,
				 ((I_all_data*ch_2)[I_all_data*ch_2!=0]-2).reshape(n_rep*ch_active[1],n_pulses).T,
				 ((I_all_data*ch_3)[I_all_data*ch_3!=0]-2).reshape(n_rep*ch_active[2],n_pulses).T,
				 ((I_all_data*ch_4)[I_all_data*ch_4!=0]-2).reshape(n_rep*ch_active[3],n_pulses).T,
				 ((I_all_data*ch_5)[I_all_data*ch_5!=0]-2).reshape(n_rep*ch_active[4],n_pulses).T,
				 ((I_all_data*ch_6)[I_all_data*ch_6!=0]-2).reshape(n_rep*ch_active[5],n_pulses).T,
				 ((I_all_data*ch_7)[I_all_data*ch_7!=0]-2).reshape(n_rep*ch_active[6],n_pulses).T,
				 ((I_all_data*ch_8)[I_all_data*ch_8!=0]-2).reshape(n_rep*ch_active[7],n_pulses).T]
			Q = [((Q_all_data*ch_1)[Q_all_data*ch_1!=0]-2).reshape(n_rep*ch_active[0],n_pulses).T,
				 ((Q_all_data*ch_2)[Q_all_data*ch_2!=0]-2).reshape(n_rep*ch_active[1],n_pulses).T,
				 ((Q_all_data*ch_3)[Q_all_data*ch_3!=0]-2).reshape(n_rep*ch_active[2],n_pulses).T,
				 ((Q_all_data*ch_4)[Q_all_data*ch_4!=0]-2).reshape(n_rep*ch_active[3],n_pulses).T,
				 ((Q_all_data*ch_5)[Q_all_data*ch_5!=0]-2).reshape(n_rep*ch_active[4],n_pulses).T,
				 ((Q_all_data*ch_6)[Q_all_data*ch_6!=0]-2).reshape(n_rep*ch_active[5],n_pulses).T,
				 ((Q_all_data*ch_7)[Q_all_data*ch_7!=0]-2).reshape(n_rep*ch_active[6],n_pulses).T,
				 ((Q_all_data*ch_8)[Q_all_data*ch_8!=0]-2).reshape(n_rep*ch_active[7],n_pulses).T]



			# I = [((I_all_data*ch_1)[I_all_data*ch_1!=0]-2),
			# 	 ((I_all_data*ch_2)[I_all_data*ch_2!=0]-2),
			# 	 ((I_all_data*ch_3)[I_all_data*ch_3!=0]-2),
			# 	 ((I_all_data*ch_4)[I_all_data*ch_4!=0]-2),
			# 	 ((I_all_data*ch_5)[I_all_data*ch_5!=0]-2),
			# 	 ((I_all_data*ch_6)[I_all_data*ch_6!=0]-2),
			# 	 ((I_all_data*ch_7)[I_all_data*ch_7!=0]-2),
			# 	 ((I_all_data*ch_8)[I_all_data*ch_8!=0]-2)]
			# Q = [((Q_all_data*ch_1)[Q_all_data*ch_1!=0]-2),
			# 	 ((Q_all_data*ch_2)[Q_all_data*ch_2!=0]-2),
			# 	 ((Q_all_data*ch_3)[Q_all_data*ch_3!=0]-2),
			# 	 ((Q_all_data*ch_4)[Q_all_data*ch_4!=0]-2),
			# 	 ((Q_all_data*ch_5)[Q_all_data*ch_5!=0]-2),
			# 	 ((Q_all_data*ch_6)[Q_all_data*ch_6!=0]-2),
			# 	 ((Q_all_data*ch_7)[Q_all_data*ch_7!=0]-2),
			# 	 ((Q_all_data*ch_8)[Q_all_data*ch_8!=0]-2)]

		elif mode == 'RAW':

			self.reset_output_data()

			#for now we consider only the one same type of acq on all adc
			mode=self.acquisition_mode.get()

			for index in range(8):

				length_vec[index] = np.unique(length_vec[index])

			getting_valid_dataset = True

			while getting_valid_dataset:

				adcdataI = [[],[],[],[],[],[],[],[]]
				adcdataQ = [[],[],[],[],[],[],[],[]]

				rep=[]

				keep_trying = True
				nb_try = 0

				self.write("SEQ:STOP")
				time.sleep(0.1)
				self.write("SEQ:START")
				time.sleep(0.1)
				i = 0

				while keep_trying:
					try:
						r = self.visa_handle.query_binary_values('OUTPUT:DATA?', datatype="h",
							is_big_endian=False)
					except:
						r=[]


					if self.debug_mode: 
						if len(r)>1:
							len_data_all +=len(r)
						print(len_data_all)

					if r == 'ERR':

						log.error('rfSoC: Instrument returned ERR!')

						# reset measurement
						rep = []
						empty_packet_count = 0
						self.write("SEQ:STOP")
						time.sleep(2)
						while True:
							try:
								junk = self.visa_handle.query_binary_values('OUTPUT:DATA?', datatype="h",
									is_big_endian=False)
							except:
								junk=[]

							# print(junk)
							time.sleep(0.1)
							if junk == []:
								break
						junk = []
						nb_try +=1
						self.write("SEQ:START")
						time.sleep(0.1)

						continue

					if len(r)>1:

						rep = rep+r
						nb_try = 0

					elif r==[]:
						try:
							junk = self.visa_handle.query_binary_values('OUTPUT:DATA?', datatype="h",
							is_big_endian=False)
						except:
								junk = []
						if junk == []:
							junk = []
							keep_trying = False

				self.write("SEQ:STOP")
				nb_try = 0

				i=0
				TSMEM=0
				while (i + 8 )<= len(rep) : # at least one header left

					entete = np.array(rep[i:i+8])
					X =entete.astype('int16').tobytes()
					V = X[0]-1 # channel (1 to 8)
					DSPTYPE = X[1]
					#N does not have the same meaning depending on DSTYPE
					N = struct.unpack('I',X[2:6])[0]
					#number of acquisition points in continuous
					#depends on the point length
					NpCont = X[7]*256 + X[6]
					TS= struct.unpack('Q',X[8:16])[0]

					# print the header for each packet
					# print("Channel={}; N={}; DSP_type={}; TimeStamp={}; Np_Cont={}; Delta_TimeStamp={}".format(V,N,DSPTYPE,TS,NpCont,TS-TSMEM))

					TSMEM=TS

					iStart=i+8
					# if not in continuous acq mode
					if ((DSPTYPE &  0x2)!=2):
						# raw adcdata for each Np points block
						if ((DSPTYPE  &  0x1)==0):
							Np=N
							adcdataI[V]=np.concatenate((adcdataI[V], np.right_shift(rep[iStart:iStart+Np],4)*0.3838e-3))

						#in the accumulation mode, only 1 I and Q point even w mixer OFF
						#mixer ON or OFF
						if ((DSPTYPE  & 0x01)==0x1):
							Np=8
							D=np.array(rep[iStart:iStart+Np])
							X = D.astype('int16').tobytes()

							#I  dvided N and 2 bcse signed 63 bits aligned to the left
							# mod div by 4 to fix amplitude -Arpit, Martina
							I=  struct.unpack('q',X[0:8])[0]*(0.3838e-3)/(N*2*4)
							Q=  struct.unpack('q',X[8:16])[0]*(0.3838e-3)/(N*2*4)

							# print the point
							# print("I/Q:",I,Q,"Amplitude:",np.sqrt(I*I+Q*Q),"Phase:",180*np.arctan2(I,Q)/np.pi)

							adcdataI[V]=np.append(adcdataI[V], I)
							adcdataQ[V]=np.append(adcdataQ[V], Q)


					#in our case we dont need the continuous mode for now
					# continuoous acquisition mode with accumulation (reduce the flow of data)
					elif ((DSPTYPE &  0x3)==0x3):
						# mixer OFF : onlyI @2Gs/s or 250Ms/s
						if ((DSPTYPE  & 0x20)==0x0):
							# points are already averaged in the PS part
							# format : 16int
							Np = NpCont
							adcdataI[V]=np.concatenate((adcdataI[V], np.right_shift(rep[iStart:iStart+Np],4)*0.3838e-3))

						# mixer ON : I and Q present
						elif ((DSPTYPE  & 0x20)==0x20):
							Np = NpCont
							adcdataI[V]=np.concatenate((adcdataI[V],np.right_shift(rep[iStart:Np:2],4)*0.3838e-3))
							adcdataQ[V]=np.concatenate((adcdataQ[V], np.right_shift(rep[iStart+1:Np:2],4)*0.3838e-3))


					i = iStart+Np # index of the new data block, new header

				# print("********************************************************************")
				# print(len(rep),"Pts treated in ",time.perf_counter()-tstart,"seconds")
				# print("********************************************************************")

				#reshaping results

				points_rec = 0
				points_expected = 0

				for index in range(8):

					if len(adcdataI[index]) > 0:

						points_rec += adcdataI[index].size

					points_expected += int(n_rep * np.sum(length_vec[index],dtype=int))

				if points_rec == points_expected:

					getting_valid_dataset = False
					if self.debug_mode:
						print(adcdataI)

					adcdataI=[np.array(adcdataI[v]).reshape(n_rep,np.sum(length_vec[v],dtype=int)) for v in range(8)]

					I,Q = adcdataI,adcdataQ

				else:

					log.error('Data curruption: rfSoC did not send all data points({}/'.format(points_rec)+str(points_expected)+').')

		else:

			log.error('rfSoC: Instrument mode not recognized.')

		return I,Q

 
	def dump_raw_readout_pulse(self):
		'''
		Legacy function, may not work with the current driver.
		This function dumps raw data to drive to avoid RAM clogging.
		'''
		self.reset_output_data()
		mode = self.acquisition_mode()
		n_rep = self.n_points_total
		length_vec = self.length_vec
		ch_vec = self.ADC_events
		N_adc_events = len(ch_vec)
		n_pulses = len(length_vec[0])
		location = self.raw_dump_location

		ch_active = self.ADC_ch_active

		if mode == 'INT':
			'''
				Get data
			'''

			count_meas = 0
			empty_packet_count = 0
			run_num = 0

			self.write("SEQ:START")
			time.sleep(0.1)

			getting_valid_dataset = True

			if self.display_IQ_progress:

				self.display_IQ_progress_bar = IntProgress(min=0, max=n_rep) # instantiate the bar
				display(self.display_IQ_progress_bar) # display the bar

			while getting_valid_dataset:

				while (count_meas//(16*N_adc_events))<n_rep:

					a = datetime.datetime.now()

					try:
						r = self.visa_handle.query_binary_values('OUTPUT:DATA?', datatype="h",
							is_big_endian=False)
					except:
						r=[]
					# lenr = len(r)
					# if len(r)>1:
					#     len_data_all +=len(r)
					# print(r[0], len_data_all, N_adc_events*16)

					b = datetime.datetime.now()
					# print('\nget data: ',b-a)

					if r == 'ERR':

						log.error('rfSoC: Instrument returned ERR!')

						# reset measurement
						count_meas = 0
						empty_packet_count = 0
						run_num = 0
						self.write("SEQ:STOP")
						time.sleep(2)
						while True:
							try:
								junk = self.visa_handle.quejunky_binajunky_values('OUTPUT:DATA?', datatype="h",
									is_big_endian=False)
							except:
								junk=[]

							# print(junk)
							time.sleep(0.1)
							if junk == [3338] or junk == [2573] or junk == []:
								break
						junk = []
						nb_try +=1
						self.write("SEQ:START")
						time.sleep(0.1)

						continue

					elif len(r)>1:

						a = datetime.datetime.now()
						empty_packet_count = 0
						r_size = len(r)
						pk.dump(r, open(location+"/raw_"+str(run_num)+".pkl","wb"))
						count_meas += r_size
						if self.display_IQ_progress:
							self.display_IQ_progress_bar.value = count_meas//(16*N_adc_events)
						run_num += 1
						b = datetime.datetime.now()
						# print('end storing: ',b-a)

						time.sleep(0.01)


					elif r == [3338] or r == [2573] or r==[]: # new empty packet?

						time.sleep(0.1)
						self.write("SEQ:STOP")
						try:
							junk = self.visa_handle.query_binary_values('OUTPUT:DATA?', datatype="h",
							is_big_endian=False)
						except:
							junk = []
						time.sleep(0.1)
						if junk == [3338] or junk == [2573] or junk ==[]:
							break
						junk = []

						self.write("SEQ:START")
						time.sleep(0.1)
						continue



				if count_meas//(16*N_adc_events) == n_rep:

					getting_valid_dataset = False

				else:

					log.error('Data curruption: rfSoC did not send all data points({}/'.format(count_meas//(16*N_adc_events))+str(n_rep)+').')

					# reset measurement
					data_unsorted = []
					count_meas = 0
					empty_packet_count = 0
					self.write("SEQ:STOP")
					time.sleep(2)
					while True:
						try:
							junk = self.visa_handle.query_binary_values('OUTPUT:DATA?', datatype="h",
							is_big_endian=False)
						except:
							junk = []

						time.sleep(0.1)
						if junk == [3338] or junk == [2573] or junk == []:
							break
					junk = []
					self.write("SEQ:START")
					time.sleep(0.1)

			self.write("SEQ:STOP")

		return run_num


	def transfer_speed(self, block_size=100):

		"""
		Legacy function, may not work with the current driver
		"""

		block_n = int(block_size/10)
		a = datetime.datetime.now()
		for i in bar(range(block_n)):
			try:
				data = self.visa_handle.query_binary_values('OUTPUT:DATATEST?', datatype="h",
							is_big_endian=False)
			except:
				pass

		b = datetime.datetime.now()
		del_t = (b-a).seconds
		speed = round(10*block_n/del_t,2)
		event_rate = round(1000*speed/32,2)
		pulse_length = round(1000/event_rate,2)
		print('Transfer speed: '+str(speed)+' MBps')
		print('Event rate: '+str(event_rate)+' K/s')
		print('Minimum size of one ADC pulse: '+str(pulse_length)+' us per active channel')


	# def ask_raw(self, cmd: str) -> str:
	# 		"""
	# 		Legacy function, may not work with the current driver

	# 		Overwriting the ask_ray qcodes native function to query binary
	# 		Low-level interface to ``visa_handle.ask``.
	# 		Args:
	# 			cmd: The command to send to the instrument.
	# 		Returns:
	# 			str: The instrument's response.
	# 		"""
	# 		with DelayedKeyboardInterrupt():
	# 			keep_trying = True
	# 			count = 0
	# 			while keep_trying:
	# 				count += 1
	# 				self.visa_log.debug(f"Querying: {cmd}")
	# 				try:
	# 					response = self.visa_handle.query_binary_values(cmd, datatype="h", is_big_endian=False)
	# 					self.visa_log.debug(f"Response: {response}")
	# 					if len(response) > 1:
	# 						i = 0
	# 				except:
	# 					try:        # try to read the data as a single point of data in case buffer is empty
	# 						response = self.visa_handle.query_binary_values(cmd, datatype="h", is_big_endian=False, data_points=1, header_fmt='ieee', expect_termination=False)
	# 					except:
	# 						response = 'ERR'
	# 				if response != 'ERR' and response != [3338]:
	# 					keep_trying = False
	# 				if count>10:
	# 					keep_trying = False

	# 		return response


	def int_0(self):
		return 0

	def int_8(self):
		return 8

	def get_idn(self):
		return {'vendor': 'Quantum Coherence Team', 'model': 'rfSoC gen01',
				'serial': 'NA', 'firmware': None}


	def time_conversion(self, t):

		"""
		This function ensure that any time is compatible with the FPGA_clock speed. 
		I.E. is a multiple of 4ns.
		It should be used to define any time to ensure phase stability. 

		Parameter: 
		- t -- (float) time to be converted 
		- t_conv -- (float) time compatible with the RFSoC 
		"""

		nb_clock_per_us = self.FPGA_clock*1e-6
		t_conv = np.round(t*nb_clock_per_us)/nb_clock_per_us
		return t_conv 













	'''
	
	Functions to help with debugging

	'''

	def seq_DAC_ADC_int_to_bin(self,cmd_int):

		cmd_bin = bin(cmd_int)[2:].zfill(32)

		print('ADC command ',cmd_bin[0:8])
		print()

		for i in range(8):

			print('DAC ' + str(i+1) + ' command ',cmd_bin[8+3*i:8+3+3*i])



	def show_sequence_commands(self):

		for i in range(int(len(self.global_sequence.astype(int))/2)):

			print(self.global_sequence_info[i])
			print(self.global_sequence.astype(int)[2*i],self.global_sequence.astype(int)[2*i+1])
			print()