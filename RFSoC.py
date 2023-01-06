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

        self.pulses = pd.DataFrame()
        self.ADC_ch_active = np.zeros(8)
        self.length_vec = [[],[],[],[],[],[],[],[]]
        self.ch_vec = []

        self.display_sequence = True
        self.display_IQ_progress = True
        self.debug_mode = False
        self.debug_mode_plot_waveforms = False
        self.debug_mode_waveform_string = False
        self.loop_time = False 

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



        log.info('Started sequence processing'+'  \n')

        pulses_raw_df = self.pulses

        pulses_raw_df['index'] = range(1, len(pulses_raw_df) + 1)
        pulses_raw_df = pulses_raw_df.set_index('index')

        # --- Check pulse labeling 

        if len(set(pulses_raw_df['label'])) < len(pulses_raw_df['label']):

            log.error('Duplicate Labels: Labels need to be unique for consistent identification of pulse hierarchy.')

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

                    if pulses_raw_df.loc[row['parent']]['parent'] == None:

                        pulses_raw_df.loc[index,'start'] = pulses_raw_df.loc[index,'start'] +  pulses_raw_df.loc[row['parent']]['start'] + pulses_raw_df.loc[row['parent']]['length']
                        pulses_raw_df.loc[index,'parent'] = None

            resolve_hierarchy = False
            for val in pulses_raw_df['parent']:
                if val != None:
                    resolve_hierarchy = True

        if self.debug_mode:

            print('Hierarchy resolution...')
            display(pulses_raw_df)

        # --- Initialisation 

        pulses_df = pd.DataFrame()
        time_ADC = [0,0,0,0,0,0,0,0]
        time_DAC = [0,0,0,0,0,0,0,0]
        length_vec = [[],[],[],[],[],[],[],[]]
        ch_vec = []
        # color of the displayed pulses with n_rep=1
        wait_color = int("D3D3D3", 16)
        DAC_color = int("306cc7", 16)
        ADC_color = int("db500b", 16)
        # color of the displayed pulses with n_rep>1
        wait_color_rep = int("808080", 16)
        DAC_color_rep = int("234e90", 16)
        ADC_color_rep = int("913608", 16)

        wait_count = 0 # counter used to label the wait 
        termination_time = 0 # keep track of the latest event
        ch_demod = None # either or not the ADC LUT will be used during the sequence



        # --- Beginning of the treatment (ADC)

        tmp_df = pulses_raw_df.loc[pulses_raw_df['module'] == 'ADC']


        for index, row in tmp_df.iterrows():

            rep_nb = int(row['repetitions']) # check for repetition
            start = row['start']
            length = row['length']


            # if n_rep>1 we need to know the dead_time, the waiting time between two consecutive pulses 
            if rep_nb>1: 
                dead_time = self.time_conversion(self.dead_time)
            else:
                dead_time = 0

            # create the start and stop of the pulses 
            start_vec = start + np.arange(rep_nb) * (length + dead_time)
            stop_vec = start + np.arange(1, rep_nb + 1) * length + np.arange(rep_nb) * dead_time


            for k in range(rep_nb):

                # check if the event is a wait of a pulses (here wait)
                if start_vec[k] > time_ADC[int(row['channel'])-1]:

                    # in case n_rep>1 it only compute the two first iterations in order to take the dead-time
                    if k < 2:

                        label = 'wait ' + str(wait_count)
                        start = time_ADC[int(row['channel'])-1]
                        stop = start_vec[k]
                        time = (start_vec[k] - time_ADC[int(row['channel'])-1])
                        module = row['module']
                        Channel = 'ADC ch' + str(int(row['channel']))
                        mode = 'wait'
                        param = None
                        ch_num = row['channel']
                        LUT = False
                        ch_demod = row['demodulation_channels']

                        if k == 0:
                            rep_nb_wait = 1
                            color = '#{0:06X}'.format(wait_color)
                        else:
                            rep_nb_wait = rep_nb
                            color = '#{0:06X}'.format(wait_color_rep)

                        pulses_df = pd.concat([pulses_df, pd.DataFrame.from_records([dict(label=label, start=start, stop=stop, time=time, module=module,
                                        Channel=Channel, mode=mode, color=str(color), param=param, ch_num=ch_num,
                                        rep_nb=rep_nb_wait, LUT=LUT, ch_demod=ch_demod)])], ignore_index=True)

                    else:
                        # if k > 2 it only updates the last event time
                        stop = start_vec[k]
                        break
                    

                    wait_count +=1

                # if case n_rep>1 it only compute the first iteration for computation efficiency
                if k < 1:
                    label = index
                    start = start_vec[k]
                    stop = stop_vec[k]
                    time = length
                    module = row['module']
                    Channel = 'ADC ch' + str(int(row['channel']))
                    ch_num = row['channel']

                    rep_nb_sig = rep_nb

                    if rep_nb>1:
                        color = '#{0:06X}'.format(ADC_color_rep)
                    else:
                        color = '#{0:06X}'.format(ADC_color)

                    # fill the parameters only used for the 'IQ_table' mode 
                    LUT = not(np.isnan(row['LUT']))
                    ch_demod = row['demodulation_channels']
                    param = row['param']
                    mode = row['mode']
                    start_pointer = 0 # shift in the starting pointer not implemented so far


                    pulses_df = pd.concat([pulses_df, pd.DataFrame.from_records([dict(label=label, start=start, stop=stop, time=time, module=module,
                                    Channel=Channel, mode=mode, color=str(color), param=param, ch_num=ch_num,start_pointer=start_pointer,
                                    rep_nb=rep_nb_sig, LUT=LUT, ch_demod=ch_demod)])],
                                    ignore_index=True)



                # --- Update of the acquisition parameter used to format the data in get_readout_pulse()
                    nb_points = int(time*1e-6*self.sampling_rate)
                    for chd in ch_demod:
                        ch_vec.append(chd - 1)
                        length_vec[chd-1].append(nb_points)


                else:
                    # if k>1 it only updates the last event time
                    stop = stop_vec[k]

                time_ADC[int(row['channel'])-1] = stop

                if stop>termination_time:
                # update of the termination time used to compute if a wait will be added before the next pulse
                    termination_time = stop



        # --- Beginning of the treatment (DAC)

        tmp_df = pulses_raw_df.loc[pulses_raw_df['module'] == 'DAC']
        for index, row in tmp_df.iterrows():

            rep_nb = int(row['repetitions'])
            start = row['start']
            length = row['length']

            start_vec = start + np.arange(rep_nb) * (length)
            stop_vec = start + np.arange(1, rep_nb + 1) * length

            for k in range(rep_nb):

                if start_vec[k] > time_DAC[int(row['channel'])-1]:

                    if k < 2:
                        label = 'wait ' + str(wait_count)
                        start = time_DAC[int(row['channel'])-1]
                        stop = start_vec[k]
                        time = (start_vec[k] - time_DAC[int(row['channel'])-1])
                        module = row['module']
                        Channel = 'DAC ch' + str(int(row['channel']))
                        mode = 'wait'
                        param = None
                        ch_num = row['channel']
                        LUT = False

                        # should not repeat the first waiting time before the loop
                        if k == 0:
                            rep_nb_wait = 1
                            color = '#{0:06X}'.format(wait_color)
                        else:
                            rep_nb_wait = rep_nb
                            color = '#{0:06X}'.format(wait_color_rep)

                        pulses_df = pd.concat([pulses_df, pd.DataFrame.from_records([dict(label=label, start=start, stop=stop, time=time, module=module,
                                        Channel=Channel, mode=mode, color=str(color), param=param, ch_num=ch_num,
                                        rep_nb=rep_nb_wait, LUT=LUT)])],
                                        ignore_index=True)

                    else:
                        stop = start_vec[k]
                        break

                    wait_count +=1


                if k < 1:

                    label = index
                    start = start_vec[k]
                    stop = stop_vec[k]
                    time = length
                    module = row['module']
                    Channel = 'DAC ch' + str(int(row['channel']))
                    mode = row['mode']
                    param = row['param']
                    ch_num = row['channel']
                    LUT = not(np.isnan(row['LUT']))
                    start_pointer = row['starting_pointer']

                    rep_nb_sig = rep_nb

                    if rep_nb>1:
                        color = '#{0:06X}'.format(DAC_color_rep)
                    else:
                        color = '#{0:06X}'.format(DAC_color)

                    pulses_df = pd.concat([pulses_df, pd.DataFrame.from_records([dict(label=label, start=start, stop=stop, time=time, module=module,
                                    Channel=Channel, mode=mode, color=str(color), param=param, ch_num=ch_num,
                                    rep_nb=rep_nb_sig, LUT=LUT, start_pointer=start_pointer)])],
                                    ignore_index=True)

                else:
                    stop = stop_vec[k]


                time_DAC[int(row['channel'])-1] = stop

                if stop>termination_time:

                    termination_time = stop




        # --- Pulse sequences treatment solving inconsistency regarding stop time when n_rep>1 and more than one pulses is used. 

        pulses_fin = pd.DataFrame()
        channel_list = set(list(pulses_df['Channel']))

        for idx, ch in enumerate(channel_list):
            pulses_ch =  pulses_df.loc[pulses_df['Channel']==ch]
            pulses_ch = pulses_ch.sort_values('start')
            for i in range(len(pulses_ch) - 1):
                pulses_ch.iloc[i, pulses_ch.columns.get_loc('stop')] = pulses_ch.iloc[i + 1]['start']
                pulses_ch.iloc[i, pulses_ch.columns.get_loc('time')] = pulses_ch.iloc[i]['stop'] - pulses_ch.iloc[i]['start']
            pulses_fin = pd.concat([pulses_fin, pulses_ch], ignore_index=True)

        # --- Sequence plotting 

        if self.display_sequence:

            fig = make_subplots(rows=len(channel_list), cols=1, shared_xaxes=True)            

            for idx, ch in enumerate(channel_list):
                pulses_loop =  pulses_fin.loc[pulses_fin['Channel']==ch]
                pulses_loop = pulses_loop.sort_values('start')

                fig.add_trace(go.Bar(x=pulses_loop.time, y=pulses_loop.Channel, orientation='h', text=pulses_loop.label, marker=dict(color=pulses_loop.color), 
                                        customdata=np.dstack((pulses_loop.label, pulses_loop.start, pulses_loop.stop))[0],
                                        hovertemplate='label: %{customdata[0]: s} <br>start: %{customdata[1]:.3f} <br>stop: %{customdata[2]:.3f}'), idx + 1, 1)

            fig.update_layout(showlegend=False)
            fig.show()


        # store the length of the pulses and the channel used, it is use for data shaping in get_readout_pulse()
        self.length_vec = length_vec
        self.ch_vec = ch_vec


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

            LUT_df = pulses_df.loc[pulses_df['module'] == 'DAC']

            if self.debug_mode:

                print('\n------------------------------- LUT DEBUGGING DAC ----------------------------------\n')

            event_time_list = list(dict.fromkeys(LUT_df['start']))
            event_time_list.sort()


            for event_time in event_time_list:
                event_df = LUT_df.loc[LUT_df['start'] == event_time]

                for index, row in event_df.iterrows():

                    ch_num = int(row['ch_num'])

                    if not(np.isnan(row['start_pointer'])):
                        if row['LUT']:
                            pulse_addr = round((len(DAC_pulses_array[ch_num-1]) + row['start_pointer']*1e-6*self.sampling_rate)/11)
                            last_pointer = pulse_addr

                        else:
                            pulse_addr =  last_pointer + round((row['start_pointer']*1e-6*self.sampling_rate)/8)
                            
                            
                        DAC_pulses_pointer[ch_num-1].append(pulse_addr)

                        if self.debug_mode:
                            print('For pulse ' + row['label'])
                            print('Saved pointer: ' + str(pulse_addr))

                    if row['LUT']:

                        SCPI_command = self.pulse_gen_SCPI(row['mode'],row['param'],row['time'],ch_num, mode='DAC')
                        # adding pulse to waveform
                        DAC_pulses_array[ch_num-1] = np.append(DAC_pulses_array[ch_num-1],SCPI_command)


            if self.debug_mode:

                print('Pointer tab:')
                print(DAC_pulses_pointer)


            for i in range(8):

                if len(DAC_pulses_array[i])>0:

                    self.write('DAC:DATA:CH{}:CLEAR'.format(str(i+1)))

                    if self.debug_mode and self.debug_mode_plot_waveforms:

                        fig = plt.figure(figsize=(8,5))
                        plt.plot(range(len(DAC_pulses_array[i])),DAC_pulses_array[i])
                        plt.grid()
                        plt.legend(fontsize = 14)
                        plt.show()

                    DAC_SCPI_cmd = 'DAC:DATA:CH' + str(i+1) + ' 0,' + ','.join((DAC_pulses_array[i].astype(int)).astype(str)) 

                    # if self.debug_mode and self.debug_mode_waveform_string:
                    if self.debug_mode_waveform_string:

                        print('DAC sequence for CH '+str(i+1)+': ',DAC_SCPI_cmd)

                    log.info('Writing waveform for CH'+str(i+1)+'  \n')
                    self.write(DAC_SCPI_cmd)

            log.info('Waveform processing complete' + '\n')

        if ADC:

            LUT_df = pulses_df.loc[pulses_df['module'] == 'ADC']

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

                print('\n------------------------------- LUT DEBUGGING ADC ----------------------------------\n')
                display(LUT_df)


            event_time_list = list(dict.fromkeys(LUT_df['start']))
            event_time_list.sort()

            for event_time in event_time_list:
                event_df = LUT_df.loc[LUT_df['start'] == event_time]

                for index, row in event_df.iterrows():

                    if not(np.isnan(row['start_pointer'])):

                        ch_num = int(row['ch_num'])
                        ch_demod = row['ch_demod']

                        for idx, chd in enumerate(ch_demod):
                            if row['LUT']:

                                if self.debug_mode:
                                    print('Pulse not found in ADC memory')
                                    print('Starting pointer set at: %i'%round(len(ADC_I_pulses_array[chd-1])/8))

                                pulse_addr_start = len(ADC_I_pulses_array[chd-1])//8
                                pulse_addr_loop = pulse_addr_start


                            else:
                                pulse_addr_start = ADC_pulses_pointer[chd-1][0][-1]
                                pulse_addr_stop = ADC_pulses_pointer[chd-1][2][-1]
                                pulse_addr_loop = ADC_pulses_pointer[chd-1][1][-1]


                            if row['LUT']:
                                if self.debug_mode:
                                    print('The pulses saved in the LUT are:')
                                    print(row['mode'], row['LUT'], row['label'])

                                if row['mode']=='sin':
                                    if self.debug_mode:
                                        print(row['param'])

                                    param_I = row['param'][idx]
                                    freq_demod = param_I['freq']

                                    if 1e3%freq_demod!=0:
                                        log.error('Demodulation frequency is not a multiple of the sampling_rate. \
                                                   As a result, you are not using the loop pointer and the memory \
                                                   usage of the ADC LUT is suboptimal.')
                                        time_vec = row['time']


                                    else:
                                        period_demod = int(1e3/freq_demod) # in ns
                                        period_loop = int(1e9/self.sampling_rate*8) # in ns
                                        time_vec = np.lcm(period_demod, period_loop)*1e-3 # in µs

                                        if self.debug_mode:
                                            print('Demodulation period/ Loop period/ Time vec: ')
                                            print(period_demod, period_loop, time_vec)


                                    param_Q = param_I

                                    param_Q = {**param_I, 'phase_offset': param_I['phase_offset'] + np.pi/2}


                                    SCPI_command_I = self.pulse_gen_SCPI(row['mode'], param_I, time_vec, chd, mode='ADC')
                                    SCPI_command_Q = self.pulse_gen_SCPI(row['mode'], param_Q, time_vec, chd, mode='ADC')

                                    if self.debug_mode:
                                        print(SCPI_command_I)

                                    ADC_I_pulses_array[chd-1] = np.append(ADC_I_pulses_array[chd-1], SCPI_command_I)
                                    ADC_Q_pulses_array[chd-1] = np.append(ADC_Q_pulses_array[chd-1], SCPI_command_Q)



                                else:
                                    log.error('Other functions than sin are not supported for demodulation for now')



                            pulse_addr_stop = len(ADC_I_pulses_array[chd-1])//8 - 1


                            ADC_pulses_pointer[chd-1][0].append(pulse_addr_start)
                            ADC_pulses_pointer[chd-1][2].append(pulse_addr_stop)
                            ADC_pulses_pointer[chd-1][1].append(pulse_addr_loop)

                            if self.debug_mode:
                                print('For pulse ' + row['label'])
                                print('Saved start/loop/stop pointer: ', str(pulse_addr_start), str(pulse_addr_loop), str(pulse_addr_stop))


            if self.debug_mode:
                print('Pointer tab:')
                print(ADC_pulses_pointer)


            # reset memory :

            for i in range(4):

                empty_list = str([0]*16384)[1:-1]
                empty_list.replace(' ', '')

                ADC_SCPI_cmd_I = 'ADC:I' + str(i+1) +' 0,' + empty_list
                ADC_SCPI_cmd_Q = 'ADC:Q' + str(i+1) +' 0,' + empty_list

                self.write(ADC_SCPI_cmd_I)
                self.write(ADC_SCPI_cmd_Q)


            for i in range(4):

                if len(ADC_I_pulses_array[i])>0:


                    # not implemented for now
                    # self.write('ADC:DATA:CH{}:CLEAR'.format(str(i+1)))

                    if self.debug_mode and self.debug_mode_plot_waveforms:

                        fig, ax = plt.subplots(1, 1, figsize=(8,5))

                        ax[0].plot(range(len(ADC_I_pulses_array[i])),ADC_I_pulses_array[i])
                        ax[1].plot(range(len(ADC_Q_pulses_array[i])),ADC_Q_pulses_array[i])
                        fig.show()



                    ADC_SCPI_cmd_I = 'ADC:I' + str(i+1) +' 0,' + ','.join((ADC_I_pulses_array[i].astype(int)).astype(str))
                    ADC_SCPI_cmd_Q = 'ADC:Q' + str(i+1) +' 0,' + ','.join((ADC_Q_pulses_array[i].astype(int)).astype(str))

                    # if self.debug_mode and self.debug_mode_waveform_string:
                    if self.debug_mode_waveform_string:

                        print('ADC sequence for I'+str(i+1)+':', ADC_SCPI_cmd_I)
                        print('ADC sequence for Q'+str(i+1)+':', ADC_SCPI_cmd_Q)

                    log.info('Writing waveform for I'+str(i+1)+'  \n')
                    self.write(ADC_SCPI_cmd_I)

                    log.info('Writing waveform for Q'+str(i+1)+'  \n')
                    self.write(ADC_SCPI_cmd_Q)

            log.info('Waveform processing complete' + '\n')

        return DAC_pulses_pointer, ADC_pulses_pointer


    def process_sequencing_IQ_table(self):

        """
        This function take the stored pulses sequence and creat the SCPI command 
        that is sent to the RFSoC sequencer. 
        It also fill the ADC/DAC LUT via LUT_and_adress_filling()
        """

        # --- Charging the pulse sequence and the event list

        pulses_df = self.pulses_sequence()
        event_time_list = list(dict.fromkeys(pulses_df['start']))
        event_time_list.sort()
        # temrination time used to close the sequence 
        termination_time = np.max((np.array(list(pulses_df['stop']))))

        # --- Initialisation of the parameters


        # n_points_total is the total number of acquisitions
        self.n_points_total = 0
        # array containing the instructions played by the sequencer
        global_sequence = np.array([])
        # load the pointer
        DAC_pulses_pointer, ADC_pulses_pointer = self.LUT_and_adress_filling(pulses_df)


        # lists of boolean storing the ADC/DAC state
        ADC_state = [0,0,0,0,0,0,0,0]
        DAC_state = [0,0,0,0,0,0,0,0]
        ADC_state_prev = [0,0,0,0,0,0,0,0]
        DAC_state_prev = [0,0,0,0,0,0,0,0]
        self.ADC_ch_active = np.array([0,0,0,0,0,0,0,0])

        # lists storing the ADC/DAC pointer index 
        pointer_dac = [0,0,0,0,0,0,0,0]
        pointer_adc = [0,0,0,0]
        # boolean checking if the repetition started
        rep_started = False
        # array used to store with physical ADC is send to which ADC LUT
        mux_state = np.zeros((4, 4), dtype='int')



        # ---  If a repetitions is used it check the number of pulses per repetitions
        nb_pulses_dac = np.zeros(8)
        nb_pulses_adc = np.zeros(8)
        pulses_counter_dac = np.ones(8)
        pulses_counter_adc = np.ones(8)

        pulses_rep_all = pulses_df.loc[(pulses_df['rep_nb']>1) & (pulses_df['mode']!='wait')]
        pulses_rep_dac = pulses_rep_all.loc[(pulses_rep_all['module']=='DAC')]
        pulses_rep_adc = pulses_rep_all.loc[(pulses_rep_all['module']=='ADC')]

        # count the pulses number per dac 

        if len(pulses_rep_dac)>0:
            for index, row in pulses_rep_dac.iterrows():
                ch_num = int(row['ch_num'])
                nb_pulses_dac[ch_num - 1] +=1 

        # same for the adc 

        if len(pulses_rep_adc)>0:
            for index, row in pulses_rep_adc.iterrows():
                ch_demod = row['ch_demod']
                for chd in ch_demod: 
                    nb_pulses_adc[chd - 1] +=1
        
        # --- Start of the sequence filling 

        if self.debug_mode:
            print('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*')
            print('*-*-*-*-*-*-*-* Beggining of sequence *-*-*-*-*-*-*-*-*-*-*')
            print('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*')
            print('Events detected at: ',event_time_list)
            print('Found the following pulses:')
            display(pulses_df.sort_values('start'))
            rep_nb = np.max(pulses_df['rep_nb'])
            print('They will all be played a first time and then the pulses with rep_nb=%i will be played %i times' %(rep_nb, rep_nb))

        for event_time in event_time_list:
            if event_time>0:

                # add a waiting time any time it is needed 
                global_sequence = np.append(global_sequence,1)
                wait_time = int(round((event_time-event_time_prev)*250) - 1)
                global_sequence = np.append(global_sequence, wait_time)

                if self.debug_mode:
                    print('adding wait till this event')
                    print(1,int(round((event_time-event_time_prev)*250)))


            event_time_prev = event_time

            # take all events occuring at a given time
            tmp_df = pulses_df.loc[pulses_df['start'] == event_time]
            tmp_df = tmp_df.sort_values(by='module', ascending=False)

            for index, row in tmp_df.iterrows():

                if self.debug_mode:
                    print('_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _')
                    print(event_time,row['mode'], row['label'])
                    print('_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _')

                ch_num = int(row['ch_num'])


                if row['module'] == 'DAC':

                    if row['rep_nb']!=1 and not(rep_started):
                        # add the command indicating the start of the loop if needed
                        global_sequence = np.append(global_sequence, 257)
                        global_sequence = np.append(global_sequence, row['rep_nb'] - 1)
                        rep_started = not(rep_started) # triggers the loop start

                        if self.debug_mode:
                            print('- - - - - - - LOOP START - - - - - - - - -')


                    if not(np.isnan(row['start_pointer'])):
                        # if there is a start pointer it must be a pulse 

                        # set the starting pointer to the correct address 
                        global_sequence = np.append(global_sequence,4096+ch_num)
                        global_sequence = np.append(global_sequence,DAC_pulses_pointer[ch_num-1][pointer_dac[ch_num - 1]])

                        if self.debug_mode:
                            print('adding sequencer command to point to address of this pulse')
                            print(4096+ch_num,DAC_pulses_pointer[ch_num-1][pointer_dac[ch_num - 1]])


                        # the corresponding DAC is set to ON
                        DAC_state[7-(ch_num-1)] = 1
                        if self.debug_mode:
                            print('DAC state is :', DAC_state)


                        # the pointer index of the given channel is updated or not depending on the situation
                        if not(rep_started):
                            pointer_dac[ch_num - 1] +=1
                            if self.debug_mode:
                                print('Pointer shifted to %i'%pointer_dac[ch_num - 1])
                        else: 
                            while pulses_counter_dac[ch_num - 1]<nb_pulses_dac[ch_num - 1]:
                                pointer_dac[ch_num - 1] +=1
                                pulses_counter_dac[ch_num - 1] +=1 

                                if self.debug_mode:
                                    print('Pointer shifted to %i'%pointer_dac[ch_num - 1])
                                    print('Pulse counter shifted to %i'%pulses_counter_dac[ch_num - 1])


                    if row['mode'] == 'wait':
                        # we switch off the correponding DAC, the waiting command is set at the next iteration
                        DAC_state[7-(ch_num-1)] = 0
                        if self.debug_mode:
                            print('DAC state is :', DAC_state)


                if row['module'] == 'ADC':

                    if self.debug_mode:
                        print(row['mode'])


                    # initialisation of the ADC to ADC LUT connection for this step
                    mux_step = ['000', '000', '000', '000']

                    ch_demod = row['ch_demod']

                    if row['rep_nb']!=1 and not(rep_started):

                        # add the command indicating the start of the loop if needed
                        global_sequence = np.append(global_sequence, 257)
                        global_sequence = np.append(global_sequence, row['rep_nb'] - 1)

                        rep_started = not(rep_started)

                        if self.debug_mode:
                            print('- - - - - - - LOOP START - - - - - - - - -')

                    if self.debug_mode:
                        print(mux_state)


                    if row['mode'] != 'wait':

                        # --- Update the ADC to ADC LUT connection

                        for k in ch_demod:
                            # check that the given ADC is not ON, if not no change in the routing
                            if ch_num == 1 and ADC_state[7-(ch_num-1)] == 0:
                                # check that the given ADC LUT is not already used  
                                if k ==1 and np.sum(mux_state[:, 0])==0:
                                    mux_step[0] = '000' # binary corresponding to the routing
                                    mux_state[ch_num-1, 0] = 1 # indicate that the ADC LUT is used 
                                elif k ==2 and np.sum(mux_state[:, 1])==0:
                                    mux_step[1] = '001'
                                    mux_state[ch_num-1, 1] = 1
                                elif k ==3 and np.sum(mux_state[:, 2])==0:
                                    mux_step[2] ='010'
                                    mux_state[ch_num-1, 2] = 1
                                elif k ==4 and np.sum(mux_state[:, 3])==0:
                                    mux_step[3] ='010'
                                    mux_state[ch_num-1, 3] = 1
                                else: log.error('Incompatible mixing tables for ch%i'%ch_num)
                            elif ch_num == 2 and ADC_state[7-(ch_num-1)] == 0:
                                if k ==1 and np.sum(mux_state[:, 0])==0:
                                    mux_step[0]='001'
                                    mux_state[ch_num-1, 0] = 1
                                elif k ==2 and np.sum(mux_state[:, 1])==0:
                                    mux_step[1] = '000'
                                    mux_state[ch_num-1, 1] = 1
                                elif k ==3 and np.sum(mux_state[:, 2])==0:
                                    mux_step[2] ='011'
                                    mux_state[ch_num-1, 2] = 1
                                elif k ==4 and np.sum(mux_state[:, 3])==0:
                                    mux_step[3] ='011'
                                    mux_state[ch_num-1, 3] = 1
                                else: log.error('Incompatible mixing tables ch%i'%ch_num)
                            elif ch_num ==3  and ADC_state[7-(ch_num-1)] == 0:
                                if k ==1 : log.error('Incompatible mixing tables')
                                elif k ==2 : log.error('Incompatible mixing tables')
                                elif k ==3 and np.sum(mux_state[:, 2])==0:
                                 mux_step[2] = '000'
                                 mux_state[ch_num-1, 2] = 1
                                elif k ==4 and np.sum(mux_state[:, 3]):
                                 mux_step[3] ='001'
                                 mux_state[ch_num-1, 3] = 1
                                else: log.error('Incompatible mixing tables ch%i'%ch_num)
                            elif ch_num ==4 and ADC_state[7-(ch_num-1)] == 0:
                                if k ==1 : log.error('Incompatible mixing tables')
                                elif k ==2 : log.error('Incompatible mixing tables')
                                elif k ==3 and np.sum(mux_state[:, 2])==0:
                                 mux_step[2] ='001'
                                 mux_state[ch_num-1, 2] = 1
                                elif k ==4 and np.sum(mux_state[:, 3])==0:
                                 mux_step[3] = '000'
                                 mux_state[ch_num-1, 3] = 1
                                else: log.error('Incompatible mixing tables ch%i'%ch_num)

                            else:
                                log.error('Cannot use other channels than 1, 2, 3 or 4 for the IQ_table mode')

                            if self.debug_mode:
                                print('The Mux state:')
                                print(mux_state)
                                print('The Mux table:')
                                print(mux_step)


                        # --- Sequence filling for ADC type command

                        for idx, chd in enumerate(ch_demod):


                            #  --- creation of the bit string storing the mixer state and the start pointer 
                            bin_cmd = mux_step[chd - 1] # mixer routing
                            bin_cmd += '1' # set the mixer to ON
                            bin_cmd += '00000000000000'# used bit

                            # convert the pointer postition in binary 
                            bin_start = bin(ADC_pulses_pointer[chd - 1][0][pointer_adc[chd - 1]])[2:]

                            # fill the unused bits 
                            len_bit_start = len(bin_start)
                            bin_start_add = [0] * (14 - len_bit_start)
                            bin_start_add = str(bin_start_add)[1:-1].replace(', ', '')

                            # add the two bit staring
                            bin_cmd += bin_start_add
                            bin_cmd += bin_start

                            # add the mixer state and the start pointer for the given demod channel
                            global_sequence = np.append(global_sequence, 4128+(chd - 1)*2)
                            global_sequence = np.append(global_sequence, int(bin_cmd,2))

                            if self.debug_mode:
                                print('ADC pointer/ demod ch/ pointer ADC:')
                                print(ADC_pulses_pointer, chd, pointer_adc)
                                print('adding sequencer the mixer state and the starting_pointer:')
                                print('binary command:', bin_cmd)
                                print(4128+2*(chd - 1), int(bin_cmd,2))



                            # --- creation of the bit string storing the loop and stop pointer

                            # convert the pointer postitions in binary 
                            bin_loop = bin(ADC_pulses_pointer[chd - 1][1][pointer_adc[chd - 1]])[2:]
                            bin_stop = bin(ADC_pulses_pointer[chd - 1][2][pointer_adc[chd - 1]])[2:]


                            # fill the unused bits
                            len_bit_loop = len(bin_loop)
                            bin_loop_add = [0] * (14 - len_bit_loop)
                            bin_loop_add = str(bin_loop_add)[1:-1].replace(', ', '')
                            bin_loop = bin_loop_add + bin_loop

                            len_bit_stop = len(bin_stop)
                            bin_stop_add = [0] * (14 - len_bit_stop)
                            bin_stop_add = str(bin_stop_add)[1:-1].replace(', ', '')
                            bin_stop = bin_stop_add + bin_stop

                            # add the two bit strings 
                            bin_cmd = '00' + bin_loop + '00' + bin_stop

                            # add loop and stop pointer for the given demod channel
                            global_sequence = np.append(global_sequence, 4129+(chd - 1)*2)
                            global_sequence = np.append(global_sequence, int(bin_cmd,2))

                            if self.debug_mode:
                                print('adding sequencer the loop and stopping pointer:')
                                print('binary command:', bin_cmd)
                                print(4129+2*(chd - 1), int(bin_cmd,2))


                            # add the number of points the ADC LUT should take 
                            global_sequence = np.append(global_sequence,4107 + (chd - 1))
                            global_sequence = np.append(global_sequence,int(row['time']*1e-6*self.sampling_rate))


                            if self.debug_mode:
                                print('adding sequencer command to set acq points')
                                print(4106+(chd - 1),int(row['time']*1e-6*self.sampling_rate))


                            # change the ADC LUT to ON 
                            ADC_state[7-(chd-1)] = 1
                            self.ADC_ch_active[chd-1] = 1

                            if self.debug_mode:
                                print('ADC state is :', ADC_state)


                            # the pointer index of the given channel is updated or not depending on the situation
                            if not(rep_started):
                                pointer_adc[chd - 1] +=1
                            else: 
                                while pulses_counter_adc[chd - 1] < nb_pulses_adc[chd- 1]:
                                    pointer_adc[chd - 1] +=1
                                    pulses_counter_adc[chd - 1] +=1 

                        self.n_points_total += row['rep_nb']


                    elif row['mode'] == 'wait':
                        ch_demod = row['ch_demod']
                        for idx, chd in enumerate(ch_demod):
                            ADC_state[7-(chd-1)] = 0
                            # reset the mixer state of the switched off channel
                            mux_state[chd - 1, :] = np.zeros(4)

                            if self.debug_mode:
                                print('ADC state is :', ADC_state)


            # ---  Update of the ADC and DAC state at every time step 

            if DAC_state != DAC_state_prev or ADC_state != ADC_state_prev:

                if self.debug_mode:
                    print('ADC state updated from %s to %s'%(ADC_state_prev, ADC_state))
                    print('DAC state updated from %s to %s'%(DAC_state_prev, DAC_state))


                # update the ADC state 
                bin_adc_cmd = str(ADC_state)[1:-1].replace(', ', '')
                bin_dac_cmd = ''

                # update the DAC state
                for i in range(8):
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
                bin_trig_cmd = bin_adc_cmd + bin_dac_cmd
                global_sequence = np.append(global_sequence,4096)
                global_sequence = np.append(global_sequence,int(bin_trig_cmd,2))

                ADC_state_prev = ADC_state.copy()

                if self.debug_mode:
                    print('Bit string of the DAC state:')
                    print(bin_dac_cmd)
                    print('Bit string of the ADC state')
                    print(bin_adc_cmd)
                    print('The command is: ')
                    print(4096, int(bin_trig_cmd,2))

            else :
                if self.debug_mode:
                    print('No ADC or DAC update at this state')
                    print('ADC previous: %s  / ADC step : %s'%(ADC_state_prev, ADC_state))
                    print('DAC previous: %s  / DAC step : %s'%(DAC_state_prev, DAC_state))



        # --- Add a last waiting time if needing 

        wait_term = int(round((termination_time - event_time)*250))-1
        global_sequence = np.append(global_sequence,1)
        global_sequence = np.append(global_sequence,wait_term)

        if self.debug_mode:
            print('Terminate by wait of : %f' %(wait_term + 1))


        # Switch off all the DAC/ADC 
        global_sequence = np.append(global_sequence,4096)
        global_sequence = np.append(global_sequence,0)

        # Close the loop if it was open
        if rep_started:
            global_sequence = np.append(global_sequence, 513)
            global_sequence = np.append(global_sequence, 0)


        # --- Set the Acquisition mode 
        if self.acquisition_mode() == 'RAW':
            acq_mode = 0
        elif self.acquisition_mode() == 'INT':
            acq_mode = 286331153
        else:
            log.error('Invalid acquisition mode\n')


        # if n_points()>1 we add a global loop 
        if self.n_points() > 1:
            period_sync = int(self.FPGA_clock/self.freq_sync())
            global_sequence_str = 'SEQ 0,1,9,4106,' + str(acq_mode) + ',258,' + str(int(self.n_points()-1)) + ',' + ','.join((global_sequence.astype(int)).astype(str))  + ',514,0,0,0'

        else :
            global_sequence_str = 'SEQ 0,1,9,4106,' + str(acq_mode) + ',' + ','.join((global_sequence.astype(int)).astype(str)) + ',0,0'

        if self.debug_mode:
            print('Sequence programmer command: ',global_sequence_str)

        # --- Send the SCPI command 
        log.info('Writing global sequence' + '\n')
        self.write(global_sequence_str)

        # Update the total acquisition points 
        self.n_points_total *=self.n_points()

        # reset the pointer indexes
        pointer_adc = [0,0,0,0,0,0,0,0]
        pointer_dac = [0,0,0,0,0,0,0,0]


    def pulse_gen_SCPI(self,function,param,duration,ch, mode):

        """
        This function convert the pulses parameters into a 1D array of Volt amplitude in bit scale.
        This is used in LUT_and_adress_filling() to creat the SCPI command. 

        Parameters: 
        function -- (str) function to generate, can be 'sin+sin' for a sum of sin
        'sin' for a sin, 'trigger' for a trigger and 'DC' for a DC channel 
        param --  list containing the function parameters
        duration -- (float) pulse duration 
        ch --  (int) LUT channel
        mode -- (str) 'DAC' or 'ADC'

        return: 
        wavepoints -- 1D array of points that will be stored in the LUT 
        """

        period = 1./self.sampling_rate
        time_vec = np.arange(period,duration*1e-6+period/2,period)
        # stop vector used to separate different pulses in a channel 
        stop_vector = np.array([0,0,0,0,0,0,0,0,0,0,16383])

        if function == 'sin+sin':

            wavepoints1 = ((2**13)*self.DAC_amplitude_calib[ch-1]*param['amp1']*np.sin(-param['phase_offset1'] - 1) + 2*np.pi*param['freq1']*1e6*time_vec)
            wavepoints2 = ((2**13)*self.DAC_amplitude_calib[ch-1]*param['amp2']*np.sin(-param['phase_offset2'] - 1) + 2*np.pi*param['freq2']*1e6*time_vec)

            wavepoints = ((2**13)*self.DAC_amplitude_calib[ch-1]*param['dc_offset'] - 1)+ wavepoints1 + wavepoints2

            # check if some points are above the limit imposed by the 13 bits
            idx_of = np.nonzero(((wavepoints > 8191) | (wavepoints < -8192)))[0]
            if len(np.nonzero(wavepoints[idx_of] != 16381)[0])>0:
                if mode =='DAC':
                    log.error('Error when filling the DAC memory: maximal amplitude is over the resolution') 
                elif mode == 'ADC':
                    log.error('Error when filling the ADC memory: maximal amplitude is over the resolution') 


            if self.debug_mode and self.debug_mode_plot_waveforms:
                print('plot of sinsin mode 1')
                fig = plt.figure(figsize=(8,5))
                plt.plot(time_vec,wavepoints1)
                plt.grid()
                plt.show()
                print('plot of sinsin mode 2')
                fig = plt.figure(figsize=(8,5))
                plt.plot(time_vec,wavepoints2)
                plt.grid()
                plt.show()
                print('plot of sinsin mode total')
                fig = plt.figure(figsize=(8,5))
                plt.plot(time_vec,wavepoints)
                plt.grid()
                plt.show()


            if mode == 'DAC':

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
                wavepoints = np.append(wavepoints,stop_vector)

                



        elif function == 'sin':

            wavepoints = ((2**13)*self.DAC_amplitude_calib[ch-1]*param['dc_offset'] - 1) + ((2**13)*self.DAC_amplitude_calib[ch-1]*param['amp'] - 1)*np.sin(-param['phase_offset'] + 2*np.pi*param['freq']*1e6*time_vec)

            idx_of = np.nonzero(((wavepoints > 8191) | (wavepoints < -8192)))[0]
            if len(np.nonzero(wavepoints[idx_of] != 16381)[0])>0:
                if mode =='DAC':
                    log.error('Error when filling the DAC memory: maximal amplitude is over the resolution') 
                elif mode == 'ADC':
                    log.error('Error when filling the ADC memory: maximal amplitude is over the resolution') 

            if self.debug_mode and self.debug_mode_plot_waveforms:
                print('plot of sin mode')
                fig = plt.figure(figsize=(8,5))
                plt.plot(time_vec*1e9,wavepoints)
                plt.grid()
                plt.legend(fontsize = 14)
                plt.show()

            if mode == 'DAC':

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
                wavepoints = np.append(wavepoints,stop_vector)




        elif function == 'trigger':

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


        elif function == 'DC':

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

            log.error('Wrong waveform mode: ',function)


        if len(wavepoints) > 128000: 
            if mode == 'DAC':
                log.error('Error when filling the DAC memory : to many points, maximal number is 128000 while you are asking for %i'%len(wavepoints))
            if mode == 'ADC':
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


    def reset_PLL(self):

        self.write("DAC:RELAY:ALL 0")
        self.write("PLLINIT")
        time.sleep(5)
        self.write("DAC:RELAY:ALL 1")


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
        ch_vec = self.ch_vec
        # N_adc_events = len(ch_vec) # to be discussed with Arpit
        N_adc_events = len(np.unique(ch_vec))
        if self.debug_mode:
            len_data_all = 0
        # print(N_adc_events, ch_vec)
        #print(length_vec)
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
            num_points = np.frombuffer(np.stack((raw_IQ_data_dump_header.T[1], raw_IQ_data_dump_header.T[2]), axis=1).astype('int16').tobytes(), dtype=np.long)

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

            # I = [((I_all_data*ch_1)[I_all_data*ch_1!=0]-2).reshape(n_rep*ch_active[0],n_pulses).T,
            #      ((I_all_data*ch_2)[I_all_data*ch_2!=0]-2).reshape(n_rep*ch_active[1],n_pulses).T,
            #      ((I_all_data*ch_3)[I_all_data*ch_3!=0]-2).reshape(n_rep*ch_active[2],n_pulses).T,
            #      ((I_all_data*ch_4)[I_all_data*ch_4!=0]-2).reshape(n_rep*ch_active[3],n_pulses).T,
            #      ((I_all_data*ch_5)[I_all_data*ch_5!=0]-2).reshape(n_rep*ch_active[4],n_pulses).T,
            #      ((I_all_data*ch_6)[I_all_data*ch_6!=0]-2).reshape(n_rep*ch_active[5],n_pulses).T,
            #      ((I_all_data*ch_7)[I_all_data*ch_7!=0]-2).reshape(n_rep*ch_active[6],n_pulses).T,
            #      ((I_all_data*ch_8)[I_all_data*ch_8!=0]-2).reshape(n_rep*ch_active[7],n_pulses).T]
            # Q = [((Q_all_data*ch_1)[Q_all_data*ch_1!=0]-2).reshape(n_rep*ch_active[0],n_pulses).T,
            #      ((Q_all_data*ch_2)[Q_all_data*ch_2!=0]-2).reshape(n_rep*ch_active[1],n_pulses).T,
            #      ((Q_all_data*ch_3)[Q_all_data*ch_3!=0]-2).reshape(n_rep*ch_active[2],n_pulses).T,
            #      ((Q_all_data*ch_4)[Q_all_data*ch_4!=0]-2).reshape(n_rep*ch_active[3],n_pulses).T,
            #      ((Q_all_data*ch_5)[Q_all_data*ch_5!=0]-2).reshape(n_rep*ch_active[4],n_pulses).T,
            #      ((Q_all_data*ch_6)[Q_all_data*ch_6!=0]-2).reshape(n_rep*ch_active[5],n_pulses).T,
            #      ((Q_all_data*ch_7)[Q_all_data*ch_7!=0]-2).reshape(n_rep*ch_active[6],n_pulses).T,
            #      ((Q_all_data*ch_8)[Q_all_data*ch_8!=0]-2).reshape(n_rep*ch_active[7],n_pulses).T]



            I = [((I_all_data*ch_1)[I_all_data*ch_1!=0]-2),
                 ((I_all_data*ch_2)[I_all_data*ch_2!=0]-2),
                 ((I_all_data*ch_3)[I_all_data*ch_3!=0]-2),
                 ((I_all_data*ch_4)[I_all_data*ch_4!=0]-2),
                 ((I_all_data*ch_5)[I_all_data*ch_5!=0]-2),
                 ((I_all_data*ch_6)[I_all_data*ch_6!=0]-2),
                 ((I_all_data*ch_7)[I_all_data*ch_7!=0]-2),
                 ((I_all_data*ch_8)[I_all_data*ch_8!=0]-2)]
            Q = [((Q_all_data*ch_1)[Q_all_data*ch_1!=0]-2),
                 ((Q_all_data*ch_2)[Q_all_data*ch_2!=0]-2),
                 ((Q_all_data*ch_3)[Q_all_data*ch_3!=0]-2),
                 ((Q_all_data*ch_4)[Q_all_data*ch_4!=0]-2),
                 ((Q_all_data*ch_5)[Q_all_data*ch_5!=0]-2),
                 ((Q_all_data*ch_6)[Q_all_data*ch_6!=0]-2),
                 ((Q_all_data*ch_7)[Q_all_data*ch_7!=0]-2),
                 ((Q_all_data*ch_8)[Q_all_data*ch_8!=0]-2)]

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
        ch_vec = self.ch_vec
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


    def ask_raw(self, cmd: str) -> str:
            """
            Legacy function, may not work with the current driver

            Overwriting the ask_ray qcodes native function to query binary
            Low-level interface to ``visa_handle.ask``.
            Args:
                cmd: The command to send to the instrument.
            Returns:
                str: The instrument's response.
            """
            with DelayedKeyboardInterrupt():
                keep_trying = True
                count = 0
                while keep_trying:
                    count += 1
                    self.visa_log.debug(f"Querying: {cmd}")
                    try:
                        response = self.visa_handle.query_binary_values(cmd, datatype="h", is_big_endian=False)
                        self.visa_log.debug(f"Response: {response}")
                        if len(response) > 1:
                            i = 0
                    except:
                        try:        # try to read the data as a single point of data in case buffer is empty
                            response = self.visa_handle.query_binary_values(cmd, datatype="h", is_big_endian=False, data_points=1, header_fmt='ieee', expect_termination=False)
                        except:
                            response = 'ERR'
                    if response != 'ERR' and response != [3338]:
                        keep_trying = False
                    if count>10:
                        keep_trying = False

            return response


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