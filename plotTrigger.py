# -*- coding: utf-8 -*-
# ------------------------------------------------------------------
# Filename: plotCharacteristicFunction.py
#  Purpose: Plots the recursive STA/LTA on top of seismograms
#   Author: Miro Pütz
#    Email: miropuetz@posteo.net
# --------------------------------------------------------------------

import triggertoolbox as tool
import obspy as obs
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Frequency of the band-pass filter, thresholds and durations of the STA and LTA can be specified as command line arguments
parser = argparse.ArgumentParser(description='Set parameters')

parser.add_argument('-f', '--frequency',
                    nargs='+', dest='filterFreq',
                    help='give frequencies for bandpass filtering. Example: "python motionTrigger.py -f 1 20", for sampling between 1 Hz and 20 Hz',
                    type=int, default=[1, 10])

parser.add_argument('-t', '--threshold',
                    nargs='+',
                    dest='thr',
                    help='give thresholds for STA/LTA: "python motionTrigger.py -t 1 5", for trigger threshold of 1 and detrigger threshold of 5',
                     type=float, default=[3., 0.7])

parser.add_argument('-d', '--duration',
                    nargs='+',
                    dest='duration',
                    help='give duration of STA and LTA window in seconds: "python motionTrigger.py -d 1 10", for sta window of 1 sec and lta window of 10 sec',
                    type=float, default=[0.5, 10])

args = parser.parse_args()
freqMin, freqMax = args.filterFreq[0], args.filterFreq[1]
thrOn, thrOff = args.thr[0], args.thr[1]
staDur, ltaDur = args.duration[0], args.duration[1]


# specify the data to be read in
wholeDay = False                       # read in the whole day
start = obs.core.utcdatetime.UTCDateTime(2017, 2, 5, 15, 0 )  # if whole day == false, specify start time
end = obs.core.utcdatetime.UTCDateTime(2017, 2, 5, 16, 0)    # if whole day == false, specify end time
dayFile = '036'                        # dayFile that is read in

# select data
comp = 'Z'                             # selcet component


#### MAIN PROGRAM ####

# pathes to the data
ST_FLBP_PATH = '/home/zmaw/u300596/seismoData/2017/FLBP*/HH*.D/GR.FLBP*..HH*.D.2017.'
ST_31_PATH = '/home/zmaw/u300596/seismoData/2017/31*/HH*.D/GR.31*..HH*.D.2017.'
ST_3G13_PATH = '/home/zmaw/u300596/seismoData/2017/3G13/HH*.D/GR.3G13..HH*.D.2017.'
SYN_PATH = '/home/zmaw/u300596/Detection-algorithm-for-seismological-recordings-master/programs/synth_seism_HHZ.mseed'

paths = [ST_FLBP_PATH, ST_31_PATH, ST_3G13_PATH, SYN_PATH]


# read the data in
st, stSyn = tool.readFiles(paths, dayFile, wholeDay, start, end)

# select only one component
st = st.select(component=comp)       # select seismological data of one component
stSyn = stSyn.select(component=comp) # select synthetic data 

# sort, to have seismological data and synthetic data in the same order
st.sort()
stSyn.sort()

# traces need higher precision for multiplication constant
for trace in st:
    trace.data = trace.data.astype('float64')

# unit conversion: counts --> ground movement (m/s)
# traces are multiplied by their system sensitivity
st.select(station='FLBP1')[0].data *=  3.2667e-11
st.select(station='FLBP2')[0].data *= 2.8824e-10
st.select(station='FLBP3')[0].data *=  3.2667e-11
st.select(station='3171')[0].data *= 4.2667e-10
st.select(station='3172')[0].data *= 4.2667e-10
st.select(station='3G13')[0].data *= 4.2667e-10

# correct the data shift of station 3172 and 3G13
st.select(station='3172')[0].data -= st.select(station='3172')[0].data.mean()
st.select(station='3G13')[0].data -= st.select(station='3G13')[0].data.mean()

# unit conversion: displacement (m) --> ground movement (m/s)
for trace in stSyn:
    trace.data = np.diff(trace.data)*trace.stats.sampling_rate

# Obspy frequency filter
# freqMin and freqMax can be specified over command line parameters.
for trace in st:
    trace.filter('bandpass',  freqmin=freqMin, freqmax=freqMax, corners=4, zerophase=True)


### Single recursive STA/LTA trigger applied to multiple stations ###
streamCopy = st.copy() # keep the original data
tool.plotRecursiveStaLta(stream=streamCopy, staDur=staDur, ltaDur=ltaDur, thrOn=thrOn, thrOff=thrOff)


### recursive STA/LTA coincidence trigger

# Correct the starttimes for better results of the coincidence trigger
streamCopy = st.copy()  # keep the original data
streamCopy.select(station='3171')[0].stats.starttime -= 2.7 # in seconds
streamCopy.select(station='3172')[0].stats.starttime -= 1.2
streamCopy.select(station='3G13')[0].stats.starttime -= 1.2

# coincidence trigger with a coincidence sum of two
coincidenceEvents = obs.signal.trigger.coincidence_trigger(trigger_type="recstalta", thr_on=thrOn,
                                                thr_off=thrOff, stream=streamCopy,
                                                thr_coincidence_sum=2, sta=staDur, lta=ltaDur, details=False)

# time and duration of all the coincidence events
time = tool.getCoincidenceEventTime(streamCopy, coincidenceEvents)

# plot
tool.plotCoincidence(streamCopy, time)

# show the plots
plt.show()

