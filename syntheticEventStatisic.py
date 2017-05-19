 
# -*- coding: utf-8 -*-
# ------------------------------------------------------------------
# Filename: syntheticEventStatistic.py
#  Purpose: A trigger statistic with synthetic events 
#   Author: Miro Pütz
#    Email: miropuetz@posteo.net
# --------------------------------------------------------------------

import triggertoolbox as tool
import obspy as obs
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

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
start = obs.core.utcdatetime.UTCDateTime(2017, 2, 5, 7, 0 )  # if whole day == false, specify start time
end = obs.core.utcdatetime.UTCDateTime(2017, 2, 5, 12, 0)    # if whole day == false, specify end time
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


### Add Synthetic Events

# Specify synthetic event times
duration = 5 # add synthetic events over a duration of 5 hours
interval = 5 # every 5 minutes

streamCopy = st.copy() # keep the original data

# 2 seconds time offset because it takes some seconds until the
# synthetic event starts. The offset is required for the event statistic.
# The event statistic expects the events at mutiples of intervall*60
timeOffset = 2

# get the times in samples where synthetic events are superimposed with the recorded seismic
# data. Add synthetic Events every 5 minutes over a duration of 5 hours
synTimes = tool.syntheticTimes(samplingRate=streamCopy[0].stats.sampling_rate, offset=timeOffset,
                              interval=5, duration=duration)

# add the synthetic events
tool.addSyntheticEvents(stream=streamCopy, syntheticStream=stSyn, synTimes=synTimes)

# METHOD 1
######################## recursive STA/LTA coincidence trigger statistic #########################################

print("EVENT STATISTIC\n")
stream1 = streamCopy.copy()  # use the unfiltered stream superimposed with synthetic events

# Obspy frequency filter
# freqMin and freqMax can be specified over command line parameters.
for trace in stream1:
    trace.filter('bandpass',  freqmin=freqMin, freqmax=freqMax, corners=4, zerophase=True)

# Correct the starttimes for better results of the coincidence trigger
stream1.select(station='3171')[0].stats.starttime -= 2.7 # in seconds
stream1.select(station='3172')[0].stats.starttime -= 1.2
stream1.select(station='3G13')[0].stats.starttime -= 1.2

# plot recursive STA/LTA triggered events
tool.plotRecursiveStaLta(stream=stream1, staDur=staDur, ltaDur=ltaDur, thrOn=thrOn, thrOff=thrOff)

# select stations. In this example only station FLBP1 so coincidence sum is one
stream1 = stream1.select(station='FLBP1')

# coincidence trigger with a coincidence sum of two
trig = obs.signal.trigger.coincidence_trigger(trigger_type="recstalta", thr_on=thrOn,
                                              thr_off=thrOff, stream=stream1,
                                              thr_coincidence_sum=1, sta=staDur, lta=ltaDur, details=False)

# time and duration of all the coincidence events
time = tool.getCoincidenceEventTime(stream1, coincidenceEvents=trig)

# plot the coincidence events
tool.plotCoincidence(stream1, time)

# event statistic
tool.eventStatistic(referenceTime=start, coincidenceEvents=trig, interval=interval, synTimes=synTimes)

plt.show()
# METHOD 2
################################ Frequency variation ############################################################
# between 1 and 30 Hz

print("EVENT STATISTIC WITH FREQUENCY VARIATION\n")
stream2 = streamCopy.copy() # use the unfiltered stream superimposed with synthetic events

# Correct the starttimes for better results of the coincidence trigger
stream2.select(station='3171')[0].stats.starttime -= 2.7 # in seconds
stream2.select(station='3172')[0].stats.starttime -= 1.2
stream2.select(station='3G13')[0].stats.starttime -= 1.2

resultRatio = np.zeros((10,30))
resultEvents = np.copy(resultRatio)

# select stations. CoincidenceTrigger with the FLBP stations and a coincidence threshold of two
stream2 = stream2.select(station='FLBP*')

for minFreq in range(1, 10, 1):
    for maxFreq in range(minFreq+1, 30, 1):
        stCopy = stream2.copy()  # copy so data is not filtered multiple times
        # filter in dependence of minFreq and highFreq
        for tr in stCopy:
                tr.filter('bandpass', freqmin=minFreq, freqmax=maxFreq, corners=4, zerophase=True)
        print 'min = '+ str(minFreq) + ' max = ' + str(maxFreq)
        trig = obs.signal.trigger.coincidence_trigger(trigger_type="recstalta", thr_on=thrOn,
                                                            thr_off=thrOff, stream=stCopy,
                                                            thr_coincidence_sum=2, sta=staDur, lta=ltaDur, details=False)
        ratio, count = tool.eventStatistic(referenceTime=start, coincidenceEvents=trig, interval=interval, synTimes=synTimes)
        resultRatio[minFreq, maxFreq] = ratio
        resultEvents[minFreq, maxFreq] = count

# plot the results
tool.plotFreqResults(resultRatio, resultEvents)

# show the plots
plt.show()


# METHOD 3
################################ High Frequency variation ############################################################
# between 20 and 95 Hz

print("EVENT STATISTIC WITH HIGH FREQUENCY VARIATION\n")
stream3 = streamCopy.copy() # use the unfiltered stream superimposed with synthetic events

# Correct the starttimes for better results of the coincidence trigger
stream3.select(station='3171')[0].stats.starttime -= 2.7 # in seconds
stream3.select(station='3172')[0].stats.starttime -= 1.2
stream3.select(station='3G13')[0].stats.starttime -= 1.2

resultRatio = np.zeros((8,15))
resultEvents = np.copy(resultRatio)

# select stations. CoincidenceTrigger with the FLBP stations and a coincidence threshold of two
stream3 = stream3.select(station='FLBP1')

minCount = -1
maxCount = -1
for minFreq in range(20, 60, 5):
    minCount += 1
    maxCount = minCount-1
    for maxFreq in range(minFreq+5, 100, 5):
        maxCount += 1
        stCopy = stream3.copy()  # copy so data is not filtered multiple times
        # filter in dependence of minFreq and highFreq
        for tr in stCopy:
                tr.filter('bandpass', freqmin=minFreq, freqmax=maxFreq, corners=4, zerophase=True)
        print 'min = '+ str(minFreq) + ' max = ' + str(maxFreq)
        trig = obs.signal.trigger.coincidence_trigger(trigger_type="recstalta", thr_on=thrOn,
                                                            thr_off=thrOff, stream=stCopy,
                                                            thr_coincidence_sum=1, sta=staDur, lta=ltaDur, details=False)
        ratio, count = tool.eventStatistic(referenceTime=start, coincidenceEvents=trig, interval=interval, synTimes=synTimes)
        resultRatio[minCount, maxCount] = ratio
        resultEvents[minCount, maxCount] = count

# plot the results
tool.plotHighFreqResults(resultRatio, resultEvents)

# show the plots
plt.show()
