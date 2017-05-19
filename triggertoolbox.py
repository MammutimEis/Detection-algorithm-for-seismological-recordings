# -*- coding: utf-8 -*-
# ------------------------------------------------------------------
# Filename: staltatoolbox.py
#  Purpose: Functions for detection algorithms
#   Author: Miro Pütz
#    Email: miropuetz@posteo.net
# --------------------------------------------------------------------

import obspy as obs
from obspy.signal.trigger import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
#mport scipy
#mport math
plt.close('all')
plt.matplotlib.rcParams.update({'font.size': 15})

def readFiles(paths, dayFile, wholeDay=True, start=None, end=None):
    """
    Read data of all stations in. FLBP1, FLBP2, FLBP3, 31

    The paths to the data files and which day (dayFile) have to be specified

    :param paths: All the paths to the dataFiles
    :type  paths: List of Strings
    :param dayFile: Dayfile that is read in
    :type  dayFile: String
    :param wholeDay: If set to true, read in the whole day otherwise specify the start and end time
    :type  wholeDay: Boolean
    :param start: starttime of the datafile
    :type  start: obspy.core.utcdatetime.UTCDateTime
    :param end: Endtime of the dayfile
    :type  end: obspy.core.utcdatetime.UTCDateTime
    :return: Returns the recorded data and the synthetic data as a stream object in a tuple
    :type return: (obspy.core.stream.Stream, obspy.core.stream.Stream)
     """
    if wholeDay:
        st = obs.read(paths[0] + dayFile)  # read in FLBP* stations
        st += obs.read(paths[1] + dayFile)     # read in 31** stations
        st += obs.read(paths[2] + dayFile)   # read in 3G13 station
        stSyn = obs.read(paths[3])  # read in synthetic data
    else:
        st = obs.read(paths[0] + dayFile, starttime=start, endtime=end) # read in FLBP* stations
        st += obs.read(paths[1]+ dayFile, starttime=start, endtime=end)    # read in 31** stations
        st += obs.read(paths[2] + dayFile, starttime=start, endtime=end)  # read in 3G13 station
        stSyn = obs.read(paths[3]) # read in synthetic data

    print("\nData is read in:\n")
    print(st.__str__(extended=True))
    print("\nsynthetic Data is read in:\n ")
    print(stSyn.__str__(extended=True))
    return st, stSyn



def spectogram(trace, traceSyn):
    '''
    Plots the amplitude spectrum of a recorded trace and a sytnthetic trace.

    Calculation with the Fast Fourier Transform (FFT) from numpy. The x-axis is converted to
    frequencies from 0 to half of the maximum sampling frequency (Nyquist frequency). The y-axis
    is converted by taking the absolute value of the FFT and multiplication by a factor of 2.
    Further the amplitude has to be normalized with the number of samples. To improve the amplitude
    spectrum the data is also multiplied with a Hanning window. Balzer explains how to apply the FFT in his
    IPython Notebook: rhttp://nbviewer.jupyter.org/github/balzer82/FFT-Python/blob/master/FFT-Tutorial.ipynb

    :param trace: The recorded trace
    :type trace: obspy.core.trace.Trace
    :param traceSyn: The synthetic trace
    :type traceSyn: obspy.core.trace.Trace
    '''
    #spektrum
    trace.data -= trace.data.mean()
    traceSyn.data -= traceSyn.data.mean()

    Fs = trace.stats.sampling_rate

    fig = plt.figure(1)
    fig.suptitle(trace.stats.station, fontsize=20)
    ax1 = plt.subplot(211)
    ax1.set_title('Recorded Signal')

    # Calculate amplitude spectrum for recorded data
    y = trace.data
    hann = np.hanning(len(y))
    N = len(y)/2 +1
    Y = np.fft.fft(y*hann)
    X = np.linspace(0, Fs/2, N, endpoint = True)

    plt.ylabel('|Ground Motion (m/s)|')
    plt.semilogx(X, (2.0*np.abs(Y[:N])/N), 'k')    # plotting the spectrum
    plt.setp(ax1.get_xticklabels(), visible=False)

    # Calculate amplitude spectrum for synthetic data
    y = traceSyn.data
    hann = np.hanning(len(y))
    N = len(y)/2 +1
    Y = np.fft.fft(y*hann)
    X = np.linspace(0, Fs/2, N, endpoint = True)

    ax2 = plt.subplot(212, sharex=ax1)
    ax2.set_title('Synthetic Signal')
    plt.semilogx(X, (2.0*np.abs(Y[:N]))/N) # plotting the spectrum
    plt.xlabel('Freq (Hz)')
    plt.ylabel('|Ground Motion (m/s)|')
    plt.yticks(np.array([0,1,2,3,4,5])*10e-8)
    plt.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    fig.subplots_adjust(top=0.85)
    plt.savefig('../spectrum.pdf', bbox_inches='tight', pad_inches=0.05)
    plt.draw()



def recursiveStaLtaCharacteristic(trace, staDur, ltaDur, thrOn, thrOff):
    '''
    Plot a trace and the corresponding characteristic function of the recursive STA/LTA trigger

    :param trace: The recorded trace
    :type trace: obspy.core.trace.Trace
    :param staDur: STA duration in seconds
    :param ltaDur: LTA duration in seconds
    :param thrOn: Trigger on threshold
    :parm thrOff: Trigger off threshold
    '''
    df = trace.stats.sampling_rate
    cft = obs.signal.trigger.recursive_sta_lta(trace.data, int(staDur*df), int(ltaDur*df))
    plot_trigger(trace, cft, thr_on=thrOn, thr_off=thrOff, show=False)
    plt.draw()


    
def getEarliestTrace(stream):
    '''
    Returns the trace in a stream with the earliest starttime

    :param stream: A stream of seismic traces
    :type  stream: obspy.core.stream.Stream
    :return: The trace with the earliest starttime
    :type return: obspy.core.stream.Stream
    '''
    earliestStarttime = stream[0].stats.starttime
    earliestTrace = stream[0]
    for trace in stream:
        if trace.stats.starttime < earliestStarttime:
            earliestStarttime = trace.stats.starttime
            earliestTrace = trace

    return earliestTrace
        


def getCoincidenceEventTime(stream, coincidenceEvents):
    '''
    The time and duration of all events is put into a list. For each event the time difference
    to the earliest starttime and the durationis of an event is put into a list. This results in
    a list for each event and one list that contatins all lists.
    
    :param stream: A stream of seismic traces
    :type  stream: obspy.core.stream.Stream
    :param coincidenceEvents: The result of the coincidence trigger, a list of coincidence events.
                              For each event a dictionary exists which provides the start time, duration
                              and which stations have triggered.
    :type  coincidenceEvents: List with dictionarys; [{}, {}]
    :return: The relative times to the earliest starttime and the durations of all coincidence events
    :type return: List in List; [[time, duration], [time,duration]]
    '''
    referenceTrace = getEarliestTrace(stream)
    eventTime = []
    for i in range(len(coincidenceEvents)):
        eventTime.append([coincidenceEvents[i]['time']-referenceTrace.stats.starttime, coincidenceEvents[i]['duration']])
    return eventTime



def plotCoincidence(stream, eventTime):
    '''
    Plots the times and durations of events on top of all traces.
    
    :param stream: A stream of recorded traces
    :type stream: obspy.core.stream.Stream
    :param eventTime: A list of the relative times to the earliest starttime of all events
    :type  eventTime: List; [] 
    '''
    numOfStations = len(stream)
    fig = stream.plot(type = 'relative', handle=True, equal_scale=True, method='full')

    for i in range(numOfStations):
        ax = fig.axes[i]
        for event in eventTime:
            #start
            ax.vlines(event[0], ax.get_ylim()[0], ax.get_ylim()[1], lw=1, color='b')
            #duration
            ax.vlines(event[0]+event[1], ax.get_ylim()[0], ax.get_ylim()[1], lw=1, color='r')

    plt.draw()
    
def plotRecursiveStaLta(stream, staDur, ltaDur, thrOn, thrOff):
    '''
    Applies the recursive STA/LTA trigger to each trace and plots the triggered events
    on top of the seismograms. If the traces do not have the same start time the difference
    to the earliest starttime has to be added to consider the time differences.
    
    :param stream: A stream of recorded traces
    :type stream: obspy.core.stream.Stream
    :param staDur: STA duration in seconds
    :param ltaDur: LTA duration in seconds
    :param thrOn: Trigger on threshold
    :parm thrOff: Trigger off threshold
    '''
    fig = stream.plot(type = 'relative', handle=True, equal_scale=True, method='full', size=(800,500))
    numOfStations = len(stream)
    referenceTrace = getEarliestTrace(stream)

    for i in range(numOfStations):
        #First calculat the characteristic function
        df = stream[i].stats.sampling_rate
        cft = obs.signal.trigger.recursive_sta_lta(stream[i].data, int(staDur*df), int(ltaDur*df))

        # Array of all triggered events in the unit samples 
        on_off = np.array(obs.signal.trigger.trigger_onset(cft, thres1=thrOn, thres2=thrOff), dtype='float64')

        # Conversion in TrigerTimes
        on_off /= df 

        ax = fig.axes[i]

        # If traces have a time offset, correction is needed
        timeDiff = stream[i].stats.starttime - referenceTrace.stats.starttime

        for event in on_off:
            #start
            ax.vlines(event[0]+timeDiff, ax.get_ylim()[0], ax.get_ylim()[1], lw=1, color='b')
            #end
            ax.vlines(event[1]+timeDiff, ax.get_ylim()[0], ax.get_ylim()[1], lw=1, color='r')

    plt.draw()



def addSyntheticEvents(stream, syntheticStream, synTimes):
    '''
    Superimpose recorded seismic traces with synthetic events.

    :param stream: A stream of recorded traces
    :type stream: obspy.core.stream.Stream
    :param syntheticStream: A stream of synthetic traces
    :type stream: obspy.core.stream.Stream
    :param synTimes: List of all synthetic times in samples where an synthetic event is added
    '''
    numOfStations = len(stream)
    for dt in synTimes[:-1]:
        for i in range(numOfStations):
            window = syntheticStream[i].stats.npts
            stream[i].data[dt:window+dt] += syntheticStream[i].data

    print("Superimposed %s synthetic events with the data" %(len(synTimes[:-1])))



def syntheticTimes(samplingRate, offset, interval, duration):
    '''
    Provide the samples where synthetic events are superimposed with the data. The inerval
    of the added synthetic events and the period in which synthetic events are added
    can be specified. All traces need to have the same sampling rate.
    
    :param samplingRate: sampling rate of the traces
    :typ   samplingRate: float
    :param offset: Offset in seconds. It takes some seconds until the synthetic event starts. 
                   The offset is required for the event statistic. It expects the events
                   at mutiples of intervall*60
    :type  offset: float, unit seconds
    :param interval: Time between two synthetic events in minutes
    :param duration: Time in hours in which synthetic events are added
    :type  duration: int
    :return: List of the times in samples where synthetic events are superimposed with 
             recorded seismic data.
    '''
    # start and max sample for synthetic events
    samplingRate = int(samplingRate)
    maxSample = int(samplingRate * 3600 * duration)      
    startSample = (interval*60-offset)*samplingRate

    # Array with the times in samples
    synTimes = np.arange(startSample, maxSample, interval*60*samplingRate)

    return synTimes



def eventStatistic(referenceTime, coincidenceEvents, interval, synTimes):
    '''
    Checks all time times when real events are triggered and looks for matches with 
    synthetic events.

    :param reference time: For the calculation of the difference between synthetic time and real
                           triggered event time. Should be old startTime.
    :type  reference time: obs.core.utcdatetime.UTCDateTime
    :type stream: obspy.core.stream.Stream
    :param coincidenceEvents: The result of the coincidence trigger, a list of coincidence events.
    :type  coincidenceEvents: List with dictionarys; [{}, {}]
    :param interval: Time between two synthetic events in minutes
    :param synTimes: List of all synthetic times in samples where an synthetic event is added
    '''
    #initialization
    count = 0       # counts matches
    j = interval*60 # the expected time of an event
    i = 0           # the index of the coincidence event
    
    while True:
        if j*200 >= synTimes[-1]:
            #break loop when expected time > last synthetic event
            break
        try:
            # difference between synthetic time und real triggered event time
            triggerDiff = coincidenceEvents[i]['time'] - (referenceTime + j)
        # index error, all real triggered events tested, break loop
        except IndexError: 
            break
        
        # match if trigger difference < 3 seconds or event duration longer than trigger difference.
        if ( abs(triggerDiff) < 3 ) or ( triggerDiff < 0 and coincidenceEvents[i]['duration'] > abs(triggerDiff) ):
            count += 1
            j += interval*60

        # synthetic event has not been triggered
        elif triggerDiff >= 3:
            j += interval*60
            i -= 1
        # synthetic event is not yet triggered
        i += 1

    print('%i events of %i found' %(count, len(synTimes)-1))
    print('%i other events found' %(len(coincidenceEvents)-count))
    try:
        ratio = (count*100.) / len(coincidenceEvents)
    except ZeroDivisionError:
        ratio = 0
    print('%i percent ratio\n' %(ratio))
    return ratio, count



def plotFreqResults(resultRatio, resultEvents):
    '''
    plots the results of different filter corner frequencies of the event statistic

    :param resultRatio: Matrix with the success rate of different corner frequencies
    :type  resultRatio: numpy.ndarray of size (10,30)
    :param resultEvents: Matrix with the synthetic events found of different corner frequencies
    :type  resultRatio: numpy.ndarray of size (10,30)
    '''
    
    fig, axes = plt.subplots(2,1, figsize=(15,15)) 
    
    axes[0].set_title('Success rate')
    cax = axes[0].matshow(resultRatio, vmin=0, vmax=100)
    fig.colorbar(cax, ax=axes[0], ticks=[0, 10,20,30,40,50,60,70,80,90,100], orientation='vertical', label='Percentage')
    axes[0].set_xticks(np.linspace(0,28,15))
    axes[0].xaxis.set_ticks_position('bottom')
    axes[0].set_ylabel('min frequency [Hz]')

    axes[1].set_title('Triggered Synthetic Events')
    cax = axes[1].matshow(resultEvents, vmin=0, vmax=59)
    fig.colorbar(cax, ax=axes[1], orientation='vertical', label='Count')
    axes[1].set_xticks(np.linspace(0,28,15))
    axes[1].xaxis.set_ticks_position('bottom')
    axes[1].set_xlabel('max frequency [Hz]')
    axes[1].set_ylabel('min frequency [Hz]')
 
    fig.subplots_adjust(hspace=0.05, top=0.7)
    plt.draw()    


    
def plotFreqHighResultsn():
    lowcount = -1
    highcount = -1
    resultMatrix = np.zeros((8,15))
    resultEvents = np.copy(resultMatrix)
    cStream = st.select(station='FLBP*')
    print cStream
    for low in range(20, 60, 5):
        lowcount += 1
        highcount = lowcount-1
        for high in range(low+5, 100, 5):
            highcount += 1
            stCopy = cStream.copy()
             # Filtern der Daten
            for tr in stCopy:
                tr.filter('bandpass', freqmin=low, freqmax=high, corners=4, zerophase=True)
            print 'i = '+ str(low) + ' j = ' + str(high)
            zTrig = coincidenceTrigger(stCopy, thrOn, thrOff)
            #zTime = getEventTime(zTrig, st)


            ########## Event Statistik ##########

            #Initialisierung
            count = 0
            j = m*60
            i = 0

            while j < maxSample:
                # Differenz zwischen Erwartung und tatsaechlicher Trigger Zeit
                try:
                    triggerDiff = zTrig[i]['time'] - (st[0].stats.starttime + j)
                # Index Error, alle Triggers durch, fertig, Schleifenabbruch
                except IndexError: 
                    print('%i events von %i gefunden' %(count, len(triggerTimes)-1))
                    print('%i andere Events gefunden' %(i-count))
                    try:
                        ratio = (count*100.)/i
                    except ZeroDivisionError:
                        ratio = 0
                    print('%i percent ratio\n' %(ratio))
                    resultMatrix[lowcount,highcount] = ratio
                    resultEvents[lowcount,highcount] = count
                    break

               # print i, triggerDiff
                # Treffer: Trigger Differenz < 3 sekunden oder Eventdauer geht ueber die Erwartung; naechster Trigger und naechste Erwartung
                if ( abs(triggerDiff) < 3 ) or ( triggerDiff < 0 and zTrig[i]['duration'] > abs(triggerDiff) ):
                    count += 1
                    j += m*60
                    #print('[+] TriggerEvent %i Uebereinstimmung' %(i))

                # Erwartung nicht getroffen, gleicher Trigger, naechste Erwartung 
                elif triggerDiff >= 3:
                    j += m*60
                    i -= 1
                # Erwartung NOCH nicht getroffen, naechster Trigger, gleiche Erwartung
                i += 1


def plotHighFreqResults(resultRatio, resultEvents):
    '''
    plots the results of different filter corner frequencies of the event statistic

    :param resultRatio: Matrix with the success rate of different corner frequencies
    :type  resultRatio: numpy.ndarray of size (8,15)
    :param resultEvents: Matrix with the synthetic events found of different corner frequencies
    :type  resultRatio: numpy.ndarray of size ()
    '''
    fig, axes = plt.subplots(2,1, figsize=(15,15)) 
    
    axes[0].set_title('Success rate')
    x_tick_list = np.arange(15)
    y_tick_list = np.arange(8)
    x_label_list = map(lambda x: str(25+5*x), x_tick_list)
    y_label_list = map(lambda x: str(20+5*x), y_tick_list)
    cax = axes[0].matshow(resultRatio, vmin=0, vmax=100)
    fig.colorbar(cax, ax=axes[0], ticks=[0, 10,20,30,40,50,60,70,80,90,100], orientation='vertical', label='Percentage')
    axes[0].set_yticks(y_tick_list)
    axes[0].set_yticklabels(y_label_list)
    axes[0].set_xticks(x_tick_list)
    axes[0].set_xticklabels(x_label_list)
    axes[0].xaxis.set_ticks_position('bottom')
    axes[0].set_ylabel('min frequency [Hz]')

    axes[1].set_title('Triggered Synthetic Events')
    x_tick_list = np.arange(15)
    y_tick_list = np.arange(8)
    x_label_list = map(lambda x: str(25+5*x), x_tick_list) 
    y_label_list = map(lambda x: str(20+5*x), y_tick_list)
    cax = axes[1].matshow(resultEvents, vmin=0, vmax=59)
    fig.colorbar(cax, ax=axes[1], orientation='vertical', label='Count')
    axes[1].set_yticks(y_tick_list)
    axes[1].set_yticklabels(y_label_list)
    axes[1].set_xticks(x_tick_list)
    axes[1].set_xticklabels(x_label_list)
    axes[1].xaxis.set_ticks_position('bottom')
    axes[1].set_xlabel('max frequency [Hz]')
    axes[1].set_ylabel('min frequency [Hz]')
 
    fig.subplots_adjust(hspace=0.05, top=0.9)
    plt.draw()    


    
def recursiveStaLtaCharacteristicParameters(trace):
    '''
    Plot a trace and the corresponding characteristic functions of the recursive STA/LTA trigger
    for different STA and LTA window lengths.

    :param trace: The recorded trace
    :type trace: obspy.core.trace.Trace
    :param staDur: STA duration in seconds
    :param ltaDur: LTA duration in seconds
    :param thrOn: Trigger on threshold
    :parm thrOff: Trigger off threshold
    '''
    df = trace.stats.sampling_rate
    npts = trace.stats.npts
    t = np.arange(npts, dtype=np.float32) / df
    cft = []

    # STA LTA duration parameter variation
    staLta = [
        (1, 3),
        (1, 10),
        (1, 20),
        (0.01, 1),
        (0.1, 10),
        (0.5, 10),
        (0.5, 5),
        (2, 20),
        (4, 20)]

    for i in range(9):
        cft.append(obs.signal.trigger.recursive_sta_lta(trace.data, int(staLta[i][0]*df), int(staLta[i][1]*df)))
    fig, axes = plt.subplots(10, 1, figsize=(20, 20), subplot_kw={'yticks': []})

    # plot recorded trace
    axes[0].plot(t, trace.data, 'k')
    axes[0].yaxis.set_ticks(np.linspace(-2e-5, 2e-5, 3))
    axes[0].set_ylim(-0.00003, 0.00003)
    axes[0].yaxis.set_major_formatter(FormatStrFormatter('%.E'))
    axes[0].set_ylabel('m/s')

    # plot chracteristic functions
    i = 0
    for ax in axes[1:]:
        #ax.add_subplot(10, 1, i)
        ax.plot(t, cft[i], 'k')
        ax.text(0.075, 0.8, 'sta='+ str(staLta[i][0]) + ', lta=' + str(staLta[i][1]), ha='center', va='center', transform=ax.transAxes)
        ax.yaxis.set_ticks(np.linspace(0, int(np.floor(np.max(cft[i]))), 3))
        ax.set_ylim(0, np.max(cft[i]+1))
        i += 1

    axes[9].set_xlabel("Time after %s [s]" % trace.stats.starttime.isoformat())

    # no x ticks
    for i in range(9):
        axes[i].get_xaxis().set_ticks([])

    fig.suptitle(trace.id)
    plt.draw()


    
