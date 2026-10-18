"""
This module contains functions to find LFEs with the
temporary stations or with the permanent stations
using the templates from Plourde et al. (2015)
"""
import obspy
from obspy import read
from obspy import read_inventory
from obspy import UTCDateTime
from obspy.core.stream import Stream
from obspy.core.trace import Trace
from obspy.signal.cross_correlation import correlate

import multiprocessing
import numpy as np
import os
import pandas as pd
import pickle

from datetime import datetime, timedelta
from functools import partial
from math import ceil, cos, floor, pi
from multiprocessing import Pool

import correlate
from get_data import get_from_IRIS, get_from_NCEDC

def clean_LFEs(index, times, meancc, dt, freq0):
    """
    This function takes all times where the
    cross-correlation is higher than a threshold
    and groups those that belongs to the same LFE

    Input:
        type index = 1D numpy array
        index = Indices where cc is higher than threshold
        type times = 1D numpy array
        times = Times where cc is computed
        type meancc = 1D numpy array
        meancc = Average cc across all channels
        type dt = float
        dt = Time step of the seismograms
        type freq0 = float
        freq0 = Maximum frequency rate of LFE occurrence
    Output:
        type time = 1D numpy array
        time = Timing of LFEs
        type cc = 1D numpy array
        cc = Maximum cc during LFE
    """
    # Initializations
    maxdiff = int(floor(1.0 / (dt * freq0)))
    list_index = [[index[0][0]]]
    list_times = [[times[index[0][0]]]]
    list_cc = [[meancc[index[0][0]]]]
    # Group LFE times that are close to each other
    for i in range(1, np.shape(index)[1]):
        if (index[0][i] - list_index[-1][-1] <= maxdiff):
            list_index[-1].append(index[0][i])
            list_times[-1].append(times[index[0][i]])
            list_cc[-1].append(meancc[index[0][i]])
        else:
            list_index.append([index[0][i]])
            list_times.append([times[index[0][i]]])
            list_cc.append([meancc[index[0][i]]])
    # Number of LFEs identified
    N = len(list_index)
    time = np.zeros(N)
    cc = np.zeros(N)
    # Timing of LFE is where cc is maximum
    for i in range(0, N):
        maxcc =  np.amax(np.array(list_cc[i]))
        imax = np.argmax(np.array(list_cc[i]))
        cc[i] = maxcc
        time[i] = list_times[i][imax]
    return(time, cc)

def rotate_data(EW, NS, dE, dN, tE, tN, result):
    """
    """
    if len(EW) == len(NS):
        # Orientation of the data
        dE = dE['azimuth'] * pi / 180.0
        dN = dN['azimuth'] * pi / 180.0
        # Orientation of the template
        tE = tE['azimuth'] * pi / 180.0
        tN = tN['azimuth'] * pi / 180.0
        for i in range(0, len(EW)):
            if (len(EW[i].data) == len(NS[i].data)):
                dataE = EW[i].data
                dataN = NS[i].data 
                if result == 'E':
                    EW[i].data = cos(dE - tE) * dataE + \
                                 cos(dN - tE) * dataN
                else:
                    NS[i].data = cos(dE - tN) * dataE + \
                                 cos(dN - tN) * dataN
    if result == 'E':
        return(EW)
    else:
        return(NS)

def download_data(staloc, Tstart, Tend, filt, nattempts, waittime, ncpu, icpu):
    """
    """
    nstations = int(ceil(len(staloc) / ncpu))
    ibegin = icpu * nstations
    iend = min((icpu + 1) * nstations, len(staloc))

    for ir in range(ibegin, iend):
        station = staloc['station'][ir]
        network = staloc['network'][ir]
        channels = staloc['channels'][ir]
        location = staloc['location'][ir]
        server = staloc['server'][ir]
        time_on = staloc['time_on'][ir]
        time_off = staloc['time_off'][ir]
        dt = staloc['dt'][ir]

        # File to write error messages
        namedir = 'error'
        if not os.path.exists(namedir):
            os.makedirs(namedir)
        errorfile = 'error/' + station + '.txt'

        # Check whether there are data for this period of time
        year_on = int(time_on[0:4])
        month_on = int(time_on[5:7])
        day_on = int(time_on[8:10])
        year_off = int(time_off[0:4])
        month_off = int(time_off[5:7])
        day_off = int(time_off[8:10])
        if ((Tstart > UTCDateTime(year=year_on, month=month_on, day=day_on)) \
            and (Tend < UTCDateTime(year=year_off, month=month_off, day=day_off))):

            # First case: we can get the data from IRIS
            if (server == 'IRIS'):
                (D, orientation) = get_from_IRIS(station, network, channels, \
                    location, Tstart, Tend, filt, dt, nattempts, waittime, \
                    errorfile)
            # Second case: we get the data from NCEDC
            elif (server == 'NCEDC'):
                (D, orientation) = get_from_NCEDC(station, network, channels, \
                    location, Tstart, Tend, filt, dt, nattempts, waittime, \
                    errorfile)
            else:
                raise ValueError('You can only download data from IRIS and NCEDC')

            # Store the data into temporary files
            if type(D) == obspy.core.stream.Stream:
                mychannels = channels.split(',')
                for channel in mychannels:
                    stream = D.select(channel=channel)
                    if type(stream) == obspy.core.stream.Stream:
                        if len(stream) > 0:
                            stream.write('tmp/' + station + '_' + channel + \
                                '.mseed', format='MSEED')
                            namefile = 'tmp/' + station + '_' + channel + '.pkl'
                            pickle.dump(orientation, open(namefile, 'wb'))

def analyze_data(families, staloc, tbegin, tend, \
    freq0, type_threshold, threshold, ncpu, icpu):
    """
    """
    nfamilies = int(ceil(len(families) / ncpu))
    ibegin = icpu * nfamilies
    iend = min((icpu + 1) * nfamilies, len(families))
    
    # Loop on families
    for i in range(ibegin, iend):

        # Create directory to store the LFEs times
        namedir = 'LFEs/' + families['family'].iloc[i]
        if not os.path.exists(namedir):
             os.makedirs(namedir)

        # File to write number of stations
        namedir = 'nstations'
        if not os.path.exists(namedir):
            os.makedirs(namedir)
        stationfile = 'nstations/' + families['family'].iloc[i] + '.txt'

        # Create dataframe to store LFE times
        df = pd.DataFrame(columns=['year', 'month', 'day', 'hour', \
            'minute', 'second', 'cc', 'nchannel'])
 
        # Read the templates
        stations = families['stations'].iloc[i].split(',')
        templates = Stream()
        orientations = []
        names = []
        for station in stations:
            subset = staloc.loc[staloc['station'] == station]
            channels = subset['channels'].iloc[0]
            mychannels = channels.split(',')
            for channel in mychannels:
                data = pickle.load(open(template_dir + '/' + \
                    families['family'].iloc[i] + '/' + station + '_' + \
                    channel + '.pkl', 'rb'))
                template = data[0]
                angle = data[1]
                templates.append(template)
                orientations.append(angle)
                names.append(station + '_' + channel)

        # Check the time step of the stations
        subset = staloc.loc[staloc['station'].isin(stations)]
        if len(subset['dt'].value_counts()) == 1:
            dt = subset['dt'].iloc[0]
        else:
            raise ValueError('All stations must have the same time step')

        # Number of hours of data to analyze 
        t1 = UTCDateTime(year=tbegin[0], month=tbegin[1], \
            day=tbegin[2], hour=tbegin[3], minute=tbegin[4], \
            second=tbegin[5])
        t2 = UTCDateTime(year=tend[0], month=tend[1], \
            day=tend[2], hour=tend[3], minute=tend[4], \
            second=tend[5])  
        nhour = int(ceil((t2 - t1) / 3600.0))
        duration = families['duration'].iloc[i]

        # To rotate components
        swap = {'E':'N', 'N':'E', '1':'2', '2':'1'}

        # Loop on hours of data
        for hour in range(0, nhour):
            Tstart = t1 + hour * 3600.0
            Tend = t1 + (hour + 1) * 3600.0 + duration
            delta = Tend - Tstart
            ndata = int(delta / dt) + 1
            
            # Get the data
            data = []
            for station in stations:
                subset = staloc.loc[staloc['station'] == station]
                channels = subset['channels'].iloc[0]
                mychannels = channels.split(',')
                for num, channel in enumerate(mychannels):
                    try:
                        D = read('tmp/' + station + '_' + channel + '.mseed')
                        D = D.slice(Tstart, Tend)

                        if (type(D) == obspy.core.stream.Stream):
                            namefile = 'tmp/' + station + '_' + channel + \
                                '.pkl'
                            orientation = pickle.load(open(namefile, 'rb')) \
                                [num]
                            index = names.index(station + '_' + channel)
                            reference = orientations[index]

                            # Rotate components
                            if (len(mychannels) > 1) and (num < 2):
                                if orientation != reference:
                                    channel_new = channel[0:2] + \
                                        swap[channel[2]]
                                    D_new = read('tmp/' + station + '_' + \
                                        channel_new + '.mseed')
                                    D_new = D_new.slice(Tstart, Tend)
                                    namefile = 'tmp/' + station + '_' + \
                                        channel_new + '.pkl'
                                    if num == 0:
                                        orientation_new = pickle.load(open( \
                                            namefile, 'rb'))[1]
                                    else:
                                        orientation_new = pickle.load(open( \
                                            namefile, 'rb'))[0]
                                    index = names.index(station + '_' + \
                                        channel_new)
                                    reference_new = orientations[index]
                                    if channel[2] in ['E', '1']:
                                        D = rotate_data(D, D_new, \
                                            orientation, orientation_new, \
                                            reference, reference_new, 'E')
                                    else:
                                        D = rotate_data(D_new, D, \
                                            orientation_new, orientation, \
                                            reference_new, reference, 'N')

                            # Append stream to data
                            data.append(D)
                    except:
                        message = 'No data available for station {}'.format( \
                            station) + ' and channel {}'.format(channel) + \
                            ' at time {}/{}/{} - {}:{}:{}\n'.format( \
                            Tstart.year, Tstart.month, Tstart.day, \
                            Tstart.hour, Tstart.minute, Tstart.second)

            # Loop on channels
            nchannel = 0
            for j in range(0, len(data)):
                subdata = data[j]
                # Check whether we have a complete one-hour-long recording
                if (len(subdata) == 1):
                    if (len(subdata[0].data) == ndata):
                        # Get the template
                        station = subdata[0].stats.station
                        channel = subdata[0].stats.channel
                        template = templates.select(station=station, \
                            channel=channel)[0]
                        # Cross correlation
                        cctemp = correlate.optimized(template, subdata[0])
                        if (nchannel > 0):
                            cc = np.vstack((cc, cctemp))
                        else:
                            cc = cctemp
                        nchannel = nchannel + 1

            # Write number of channels
            with open(stationfile, 'a') as file:
                file.write('{} {} {} {} {}\n'.format(Tstart.year, \
                    Tstart.month, Tstart.day, Tstart.hour, nchannel))

            if (nchannel > 0):   
                # Compute average cross-correlation across channels
                if len(np.shape(cc)) == 1:
                    meancc = cc
                else:
                    meancc = np.mean(cc, axis=0)
                if (type_threshold == 'MAD'):
                    MAD = np.median(np.abs(meancc - np.mean(meancc)))
                    index = np.where(meancc >= threshold * MAD)
                elif (type_threshold == 'Threshold'):
                    index = np.where(meancc >= threshold)
                else:
                    raise ValueError( \
                        'Type of threshold must be MAD or Threshold')
                times = np.arange(0.0, np.shape(meancc)[0] * dt, dt)

                # Get LFE times
                if np.shape(index)[1] > 0:
                    (time, cc) = clean_LFEs(index, times, meancc, dt, freq0)

                    # Add LFE times to dataframe
                    i0 = len(df.index)
                    for j in range(0, len(time)):
                        timeLFE = Tstart + time[j]
                        df.loc[i0 + j] = [int(timeLFE.year), \
                            int(timeLFE.month), int(timeLFE.day), \
                            int(timeLFE.hour), int(timeLFE.minute), \
                            timeLFE.second + timeLFE.microsecond / 1000000.0, \
                            cc[j], nchannel]

        # Add to pandas dataframe and save
        namefile = 'LFEs/' + families['family'].iloc[i] + '/catalog.pkl'
        if os.path.exists(namefile):
            df_all = pickle.load(open(namefile, 'rb'))
            df_all = pd.concat([df_all, df], ignore_index=True)
        else:
            df_all = df    
        df_all = df_all.astype(dtype={'year':'int32', 'month':'int32', \
            'day':'int32', 'hour':'int32', 'minute':'int32', \
            'second':'float', 'cc':'float', 'nchannel':'int32'})
        pickle.dump(df_all, open(namefile, 'wb'))

def find_LFEs(family_file, station_file, template_dir, tbegin, tend, \
    TDUR, filt, freq0, nattempts, waittime, type_threshold='MAD', \
    threshold=0.0075):
    """    
    """
    # Number of CPUs
    ncpu = multiprocessing.cpu_count()

    # Get the network, channels, and location of the stations
    staloc = pd.read_csv(station_file, sep=r'\s{1,}', header=None, \
        engine='python')
    staloc.columns = ['station', 'network', 'channels', 'location', \
        'server', 'latitude', 'longitude', 'time_on', 'time_off', 'dt']

    # Get the families, stations, and duration of the template
    families = pd.read_csv(family_file, sep=r'\s{1,}', header=None, \
        engine='python')
    families.columns = ['family', 'stations', 'duration']

    # Begin and end time of analysis
    t1 = UTCDateTime(year=tbegin[0], month=tbegin[1], \
        day=tbegin[2], hour=tbegin[3], minute=tbegin[4], \
        second=tbegin[5])
    t2 = UTCDateTime(year=tend[0], month=tend[1], \
        day=tend[2], hour=tend[3], minute=tend[4], \
        second=tend[5])
    
    # Number of hours of data to analyze
    nhour = int(ceil((t2 - t1) / 3600.0))

    # Begin and end time of downloading
    duration = np.max(families['duration'])
    Tstart = t1 - TDUR
    Tend = t2 + duration + TDUR

    # Temporary directory to store the data
    namedir = 'tmp'
    if not os.path.exists(namedir):
        os.makedirs(namedir)

    # Download the data from the stations
    map_func = partial(download_data, staloc, Tstart, Tend, filt, \
        nattempts, waittime, ncpu)
    with Pool(ncpu) as pool:
        pool.map(map_func, iter(range(0, ncpu)))

    # Analyze seismic data
    map_func = partial(analyze_data, families, staloc, tbegin, tend, \
        freq0, type_threshold, threshold, ncpu)
    with Pool(ncpu) as pool:
        pool.map(map_func, iter(range(0, ncpu)))
    
if __name__ == '__main__':

    # Set the parameters
    family_file = '../data/Ducellier/families_permanent.txt'
    station_file = '../data/Ducellier/stations_permanent.txt'
    template_dir = 'templates_v1'
    TDUR = 10.0
    filt = (1.5, 9.0)
    freq0 = 1.0
    nattempts = 10
    waittime = 10.0
    type_threshold = 'MAD'
    threshold = 8.0

    begin = datetime.now()

    # August 2010
    year = 2010
    month = 8
    for day in range(28, 32):
        for hour in range(0, 24, 12):
            tbegin = (year, month, day, hour, 0, 0)
            if (hour == 12):
                if (day == 31):
                    tend = (year, month + 1, 1, 0, 0, 0)
                else:
                    tend = (year, month, day + 1, 0, 0, 0)
            else:
                tend = (year, month, day, hour + 12, 0, 0)
            find_LFEs(family_file, station_file, template_dir, tbegin, tend, \
                TDUR, filt, freq0, nattempts, waittime, type_threshold, \
                threshold)

    families = pd.read_csv(family_file, sep=r'\s{1,}', header=None, \
         engine='python')
    families.columns = ['family', 'stations', 'duration']
    for i in range(0, len(families)):
        os.rename('LFEs/' + families['family'].iloc[i] + '/catalog.pkl', \
            'LFEs/' + families['family'].iloc[i] + \
            '/catalog_{:04d}_{:02d}'.format(year, month) + '.pkl')

    end = datetime.now()
    duration = end - begin
    duration.total_seconds()
    print(duration)
