"""
This module contains a function to take only the best LFEs from an LFE catalog,
download every one-minute time window where there is an LFE recorded,
and stack the signal over all the LFEs to get the template
"""

import obspy
from obspy import read_inventory
from obspy import UTCDateTime
from obspy.core.stream import Stream
from obspy.signal.cross_correlation import correlate

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle

from math import cos, floor, pi, sin, sqrt

import get_data
import get_responses

from get_stations import filter_stations
from stacking import linstack

def rotate_stream(D, orientation, reference):
    """
    """
    channels = []
    for trace in D:
        channels.append(trace.stats.channel)
    Drot = Stream()
    for (angle, trace) in zip(orientation, D):
        channel = trace.stats.channel
        trace_rot = trace.copy()
        if channel[2] in ['E', 'N', '1', '2']:
            az = angle['azimuth']
            az0 = reference[channel]['azimuth']
            if (az != az0):
                normal_dict = {'E':'N', 'N':'E', '1':'2', '2':'1'}
                normal = channel[0 : 2] + normal_dict[channel[2]]
                index = channels.index(normal)
                if channel[2] in ['E', '1']:
                    dE = angle['azimuth'] * pi / 180.0
                    dN = orientation[index]['azimuth'] * pi / 180.0
                    tE = reference[channel]['azimuth'] * pi / 180.0
                    tN = reference[normal]['azimuth'] * pi / 180.0
                    dataE = trace.data
                    dataN = D[index].data
                    trace_rot.data = cos(dE - tE) * dataE + cos(dN - tE) * dataN
                else:
                    dE = orientation[index]['azimuth'] * pi / 180.0
                    dN = angle['azimuth'] * pi / 180.0
                    tE = reference[normal]['azimuth'] * pi / 180.0
                    tN = reference[channel]['azimuth'] * pi / 180.0
                    dataE = D[index].data
                    dataN = trace.data
                    trace_rot.data = cos(dE - tN) * dataE + cos(dN - tN) * dataN
        Drot.append(trace_rot)
    return Drot

def compute_new_templates(family, latitude, longitude, catalog, filename, \
    directory, max_dist, max_LFEs, \
    TDUR, filt, dt, nattempts, waittime, method='RMS'):
    """
    This function take only the best LFEs from an LFE catalog,
    downloads every one-minute time window where there is an LFE recorded,
    and stacks the signal over all the LFEs to get the template

    Input:
        type filename = string
        filename = Name of the template
        type catalog = string
        catalog = Name of the catalog containing the LFEs
        type threshold = float
        threshold = Minimun value of cross correlation to keep LFE
        type stations = list of strings
        stations = Name of the stations where we want a template
        type TDUR = float
        TDUR = Time to add before and after the time window for tapering
        type filt = tuple of floats
        filt = Lower and upper frequencies of the filter
        type dt = float
        dt = Time step for resampling
        type nattempts = integer
        nattempts = Number of times we try to download data
        type waittime = positive float
        waittime = Type to wait between two attempts at downloading
        type method = string
        method = Normalization method for linear stack (RMS or Max)
    Output:
        None
    """
    # Get the time of LFE detections
    namefile = '../data/' + catalog + '/catalogs/' + family + '/' + filename
    LFEtime = pickle.load(open(namefile, 'rb'))
    LFEsort = LFEtime.sort_values(by=['cc'], ascending=False)

    # Get the network, channels, and location of the stations
    stations_BK = filter_stations('BK')
    stations_NC = filter_stations('NC')
    stations_PB = filter_stations('PB')
    stations = pd.concat([stations_BK, stations_NC, stations_PB], ignore_index=True)
#    stations = pd.concat([stations_BK, stations_NC], ignore_index=True)
#    stations = stations_PB

    # Create directory to store the waveforms
    namedir = directory + '/' + family
    if not os.path.exists(namedir):
        os.makedirs(namedir)

    # File to write error messages
    errorfile = 'error/' + family + '.txt'

    # Keep only stations close to LFE family
    a = 6378.136
    e = 0.006694470
    dx = (pi / 180.0) * a * cos(latitude * pi / 180.0) / sqrt(1.0 - e * e * \
        sin(latitude * pi / 180.0) * sin(latitude * pi / 180.0))
    dy = (3.6 * pi / 648.0) * a * (1.0 - e * e) / ((1.0 - e * e * sin(latitude * \
        pi / 180.0) * sin(latitude * pi / 180.0)) ** 1.5)
    x = dx * (stations['Lon'] - longitude)
    y = dy * (stations['Lat'] - latitude)
    stations['distance'] = np.sqrt(np.power(x, 2.0) + np.power(y, 2.0))
    mask = stations['distance'] <= max_dist
    stations = stations.loc[mask]

    mask = stations[stations['Lo'] == '2'].index
    stations.drop(mask, inplace=True)
    print(stations)

    # Loop over stations
    for ir in range(0, len(stations)):
        # Get station metadata for downloading
        network = stations['Net'].iloc[ir]
        station = stations['Stat'].iloc[ir]
        channels = stations['Cha'].iloc[ir]
        location = stations['Lo'].iloc[ir]
        starttime = stations['Start time'].iloc[ir]
        endtime = stations['End time'].iloc[ir]
        if (network == 'PB'):
            server = 'IRIS'
        else:
            server = 'NCEDC'

        # Filter LFEs for the period where the station was recording
        df = pd.DataFrame({'year': LFEsort['year'], 'month': LFEsort['month'], 'day': LFEsort['day']})
        date = pd.to_datetime(df)
        mask = (date >= starttime) & (date <= endtime)
        LFEsub = LFEsort.loc[mask]

        # Download instrument response
        if (server == 'IRIS'):
            get_responses.get_from_IRIS(station, network)
        elif (server == 'NCEDC'):
            get_responses.get_from_NCEDC(station, network)
        else:
            raise ValueError('You can only download data from IRIS and NCEDC')
        # Create streams
        cha_list = channels.split(',')
        streams = []
        for channel in cha_list:
            streams.append(Stream())
        # Create dictionary of channel orientations
        reference = dict.fromkeys(cha_list)
        # Initialization
        complete = False
        index = 0
        # Loop on LFEs
        while ((index < len(LFEsub)) and (complete == False)):
            mySecond = int(floor(LFEsub['second'].iloc[index]))
            myMicrosecond = int(1000000.0 * \
                (LFEsub['second'].iloc[index] - floor(LFEsub['second'].iloc[index])))
            Tori = UTCDateTime(year=LFEsub['year'].iloc[index], \
                month=LFEsub['month'].iloc[index], day=LFEsub['day'].iloc[index], \
                hour=LFEsub['hour'].iloc[index], minute=LFEsub['minute'].iloc[index], \
                second=mySecond, microsecond=myMicrosecond)
            Tstart = Tori - TDUR
            Tend = Tori + 60.0 + TDUR
            # First case: we can get the data from IRIS
            if (server == 'IRIS'):
                (D, orientation) = get_data.get_from_IRIS(station, network, channels, location, \
                    Tstart, Tend, filt, dt, nattempts, waittime, errorfile)
            # Second case: we get the data from NCEDC
            elif (server == 'NCEDC'):
                (D, orientation) = get_data.get_from_NCEDC(station, network, channels, location, \
                        Tstart, Tend, filt, dt, nattempts, waittime, errorfile)
            else:
                raise ValueError('You can only download data from IRIS and NCEDC')
            if (type(D) == obspy.core.stream.Stream):
                # Fill dictionary of channel orientations
                for channel in cha_list:
                    if reference[channel] == None:
                        mylocation = location
                        if (mylocation == '--'):
                            mylocation = ''
                        response = '../data/response/' + network + '_' + station + '.xml'
                        inventory = read_inventory(response, format='STATIONXML')
                        angle = inventory.get_orientation(network + '.' + \
                            station + '.' + mylocation + '.' + channel, Tori)
                        reference[channel] = angle
                # Rotation of components
                D = rotate_stream(D, orientation, reference)
                # Add to stream
                for (i, channel) in enumerate(cha_list):
                    if len(streams[i]) < max_LFEs:
                        Dselect = D.select(channel=channel).slice(Tori, Tori + 60.0)
                        if len(Dselect) == 1:
                            if Dselect[0].stats.npts == int(60.0 / dt + 1):
                                streams[i].append(Dselect[0])
            else:
                print('Failed at downloading data')
            # Update
            index = index + 1
            complete = True
            for (channel, stream) in zip(cha_list, streams):
                if len(stream) < max_LFEs:
                    complete = False
        # Stack
        for (channel, stream) in zip(cha_list, streams):
            if (len(stream) > 0):
                stack = linstack([stream], normalize=True, method=method) 
                savename = namedir + '/' + station + '_' + channel + '.pkl'
                pickle.dump([stack[0], reference[channel]], open(savename, 'wb'))

if __name__ == '__main__':

    # Set the parameters
    catalog = 'Ducellier'
    filename = 'catalog_2007_2009.pkl'
    directory = 'templates_2007_2009'
    max_dist = 100.0
    max_LFEs = 150
    TDUR = 10.0
    filt = (1.5, 9.0)
    dt = 0.05
    nattempts = 10
    waittime = 10.0    
    method = 'RMS'

    LFEloc = np.loadtxt('../data/Plourde_2015/templates_list.txt', \
        dtype={'names': ('name', 'family', 'lat', 'lon', 'depth', 'eH', \
        'eZ', 'nb'), \
             'formats': ('S13', 'S3', np.float, np.float, np.float, \
        np.float, np.float, np.int)}, \
        skiprows=1)

    for ie in range(61, 62): #len(LFEloc)):
        family = LFEloc[ie][0].decode('utf-8')
        print(family)
        latitude = LFEloc[ie][2]
        longitude = LFEloc[ie][3]
        compute_new_templates(family, latitude, longitude, catalog, filename, \
            directory, max_dist, max_LFEs, \
            TDUR, filt, dt, nattempts, waittime, method='RMS')
