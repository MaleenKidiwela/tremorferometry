"""
This module convert horizontal components into radial and transverse,
cross-correlate templates between components and between stations
and returns the relative travel times
"""

import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np
import obspy
import os
import pandas as pd
import pickle

from math import atan2, cos, pi, sin, sqrt

def compute_traveltimes(family, stations, channels, latitudes, longitudes, \
        nt, dt, lat0, lon0):
    """
    type family = String
    family = Name of the LFE family
    type stations = List of strings
    stations = Names of the stations to be used
    type channels = List of strings
    channels = Channels to be used
    type latitudes = List of floats
    latitudes = List of latitudes of the stations
    type longitudes = List of floats
    longitudes = List of longitudes of the stations
    type nt = Integer
    nt = Length of templates
    type dt = Float
    dt = Time step for templates
    type lat0 = Float
    lat0 = Latitude of the epicenter of the LFE family
    type lon0 = Float
    lon0 = Longitude of the epicenter of the LFE family
    """
    # To compute distances
    a = 6378.136
    e = 0.006694470
    # Get the number of stations with horizontal and vertical components
    nsta3 = 0
    nsta1 = 0
    for (station, channel) in zip(stations, channels):
        if len(channel.split(',')) == 3:
            nsta3 = nsta3 + 1
        else:
            nsta1 = nsta1 + 1
    # Initialization
    vertical = np.zeros((nt, nsta1 + nsta3))
    radial = np.zeros((nt, nsta3))
    transverse = np.zeros((nt, nsta3))
    sta_V = []
    sta_H = []
    index_V = 0
    index_H = 0
    df_SP = pd.DataFrame(columns=['station', 't_SP_R', 't_SP_T'])
    # Fill the arrays
    for station, channel, latitude, longitude in zip(stations, channels, latitudes, longitudes):
        channel = channel.split(',')
        # Get vertical component
        filename = 'templates_both/' + family + '/' + station + '_' + channel[-1] + '.pkl'
        trace = pickle.load(open(filename, 'rb'))[0]
        vertical[:, index_V] = trace.data
        # Increment index
        sta_V.append(station)
        index_V = index_V + 1
        # Check if there is an horizontal component
        if len(channel) == 3:
            # Distance
            dx = (pi / 180.0) * a * cos(lat0 * pi / 180.0) / \
                sqrt(1.0 - e * e * sin(lat0 * pi / 180.0) * sin(lat0 * pi / 180.0))
            dy = (3.6 * pi / 648.0) * a * (1.0 - e * e) / \
                ((1.0 - e * e * sin(lat0 * pi / 180.0) * sin(lat0 * pi / 180.0)) ** 1.5)
            x = dx * (longitude - lon0)
            y = dy * (latitude - lat0)
            theta = atan2(y, x)
            # Get horizontal component
            filename = 'templates_both/' + family + '/' + station + '_' + channel[0] + '.pkl'
            trace_EW = pickle.load(open(filename, 'rb'))[0]
            filename = 'templates_both/' + family + '/' + station + '_' + channel[1] + '.pkl'
            trace_NS = pickle.load(open(filename, 'rb'))[0]
            filename = 'templates_both/' + family + '/' + station + '_' + channel[2] + '.pkl'
            trace_UD = pickle.load(open(filename, 'rb'))[0]
            # Rotate horizontal component
            radial[:, index_H] = trace_EW.data * cos(theta) + trace_NS.data * sin(theta)
            transverse[:, index_H] = - trace_EW.data * sin(theta) + trace_NS.data * cos(theta)
            # Compute the S-P travel times
            cc_R = np.correlate(trace_UD.data, radial[:, index_H], mode='full')
            cc_T = np.correlate(trace_UD.data, transverse[:, index_H], mode='full')
            time_R = dt * (nt - np.argmax(np.abs(cc_R)) - 1)
            time_T = dt * (nt - np.argmax(np.abs(cc_T)) - 1)
            df_SP.loc[len(df_SP.index)] = [station, time_R, time_T]
            # Increment index
            sta_H.append(station)
            index_H = index_H + 1
            #
    # Compute the P-P travel times
#    df_P = pd.DataFrame(columns=['station1', 'station2', 't_PP'])
#    for (index1, station1) in enumerate(sta_V):
#        for (index2, station2) in enumerate(sta_V):
#            if index1 < index2:
#                cc = np.correlate(vertical[:, index1], vertical[:, index2], mode='full')
#                time = dt * (np.argmax(np.abs(cc)) - nt)
#                df_P.loc[len(df_P.index)] = [station1, station2, time]
    # Compute the S-S travel times
#    df_S = pd.DataFrame(columns=['station1', 'station2', 't_SS_R', 't_SS_T'])
#    for (index1, station1) in enumerate(sta_H):
#        for (index2, station2) in enumerate(sta_H):
#            if index1 < index2:
#                cc_R = np.correlate(radial[:, index1], radial[:, index2], mode='full')
#                cc_T = np.correlate(transverse[:, index1], transverse[:, index2], mode='full')
#                time_R = dt * (np.argmax(np.abs(cc_R)) - nt)
#                time_S = dt * (np.argmax(np.abs(cc_T)) - nt)
#                df_S.loc[len(df_S.index)] = [station1, station2, time_R, time_S]
    # Save output files
#    namedir = 'traveltimes/' + family
#    if not os.path.exists(namedir):
#        os.makedirs(namedir)
#    pickle.dump(df_P, open(namedir + '/t_PP.pkl', 'wb'))
#    pickle.dump(df_S, open(namedir + '/t_SS.pkl', 'wb'))
#    pickle.dump(df_SP, open(namedir + '/t_SP.pkl', 'wb'))
    return(cc_R, cc_T)

def compute_plot_PP(family, stations, channels, ibegins, nt, dt, nkeep):
    """
    """
    namedir = 'traveltimes/' + family
    if not os.path.exists(namedir):
        os.makedirs(namedir)

    time = dt * np.arange(0, nt)

    params = {'legend.fontsize': 24, \
              'xtick.labelsize':24, \
              'ytick.labelsize':24}
    pylab.rcParams.update(params)

    df_P = pd.DataFrame(columns=['station1', 'station2', 't_PP'])

    for i in range(0, len(stations) - 1):
        for j in range(i + 1, len(stations)):
            filename1 = 'templates_both/' + family + '/' + stations[i] + '_' + channels[i] + '.pkl'
            filename2 = 'templates_both/' + family + '/' + stations[j] + '_' + channels[j] + '.pkl'
            trace1 = pickle.load(open(filename1, 'rb'))[0]
            trace2 = pickle.load(open(filename2, 'rb'))[0]
            data1 = trace1.data[ibegins[i] : ibegins[i] + nkeep + 1]
            data2 = trace2.data[ibegins[j] : ibegins[j] + nkeep + 1]
            cc = np.correlate(data1, data2, mode='full')
            time_cc = dt * np.arange(- nkeep, nkeep + 1)
            tlag = dt * (ibegins[j] - ibegins[i] + nkeep - np.argmax(cc))
            df_P.loc[len(df_P.index)] = [stations[i], stations[j], tlag]

            # Figure
            plt.figure(1, figsize=(8, 15))
            plt.subplot2grid((3, 1), (2, 0))
            plt.plot(time, trace1.data)
            plt.axvline(time[ibegins[i]], color='grey', linewidth=2)
            plt.axvline(time[ibegins[i] + nkeep + 1], color='grey', linewidth=2)
            plt.xlabel('Time (s)', fontsize=24)
            plt.ylabel(stations[i] + ' - ' + channels[i], fontsize=24)
            plt.subplot2grid((3, 1), (1, 0))
            plt.plot(time, trace2.data)
            plt.axvline(time[ibegins[j]], color='grey', linewidth=2)
            plt.axvline(time[ibegins[j] + nkeep + 1], color='grey', linewidth=2)
            plt.xlabel('Time (s)', fontsize=24)
            plt.ylabel(stations[j] + ' - ' + channels[j], fontsize=24)
            plt.subplot2grid((3, 1), (0, 0))
            plt.plot(time_cc, cc)
            plt.axvline(time_cc[np.argmax(cc)], color='grey', linewidth=2)
            plt.xlabel('Time (s)', fontsize=24)
            plt.ylabel('Cross correlation', fontsize=24)
            plt.suptitle('{} - {} - Time lag = {:.2f} seconds'.format(stations[i], stations[j], tlag), fontsize=24)
            plt.tight_layout()
            plt.savefig(namedir + '/' + stations[i] + '_' + stations[j] + '_PP.eps', format='eps')
            plt.close(1)

    pickle.dump(df_P, open(namedir + '/t_PP.pkl', 'wb'))

def compute_plot_SP(family, stations, channels, latitudes, longitudes, ibegins_P, ibegins_S, nt, dt, nkeep, lat0, lon0):
    """
    """
    # To compute distances
    a = 6378.136
    e = 0.006694470

    namedir = 'traveltimes/' + family
    if not os.path.exists(namedir):
        os.makedirs(namedir)

    time = dt * np.arange(0, nt)

    params = {'legend.fontsize': 24, \
              'xtick.labelsize':24, \
              'ytick.labelsize':24}
    pylab.rcParams.update(params)

    df_SP = pd.DataFrame(columns=['station', 't_SP'])

    for i in range(0, len(stations)):
        channel = channels[i].split(',')
        filename_EW = 'templates_both/' + family + '/' + stations[i] + '_' + channel[0] + '.pkl'
        filename_NS = 'templates_both/' + family + '/' + stations[i] + '_' + channel[1] + '.pkl'
        filename_UD = 'templates_both/' + family + '/' + stations[i] + '_' + channel[2] + '.pkl'
        trace_EW = pickle.load(open(filename_EW, 'rb'))[0]
        trace_NS = pickle.load(open(filename_NS, 'rb'))[0]
        trace_UD = pickle.load(open(filename_UD, 'rb'))[0]
        # Distance
        dx = (pi / 180.0) * a * cos(lat0 * pi / 180.0) / \
            sqrt(1.0 - e * e * sin(lat0 * pi / 180.0) * sin(lat0 * pi / 180.0))
        dy = (3.6 * pi / 648.0) * a * (1.0 - e * e) / \
            ((1.0 - e * e * sin(lat0 * pi / 180.0) * sin(lat0 * pi / 180.0)) ** 1.5)
        x = dx * (longitudes[i] - lon0)
        y = dy * (latitudes[i] - lat0)
        theta = atan2(y, x)
        # Rotate horizontal component
        radial = trace_EW.data * cos(theta) + trace_NS.data * sin(theta)
        data_radial = radial[ibegins_S[i] : ibegins_S[i] + nkeep + 1]
        data_UD = trace_UD.data[ibegins_P[i] : ibegins_P[i] + nkeep + 1]
        cc = np.correlate(data_UD, data_radial, mode='full')
        time_cc = dt * np.arange(- nkeep, nkeep + 1)
        tlag = dt * (ibegins_S[i] - ibegins_P[i] + nkeep - np.argmax(cc))
        df_SP.loc[len(df_SP.index)] = [stations[i], tlag]

        # Figure
        plt.figure(1, figsize=(8, 10))
        plt.subplot2grid((3, 1), (2, 0))
        plt.plot(time, trace_UD.data)
        plt.axvline(time[ibegins_P[i]], color='grey', linewidth=2)
        plt.axvline(time[ibegins_P[i] + nkeep + 1], color='grey', linewidth=2)
        plt.xlabel('Time (s)', fontsize=24)
        plt.ylabel(stations[i] + ' - ' + channel[2], fontsize=24)
        plt.subplot2grid((3, 1), (1, 0))
        plt.plot(time, radial)
        plt.axvline(time[ibegins_S[i]], color='grey', linewidth=2)
        plt.axvline(time[ibegins_S[i] + nkeep + 1], color='grey', linewidth=2)
        plt.xlabel('Time (s)', fontsize=24)
        plt.ylabel(stations[i] + ' - Radial', fontsize=24)
        plt.subplot2grid((3, 1), (0, 0))
        plt.plot(time_cc, cc)
        plt.axvline(time_cc[np.argmax(cc)], color='grey', linewidth=2)
        plt.xlabel('Time (s)', fontsize=24)
        plt.ylabel('Cross correlation', fontsize=24)
        plt.suptitle('{} - Time lag = {:.2f} seconds'.format(stations[i], tlag), fontsize=24)
        plt.tight_layout()
        plt.savefig(namedir + '/' + stations[i] + '_SP.eps', format='eps')
        plt.close(1)

    pickle.dump(df_SP, open(namedir + '/t_SP.pkl', 'wb'))

if __name__ == '__main__':

    family = '080401.05.050'
    stations = ['GASB', 'KBN', 'KHBB', 'LRB', 'ME08', 'ME12', 'ME27', \
        'ME28', 'ME37', 'ME39', 'ME41', 'ME57', 'WDC']
    channels = ['HHE,HHN,HHZ', 'SHZ', 'HHE,HHN,HHZ', 'SHZ', 'BHE,BHN,BHZ', \
        'BHE,BHN,BHZ', 'BHE,BHN,BHZ', 'BHE,BHN,BHZ', 'BHE,BHN,BHZ', \
        'BHE,BHN,BHZ', 'BHE,BHN,BHZ', 'BHE,BHN,BHZ', 'HHE,HHN,HHZ']
    latitudes = [39.65471, 39.89237, 40.65990, 40.14323, \
        40.222, 40.104, 40.453, 40.327, 40.285, 40.188202, \
        39.884998, 39.9118, 40.57988]
    longitudes =[-122.71595, -123.19503, -123.21966, -122.55772, \
         -123.305, -122.498, -123.155, -122.471, -123.653999, -123.594299, \
         -123.361, -122.5676, -122.54113]
    nt = 1201
    dt = 0.05
    lat0 = 40.09
    lon0 = -122.87
#    compute_traveltimes(family, stations, channels, latitudes, longitudes, \
#        nt, dt, lat0, lon0)

    stations = ['KBN', 'KHBB', 'LRB', 'ME08', 'ME12', 'ME28', 'ME37', 'ME39', 'ME41', 'ME57', 'WDC']
    channels = ['SHZ', 'HHZ', 'SHZ', 'BHZ', 'BHZ', 'BHZ', 'BHZ', 'BHZ', 'BHZ', 'BHZ', 'HHZ']
    ibegins = [360, 420, 160, 360, 200, 200, 480, 440, 400, 200, 280]
    nkeep = 150
    compute_plot_PP(family, stations, channels, ibegins, nt, dt, nkeep)

    stations = ['KHBB', 'ME08', 'ME12', 'ME28', 'ME37', 'ME39', 'ME41', 'ME57', 'WDC']
    channels = ['HHE,HHN,HHZ', 'BHE,BHN,BHZ', 'BHE,BHN,BHZ', 'BHE,BHN,BHZ', 'BHE,BHN,BHZ', \
        'BHE,BHN,BHZ', 'BHE,BHN,BHZ', 'BHE,BHN,BHZ', 'BHE,BHN,BHZ', 'HHE,HHN,HHZ']
    latitudes = [40.65990, 40.222, 40.104, 40.327, 40.285, 40.188202, 39.884998, 39.9118, 40.57988]
    longitudes =[-123.21966, -123.305, -122.498, -122.471, -123.653999, -123.594299, -123.361, -122.5676, -122.54113]
    ibegins_P = [400, 360, 140, 180, 480, 430, 400, 160, 240]
    ibegins_S = [700, 600, 340, 380, 800, 740, 640, 360, 440]
    nkeep = 200
    compute_plot_SP(family, stations, channels, latitudes, longitudes, ibegins_P, ibegins_S, nt, dt, nkeep, lat0, lon0)
