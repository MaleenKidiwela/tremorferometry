import obspy

import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle

from datetime import datetime
from math import cos, pi, sin, sqrt

def refine_templates(directory, family, stations, tbegin, tend):
    """
    """
    staloc = pd.read_csv('../data/Ducellier/stations_permanent.txt', \
        sep=r'\s{1,}', header=None, engine='python')
    staloc.columns = ['station', 'network', 'channels', 'location', \
        'server', 'latitude', 'longitude', 'time_on', 'time_off', 'dt']

    params = {'xtick.labelsize':20, \
              'ytick.labelsize':20}
    pylab.rcParams.update(params)
    
    # Create directory to store the waveforms
    namedir = directory + '/' + family
    if not os.path.exists(namedir):
        os.makedirs(namedir)

    stations = stations.split(',')

    for station in stations:
        channels = []
        for ir in range(0, len(staloc)):
            if staloc['station'][ir] == station:
                time_on = staloc['time_on'][ir]
                time_off = staloc['time_off'][ir]
                starttime = datetime(int(time_on[0:4]), int(time_on[5:7]), int(time_on[8:10]))
                endtime = datetime(int(time_off[0:4]), int(time_off[5:7]), int(time_off[8:10]))      
                if (endtime >= datetime(2007, 7, 1)) and (starttime <= datetime(2009, 7, 1)):
                    channels = staloc['channels'][ir]
                    channels = channels.split(',')
        for channel in channels:
            filename = 'templates_2007_2009/' + family + '/' + station + '_' + channel + '.pkl'
            data = pickle.load(open(filename, 'rb'))
            stack = data[0]
            stack = stack.slice(stack.stats.starttime + tbegin, stack.stats.starttime + tend)
            dt = stack.stats.delta
            nt = stack.stats.npts
            t = dt * np.arange(0, nt)

            # Plot
            plt.figure(1, figsize=(10, 5))
            plt.plot(t, stack.data, 'k')
            plt.xlim([np.min(t), np.max(t)])
            plt.title('{} - {}'.format(station, channel), fontsize=20)
            plt.xlabel('Time (s)', fontsize=20)
            plt.tight_layout()
            plt.savefig(namedir + '/' + station + '_' + channel + '.eps', format='eps')
            plt.close(1)

            # Save
            stack.stats.station = station
            stack.stats.channel = channel
            data[0] = stack
            savename = namedir + '/' + station + '_' + channel + '.pkl'
            pickle.dump(data, open(savename, 'wb'))

if __name__ == '__main__':

    # Set the parameters
    directory = 'templates_v1'

    families = pd.read_csv('../data/Ducellier/families_permanent.txt', \
        sep=r'\s{1,}', header=None, engine='python')
    families.columns = ['family', 'stations', 'tbegin', 'tend']

    for i in range(0, len(families)):
        family = families['family'][i]
        stations = families['stations'][i]
        tbegin = families['tbegin'][i]
        tend = families['tend'][i]
        refine_templates(directory, family, stations, tbegin, tend)
