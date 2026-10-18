"""
Script to read monthly catalogs and put them
into a single file
"""
import numpy as np
import os
import pandas as pd
import pickle

# List of LFE families
templates = np.loadtxt('../Plourde_2015/templates_list.txt', \
    dtype={'names': ('name', 'family', 'lat', 'lon', 'depth', 'eH', \
    'eZ', 'nb'), \
         'formats': ('S13', 'S3', np.float, np.float, np.float, \
    np.float, np.float, np.int)}, \
    skiprows=1)

for i in range(0, np.shape(templates)[0]):

    # Create directory to store the catalog
    namedir = 'catalogs/' + templates[i][0].astype(str)
    if not os.path.exists(namedir):
        os.makedirs(namedir)

    # Read monthly catalogs
    df_2004_01 = pickle.load(open('new/catalog_2004_01/' + templates[i][0].astype(str) + '/catalog_2004_01.pkl', 'rb'))
    df_2004_02 = pickle.load(open('new/catalog_2004_02/' + templates[i][0].astype(str) + '/catalog_2004_02.pkl', 'rb'))
    df_2004_03 = pickle.load(open('new/catalog_2004_03/' + templates[i][0].astype(str) + '/catalog_2004_03.pkl', 'rb'))
    df_2004_04 = pickle.load(open('new/catalog_2004_04/' + templates[i][0].astype(str) + '/catalog_2004_04.pkl', 'rb'))
    df_2004_05 = pickle.load(open('new/catalog_2004_05/' + templates[i][0].astype(str) + '/catalog_2004_05.pkl', 'rb'))
    df_2004_06 = pickle.load(open('new/catalog_2004_06/' + templates[i][0].astype(str) + '/catalog_2004_06.pkl', 'rb'))
    df_2004_07 = pickle.load(open('new/catalog_2004_07/' + templates[i][0].astype(str) + '/catalog_2004_07.pkl', 'rb'))
    df_2004_08 = pickle.load(open('new/catalog_2004_08/' + templates[i][0].astype(str) + '/catalog_2004_08.pkl', 'rb'))
    df_2004_09 = pickle.load(open('new/catalog_2004_09/' + templates[i][0].astype(str) + '/catalog_2004_09.pkl', 'rb'))
    df_2004_10 = pickle.load(open('new/catalog_2004_10/' + templates[i][0].astype(str) + '/catalog_2004_10.pkl', 'rb'))
    df_2004_11 = pickle.load(open('new/catalog_2004_11/' + templates[i][0].astype(str) + '/catalog_2004_11.pkl', 'rb'))
    df_2004_12 = pickle.load(open('new/catalog_2004_12/' + templates[i][0].astype(str) + '/catalog_2004_12.pkl', 'rb'))
    df_2005_01 = pickle.load(open('new/catalog_2005_01/' + templates[i][0].astype(str) + '/catalog_2005_01.pkl', 'rb'))
    df_2005_02 = pickle.load(open('new/catalog_2005_02/' + templates[i][0].astype(str) + '/catalog_2005_02.pkl', 'rb'))
    df_2005_03 = pickle.load(open('new/catalog_2005_03/' + templates[i][0].astype(str) + '/catalog_2005_03.pkl', 'rb'))
    df_2005_04 = pickle.load(open('new/catalog_2005_04/' + templates[i][0].astype(str) + '/catalog_2005_04.pkl', 'rb'))
    df_2005_05 = pickle.load(open('new/catalog_2005_05/' + templates[i][0].astype(str) + '/catalog_2005_05.pkl', 'rb'))
    df_2005_06 = pickle.load(open('new/catalog_2005_06/' + templates[i][0].astype(str) + '/catalog_2005_06.pkl', 'rb'))
    df_2005_07 = pickle.load(open('new/catalog_2005_07/' + templates[i][0].astype(str) + '/catalog_2005_07.pkl', 'rb'))
    df_2005_08 = pickle.load(open('new/catalog_2005_08/' + templates[i][0].astype(str) + '/catalog_2005_08.pkl', 'rb'))
    df_2005_09 = pickle.load(open('new/catalog_2005_09/' + templates[i][0].astype(str) + '/catalog_2005_09.pkl', 'rb'))
    df_2005_10 = pickle.load(open('new/catalog_2005_10/' + templates[i][0].astype(str) + '/catalog_2005_10.pkl', 'rb'))
    df_2005_11 = pickle.load(open('new/catalog_2005_11/' + templates[i][0].astype(str) + '/catalog_2005_11.pkl', 'rb'))
    df_2005_12 = pickle.load(open('new/catalog_2005_12/' + templates[i][0].astype(str) + '/catalog_2005_12.pkl', 'rb'))
    df_2006_01 = pickle.load(open('new/catalog_2006_01/' + templates[i][0].astype(str) + '/catalog_2006_01.pkl', 'rb'))
    df_2006_02 = pickle.load(open('new/catalog_2006_02/' + templates[i][0].astype(str) + '/catalog_2006_02.pkl', 'rb'))
    df_2006_03 = pickle.load(open('new/catalog_2006_03/' + templates[i][0].astype(str) + '/catalog_2006_03.pkl', 'rb'))
    df_2006_04 = pickle.load(open('new/catalog_2006_04/' + templates[i][0].astype(str) + '/catalog_2006_04.pkl', 'rb'))
    df_2006_05 = pickle.load(open('new/catalog_2006_05/' + templates[i][0].astype(str) + '/catalog_2006_05.pkl', 'rb'))
    df_2006_06 = pickle.load(open('new/catalog_2006_06/' + templates[i][0].astype(str) + '/catalog_2006_06.pkl', 'rb'))
    df_2006_07 = pickle.load(open('new/catalog_2006_07/' + templates[i][0].astype(str) + '/catalog_2006_07.pkl', 'rb'))
    df_2006_08 = pickle.load(open('new/catalog_2006_08/' + templates[i][0].astype(str) + '/catalog_2006_08.pkl', 'rb'))
    df_2006_09 = pickle.load(open('new/catalog_2006_09/' + templates[i][0].astype(str) + '/catalog_2006_09.pkl', 'rb'))
    df_2006_10 = pickle.load(open('new/catalog_2006_10/' + templates[i][0].astype(str) + '/catalog_2006_10.pkl', 'rb'))
    df_2006_11 = pickle.load(open('new/catalog_2006_11/' + templates[i][0].astype(str) + '/catalog_2006_11.pkl', 'rb'))
    df_2006_12 = pickle.load(open('new/catalog_2006_12/' + templates[i][0].astype(str) + '/catalog_2006_12.pkl', 'rb'))
    df_2007_01 = pickle.load(open('new/catalog_2007_01/' + templates[i][0].astype(str) + '/catalog_2007_01.pkl', 'rb'))
    df_2007_02 = pickle.load(open('new/catalog_2007_02/' + templates[i][0].astype(str) + '/catalog_2007_02.pkl', 'rb'))
    df_2007_03 = pickle.load(open('new/catalog_2007_03/' + templates[i][0].astype(str) + '/catalog_2007_03.pkl', 'rb'))
    df_2007_04 = pickle.load(open('new/catalog_2007_04/' + templates[i][0].astype(str) + '/catalog_2007_04.pkl', 'rb'))
    df_2007_05 = pickle.load(open('new/catalog_2007_05/' + templates[i][0].astype(str) + '/catalog_2007_05.pkl', 'rb'))
    df_2007_06 = pickle.load(open('new/catalog_2007_06/' + templates[i][0].astype(str) + '/catalog_2007_06.pkl', 'rb'))
    df_2007_07 = pickle.load(open('new/catalog_2007_07/' + templates[i][0].astype(str) + '/catalog_2007_07.pkl', 'rb'))
    df_2007_08 = pickle.load(open('new/catalog_2007_08/' + templates[i][0].astype(str) + '/catalog_2007_08.pkl', 'rb'))
    df_2007_09 = pickle.load(open('new/catalog_2007_09/' + templates[i][0].astype(str) + '/catalog_2007_09.pkl', 'rb'))
    df_2007_10 = pickle.load(open('new/catalog_2007_10/' + templates[i][0].astype(str) + '/catalog_2007_10.pkl', 'rb'))
    df_2007_11 = pickle.load(open('new/catalog_2007_11/' + templates[i][0].astype(str) + '/catalog_2007_11.pkl', 'rb'))
    df_2007_12 = pickle.load(open('new/catalog_2007_12/' + templates[i][0].astype(str) + '/catalog_2007_12.pkl', 'rb'))
    df_2008_01 = pickle.load(open('new/catalog_2008_01/' + templates[i][0].astype(str) + '/catalog_2008_01.pkl', 'rb'))
    df_2008_02 = pickle.load(open('new/catalog_2008_02/' + templates[i][0].astype(str) + '/catalog_2008_02.pkl', 'rb'))
    df_2008_03 = pickle.load(open('new/catalog_2008_03/' + templates[i][0].astype(str) + '/catalog_2008_03.pkl', 'rb'))
    df_2008_04 = pickle.load(open('new/catalog_2008_04/' + templates[i][0].astype(str) + '/catalog_2008_04.pkl', 'rb'))
    df_2008_05 = pickle.load(open('new/catalog_2008_05/' + templates[i][0].astype(str) + '/catalog_2008_05.pkl', 'rb'))
    df_2008_06 = pickle.load(open('new/catalog_2008_06/' + templates[i][0].astype(str) + '/catalog_2008_06.pkl', 'rb'))
    df_2008_07 = pickle.load(open('new/catalog_2008_07/' + templates[i][0].astype(str) + '/catalog_2008_07.pkl', 'rb'))
    df_2008_08 = pickle.load(open('new/catalog_2008_08/' + templates[i][0].astype(str) + '/catalog_2008_08.pkl', 'rb'))
    df_2008_09 = pickle.load(open('new/catalog_2008_09/' + templates[i][0].astype(str) + '/catalog_2008_09.pkl', 'rb'))
    df_2008_10 = pickle.load(open('new/catalog_2008_10/' + templates[i][0].astype(str) + '/catalog_2008_10.pkl', 'rb'))
    df_2008_11 = pickle.load(open('new/catalog_2008_11/' + templates[i][0].astype(str) + '/catalog_2008_11.pkl', 'rb'))
    df_2008_12 = pickle.load(open('new/catalog_2008_12/' + templates[i][0].astype(str) + '/catalog_2008_12.pkl', 'rb'))
    df_2009_01 = pickle.load(open('new/catalog_2009_01/' + templates[i][0].astype(str) + '/catalog_2009_01.pkl', 'rb'))
    df_2009_02 = pickle.load(open('new/catalog_2009_02/' + templates[i][0].astype(str) + '/catalog_2009_02.pkl', 'rb'))
    df_2009_03 = pickle.load(open('new/catalog_2009_03/' + templates[i][0].astype(str) + '/catalog_2009_03.pkl', 'rb'))
    df_2009_04 = pickle.load(open('new/catalog_2009_04/' + templates[i][0].astype(str) + '/catalog_2009_04.pkl', 'rb'))
    df_2009_05 = pickle.load(open('new/catalog_2009_05/' + templates[i][0].astype(str) + '/catalog_2009_05.pkl', 'rb'))
    df_2009_06 = pickle.load(open('new/catalog_2009_06/' + templates[i][0].astype(str) + '/catalog_2009_06.pkl', 'rb'))
    df_2009_07 = pickle.load(open('new/catalog_2009_07/' + templates[i][0].astype(str) + '/catalog_2009_07.pkl', 'rb'))
    df_2009_08 = pickle.load(open('new/catalog_2009_08/' + templates[i][0].astype(str) + '/catalog_2009_08.pkl', 'rb'))
    df_2009_09 = pickle.load(open('new/catalog_2009_09/' + templates[i][0].astype(str) + '/catalog_2009_09.pkl', 'rb'))
    df_2009_10 = pickle.load(open('new/catalog_2009_10/' + templates[i][0].astype(str) + '/catalog_2009_10.pkl', 'rb'))
    df_2009_11 = pickle.load(open('new/catalog_2009_11/' + templates[i][0].astype(str) + '/catalog_2009_11.pkl', 'rb'))
    df_2009_12 = pickle.load(open('new/catalog_2009_12/' + templates[i][0].astype(str) + '/catalog_2009_12.pkl', 'rb'))
    df_2010_01 = pickle.load(open('new/catalog_2010_01/' + templates[i][0].astype(str) + '/catalog_2010_01.pkl', 'rb'))
    df_2010_02 = pickle.load(open('new/catalog_2010_02/' + templates[i][0].astype(str) + '/catalog_2010_02.pkl', 'rb'))
    df_2010_03 = pickle.load(open('new/catalog_2010_03/' + templates[i][0].astype(str) + '/catalog_2010_03.pkl', 'rb'))
    df_2010_04 = pickle.load(open('new/catalog_2010_04/' + templates[i][0].astype(str) + '/catalog_2010_04.pkl', 'rb'))
    df_2010_05 = pickle.load(open('new/catalog_2010_05/' + templates[i][0].astype(str) + '/catalog_2010_05.pkl', 'rb'))
    df_2010_06 = pickle.load(open('new/catalog_2010_06/' + templates[i][0].astype(str) + '/catalog_2010_06.pkl', 'rb'))
    df_2010_07 = pickle.load(open('new/catalog_2010_07/' + templates[i][0].astype(str) + '/catalog_2010_07.pkl', 'rb'))
    df_2010_08 = pickle.load(open('new/catalog_2010_08/' + templates[i][0].astype(str) + '/catalog_2010_08.pkl', 'rb'))
    df_2010_09 = pickle.load(open('new/catalog_2010_09/' + templates[i][0].astype(str) + '/catalog_2010_09.pkl', 'rb'))
    df_2010_10 = pickle.load(open('new/catalog_2010_10/' + templates[i][0].astype(str) + '/catalog_2010_10.pkl', 'rb'))
    df_2010_11 = pickle.load(open('new/catalog_2010_11/' + templates[i][0].astype(str) + '/catalog_2010_11.pkl', 'rb'))
    df_2010_12 = pickle.load(open('new/catalog_2010_12/' + templates[i][0].astype(str) + '/catalog_2010_12.pkl', 'rb'))
    df_2011_01 = pickle.load(open('new/catalog_2011_01/' + templates[i][0].astype(str) + '/catalog_2011_01.pkl', 'rb'))
    df_2011_02 = pickle.load(open('new/catalog_2011_02/' + templates[i][0].astype(str) + '/catalog_2011_02.pkl', 'rb'))
    df_2011_03 = pickle.load(open('new/catalog_2011_03/' + templates[i][0].astype(str) + '/catalog_2011_03.pkl', 'rb'))
    df_2011_04 = pickle.load(open('new/catalog_2011_04/' + templates[i][0].astype(str) + '/catalog_2011_04.pkl', 'rb'))
    df_2011_05 = pickle.load(open('new/catalog_2011_05/' + templates[i][0].astype(str) + '/catalog_2011_05.pkl', 'rb'))
    df_2011_06 = pickle.load(open('new/catalog_2011_06/' + templates[i][0].astype(str) + '/catalog_2011_06.pkl', 'rb'))
    df_2011_07 = pickle.load(open('new/catalog_2011_07/' + templates[i][0].astype(str) + '/catalog_2011_07.pkl', 'rb'))
    df_2011_08 = pickle.load(open('new/catalog_2011_08/' + templates[i][0].astype(str) + '/catalog_2011_08.pkl', 'rb'))
    df_2011_09 = pickle.load(open('new/catalog_2011_09/' + templates[i][0].astype(str) + '/catalog_2011_09.pkl', 'rb'))
    df_2011_10 = pickle.load(open('new/catalog_2011_10/' + templates[i][0].astype(str) + '/catalog_2011_10.pkl', 'rb'))
    df_2011_11 = pickle.load(open('new/catalog_2011_11/' + templates[i][0].astype(str) + '/catalog_2011_11.pkl', 'rb'))
    df_2011_12 = pickle.load(open('new/catalog_2011_12/' + templates[i][0].astype(str) + '/catalog_2011_12.pkl', 'rb'))

    # Concatenate catalogs
    df_2004_2011 = pd.concat([ \
        df_2004_01, df_2004_02, df_2004_03, df_2004_04, df_2004_05, df_2004_06, \
        df_2004_07, df_2004_08, df_2004_09, df_2004_10, df_2004_11, df_2004_12, \
        df_2005_01, df_2005_02, df_2005_03, df_2005_04, df_2005_05, df_2005_06, \
        df_2005_07, df_2005_08, df_2005_09, df_2005_10, df_2005_11, df_2005_12, \
        df_2006_01, df_2006_02, df_2006_03, df_2006_04, df_2006_05, df_2006_06, \
        df_2006_07, df_2006_08, df_2006_09, df_2006_10, df_2006_11, df_2006_12, \
        df_2007_01, df_2007_02, df_2007_03, df_2007_04, df_2007_05, df_2007_06, \
        df_2007_07, df_2007_08, df_2007_09, df_2007_10, df_2007_11, df_2007_12, \
        df_2008_01, df_2008_02, df_2008_03, df_2008_04, df_2008_05, df_2008_06, \
        df_2008_07, df_2008_08, df_2008_09, df_2008_10, df_2008_11, df_2008_12, \
        df_2009_01, df_2009_02, df_2009_03, df_2009_04, df_2009_05, df_2009_06, \
        df_2009_07, df_2009_08, df_2009_09, df_2009_10, df_2009_11, df_2009_12, \
        df_2010_01, df_2010_02, df_2010_03, df_2010_04, df_2010_05, df_2010_06, \
        df_2010_07, df_2010_08, df_2010_09, df_2010_10, df_2010_11, df_2010_12, \
        df_2011_01, df_2011_02, df_2011_03, df_2011_04, df_2011_05, df_2011_06, \
        df_2011_07, df_2011_08, df_2011_09, df_2011_10, df_2011_11, df_2011_12 \
        ], ignore_index=True)

    # Drop duplicates
    df_2004_2011.drop_duplicates(inplace=True, ignore_index=True)
    
    # Write catalog into file
    namefile = namedir + '/catalog_2004_2011.pkl'
    pickle.dump(df_2004_2011, open(namefile, 'wb'))
#    tfile = open(namedir + '/catalog_2007_2011.txt', 'w')
#    tfile.write(df_2007_2011.to_string())
#    tfile.close()

# To merge catalogs from two different instances
#for i in range(0, np.shape(templates)[0]):
#    df1 = pickle.load(open('catalog_2009_04_1/' + templates[i][0].astype(str) + '/catalog.pkl', 'rb'))
#    df2 = pickle.load(open('catalog_2009_04/' + templates[i][0].astype(str) + '/catalog_2009_04.pkl', 'rb'))
#    df = pd.concat([df1, df2])
#    df.reset_index(drop=True, inplace=True)
#    pickle.dump(df, open('catalog_2009_04/' + templates[i][0].astype(str) + '/catalog_2009_04.pkl', 'wb'))
