"""
Analyse GalNest output
A. Malyali, M. Rivi 2019
"""
from __future__ import absolute_import, unicode_literals
import os
import sys
import argparse
import numpy as np
import pymultinest
import pandas as pd
from sklearn.cluster import MeanShift

parser = argparse.ArgumentParser(description='Modes-analyse')
parser.add_argument('prefix', help='Prefix MultiNest output') 
parser.add_argument('npar', type=int, help='Number of parameters')
parser.add_argument('threshold', type=float, help='Flux threshold [uJy]')
args = parser.parse_args(sys.argv[1:])

def identify_clusters(df_modes):
    """
    Use meanshift algorithm to identify spatial clustering of modes returned by MultiNest.
    :param df_modes:
    :return: input df with each mode assigned a cluster id.
    """
    positions = []
    for x, y in zip(df_modes['mean'].str[0].values, df_modes['mean'].str[1].values):
        positions.append((x, y))

    positions = np.asarray(positions)

    # Break each measurement run into clusters of points by position
    ms = MeanShift(bandwidth=0.0001, bin_seeding=True)
    ms.fit(positions)

    # cluster_centers = ms.cluster_centers_
    labels = ms.labels_
    df_modes['cluster_id'] = labels
    return df_modes

# Fetch stats
a = pymultinest.Analyzer(n_params=args.npar, outputfiles_basename=args.prefix)
s = a.get_stats()

mode_stats = a.get_mode_stats()
df = pd.DataFrame.from_dict(mode_stats['modes'])

df_modes = identify_clusters(df)
df_modes['final'] = df_modes.groupby(['cluster_id'])['local log-evidence'].transform(max) == df_modes['local log-evidence']

l=[]
errl=[]
m=[]
errm=[]
flux=[]
errf=[]
scale=[]
errscale=[]
e1=[]
e2=[]
ev=[]
threshold = args.threshold
for index, row in df_modes.iterrows():
    S = row['mean'][2]*1e+6 #uJy
    R = row['mean'][3]
    if row['final'] and (S/(R*R) > args.threshold):
        l.append(row['mean'][0])
        errl.append(row['sigma'][0])
        m.append(row['mean'][1])
        errm.append(row['sigma'][1])
        flux.append(S)
        errf.append(row['sigma'][2])
        scale.append(R)
        errscale.append(row['sigma'][3])
        e1.append(row['mean'][4])
	e2.append(row['mean'][5])
	ev.append(row['local log-evidence'])

print len(l), " final modes"
results = np.empty((len(l),11))
results[:,0] = np.array(l)
results[:,1] = np.array(errl)
results[:,2] = np.array(m)
results[:,3] = np.array(errm)
results[:,4] = np.array(flux)
results[:,5] = np.array(errf)
results[:,6] = np.array(scale)
results[:,7] = np.array(errscale)
results[:,8] = np.array(e1)
results[:,9] = np.array(e2)
results[:,10] = np.array(ev)
np.savetxt("results.txt",results)
