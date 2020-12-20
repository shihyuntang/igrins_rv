# import all need modules...
import sys, argparse, os, ast, re, logging, nlopt
from os      import listdir
from os.path import isfile, join, isdir
# Check python version
#-------------------------------------------------------------------------------
ver     = sys.version_info # Get Python version
version = ver.major + 0.1*ver.minor
if version < 3.7:
    sys.exit(f'Python 3.7 or later is required! You are using py{version}')

# Astropy -----------------------------------------------------
from astropy.table import Table
from astropy.time  import Time
from astropy.io    import fits
from astropy.coordinates import SkyCoord, solar_system, EarthLocation, ICRS
from astropy       import units
# Others  -----------------------------------------------------
from scipy.interpolate import interp1d

import numpy as np
import pandas as pd
import time

from itertools import groupby
import more_itertools as mit
from operator  import itemgetter
from functools import partial
from datetime  import datetime

import multiprocessing as mp

# matplotlib  -----------------------------------------------------
import matplotlib
from matplotlib import gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# font = {'size'   : 6,
#         'style'  : 'normal',
#         'family' : 'sans-serif'   }
#
# matplotlib.rc('font', **font)
#
# matplotlib.rcParams['figure.dpi'] = 300
# matplotlib.rcParams['figure.facecolor'] = 'white'
#
# #matplotlib.rcParams['axes.linewidth'] = 2.0
# matplotlib.rcParams['axes.labelsize'] = 6
#
# matplotlib.rcParams['xtick.direction'] = 'in'
# matplotlib.rcParams['ytick.direction'] = 'in'
#
# matplotlib.rcParams['xtick.major.top'] = True
# matplotlib.rcParams['ytick.major.right'] = True
#
# matplotlib.rcParams['xtick.minor.top'] = True
# matplotlib.rcParams['ytick.minor.right'] = True
#
# matplotlib.rcParams['xtick.major.width'] = .6
# matplotlib.rcParams['ytick.major.width'] = .6
#
# matplotlib.rcParams['xtick.labelsize'] = 6
# matplotlib.rcParams['ytick.labelsize'] = 6

# -------------------------------------------------------------
import warnings
warnings.filterwarnings("ignore")
np.seterr(divide='ignore', invalid='ignore')

# -------------------------------------------------------------
def read_prepdata(args):
        ## Collect relevant file information from Predata files
        if 'igrins' in os.getcwd().split('/')[-1]:
            A0data   = Table.read('./Input/Prepdata/Prepdata_A0_{}.txt'.format(args.targname), format='ascii')
        else:
            A0data   = Table.read('../Input/Prepdata/Prepdata_A0_{}.txt'.format(args.targname), format='ascii')
        A0nights = np.array(A0data['night'],dtype='str')
        ams0     = np.array(A0data['airmass'])
        obs0     = {str(k):str(v) for k,v in zip(A0data['night'],A0data['obs'])}

        if 'igrins' in os.getcwd().split('/')[-1]:
            targdata = Table.read('./Input/Prepdata/Prepdata_targ_{}.txt'.format(args.targname), format='ascii')
        else:
            targdata = Table.read('../Input/Prepdata/Prepdata_targ_{}.txt'.format(args.targname), format='ascii')
        Tnights = np.array(targdata['night'],dtype='str')
        tags0   = np.array(targdata['tag'], dtype='int')
        beams0  = np.array(targdata['beam'],dtype='str')
        mjds0   = np.array(targdata['mjd'],dtype=float)
        bvcs0   = np.array(targdata['bvc'])
        ams     = np.array(targdata['airmass'])

        if 'igrins' in os.getcwd().split('/')[-1]:
            bounddata = Table.read('./Input/UseWv/XRegions_{}_{}.csv'.format(args.WRegion, args.band), format='csv')
        else:
            bounddata = Table.read('../Input/UseWv/XRegions_{}_{}.csv'.format(args.WRegion, args.band), format='csv')
        starts  = np.array(bounddata['start'])
        ends    = np.array(bounddata['end'])
        orders  = np.array(bounddata['order'], dtype=int)
        masks    = np.array(bounddata['masks'])
        xbounddict = {orders[i]:np.array([starts[i],ends[i]]) for i in range(len(starts))}
        maskdict = {orders[i]:masks[i] for i in range(len(starts))}

        # Attribute A and B exposures to right file numbers
        tagsA = {}; tagsB = {}; mjds = {}; bvcs = {};
        night_orig = Tnights[0]; tagsA0 = []; tagsB0 = [];

        nights_unique = np.unique(Tnights)
        for hrt in range(len(nights_unique)):
            jdset = mjds0[(Tnights == nights_unique[hrt])]
            mjds[nights_unique[hrt]] = np.nanmean(jdset)

        for hrt in range(len(Tnights)):
            tag1 = '{:04d}'.format(tags0[hrt])

            bvcs[str(Tnights[hrt])+str(tag1)] = float(bvcs0[hrt])

            if Tnights[hrt] == night_orig:
                if beams0[hrt] == 'A':
                    tagsA0.append(tag1)
                else:
                    tagsB0.append(tag1)
            else:
                tagsA[Tnights[hrt-1]] = tagsA0
                tagsB[Tnights[hrt-1]] = tagsB0
                tagsA0 = []; tagsB0 = [];
                if beams0[hrt] == 'A':
                    tagsA0.append(tag1)
                else:
                    tagsB0.append(tag1)
                night_orig = Tnights[hrt].copy()

        tagsA[Tnights[-1]] = tagsA0
        tagsB[Tnights[-1]] = tagsB0

        nightsFinal = np.array(list(sorted(set(Tnights))))

        obs = np.array([obs0[n[:8]] for n in nightsFinal])

        return xbounddict, maskdict, tagsA, tagsB, mjds, bvcs, nightsFinal, orders, obs
