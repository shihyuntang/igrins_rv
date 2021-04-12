# import all need modules...
import sys, argparse, os, ast, re, logging, nlopt
import logging.handlers
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
from functools import partial, wraps
from datetime  import datetime

import multiprocessing as mp
from pqdm.processes import pqdm
# matplotlib  -----------------------------------------------------
import matplotlib
from matplotlib import gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# -------------------------------------------------------------
import warnings
warnings.filterwarnings("ignore")
np.seterr(divide='ignore', invalid='ignore')

# -------------------------------------------------------------

def suppress_stdout(f, *args, **kwargs):
    """
    A simple decorator to suppress function print outputs.
    Borrowed from the lightkurve pkg @ https://github.com/lightkurve/lightkurve
    """

    @wraps(f)
    def wrapper(*args, **kwargs):
        # redirect output to `null`
        with open(os.devnull, "w") as devnull:
            old_out = sys.stdout
            sys.stdout = devnull
            try:
                return f(*args, **kwargs)
            # restore to default
            finally:
                sys.stdout = old_out

    return wrapper

# -------------------------------------------------------------

def log_warning_id(file, start_t):
    """
    Under silent mode, check the .log file to see if any logger.warning.
    """
    file1 = open(file, 'r')
    Lines = file1.readlines()
    loop_range = np.arange(len(Lines)-1, -1, -1)

    # find lines (lidx) logging for this run
    for lidx in loop_range:
        line_str = Lines[lidx]
        try:
            int(line_str[:4])
        except ValueError:
            continue
        date_str = line_str[:19] # extract the date, e.g., '2021-04-11 08:29:50'
        datetemp = datetime.strptime(date_str, '%Y-%m-%d %H:%M:%S')

        if start_t > datetemp:
            start_lidx = lidx
            break
        if lidx == loop_range[-1]: # if loop to the first row
            start_lidx = lidx

    this_run = Lines[start_lidx:]

    for i in this_run:
        if 'WARNING' in i:
            return True
    return False


# -------------------------------------------------------------
def read_prepdata(args):
        '''

        Collect relevant file information from Predata files
        
        Inputs:
        args  : Information specified by user at command line 

        Outputs:
        xbounddict  : Dictionary of pixel ranges to be analyzed, referenced by order
        maskdict    : Dictionary of pixel ranges to be masked from analysis, referenced by order
        tagsA       : Dictionary of A frame file numbers, referenced by night
        tagsB       : Ditto but for B frames
        jds         : Julian Dates
        bvcs        : Barycentric velocity corrections
        nightsFinal : Dates of observations in YYYYMMDD
        orders      : Echelle orders, as characterized by file index (as opposed to m number; for conversion between the two, see Stahl et al. 2021) 
        obs         : Dictionary of observatory corresponding to observation, referencedby night
        '''
        
        ## 
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
        jds0   = np.array(targdata['jd'],dtype=float)
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
        tagsA = {}; tagsB = {}; jds = {}; bvcs = {};
        night_orig = Tnights[0]; tagsA0 = []; tagsB0 = [];

        nights_unique = np.unique(Tnights)
        for hrt in range(len(nights_unique)):
            jdset = jds0[(Tnights == nights_unique[hrt])]
            jds[nights_unique[hrt]] = np.nanmean(jdset)

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

        return xbounddict, maskdict, tagsA, tagsB, jds, bvcs, nightsFinal, orders, obs

# -------------------------------------------------------------
#---- logger init ---

# def listener_configurer():
#     root = logging.getLogger()
#     h = logging.handlers.RotatingFileHandler('mptest.log', 'a', 300, 10)
#     f = logging.Formatter('%(asctime)s: %(module)s.py: %(levelname)s--> %(message)s')
#     h.setFormatter(f)
#     root.addHandler(h)
#
# # This is the listener process top-level loop: wait for logging events
# # (LogRecords)on the queue and handle them, quit when you get a None for a
# # LogRecord.
# def listener_process(queue, configurer):
#     configurer()
#     while True:
#         try:
#             record = queue.get()
#             if record is None:  # We send this as a sentinel to tell the listener to quit.
#                 break
#             logger = logging.getLogger(record.name)
#             logger.handle(record)  # No level or filter logic applied - just do it!
#         except Exception:
#             import sys, traceback
#             print('Whoops! Problem:', file=sys.stderr)
#             traceback.print_exc(file=sys.stderr)
#
#
# # The worker configuration is done at the start of the worker process run.
# # Note that on Windows you can't rely on fork semantics, so each process
# # will run the logging configuration code when it starts.
# def worker_configurer(queue):
#     h = logging.handlers.QueueHandler(queue)  # Just the one handler needed
#     root = logging.getLogger()
#     formatter = logging.Formatter('%(asctime)s: %(module)s.py: %(levelname)s--> %(message)s')
#     root.addHandler(h)
#     # send all messages, for demo; no other level or filter logic applied.
#     root.setLevel(logging.DEBUG)
#
# # This is the worker process top-level loop, which just logs ten events with
# # random intervening delays before terminating.
# # The print messages are just so you know it's doing something!
# def worker_process(queue, configurer):
#     configurer(queue)
#     name = multiprocessing.current_process().name
#     print('Worker started: %s' % name)
#     for i in range(10):
#         time.sleep(random())
#         logger = logging.getLogger(choice(LOGGERS))
#         level = choice(LEVELS)
#         message = choice(MESSAGES)
#         logger.log(level, message)
#     print('Worker finished: %s' % name)
