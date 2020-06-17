# import all need modules...
# Check python version
import sys, argparse, os
ver = sys.version_info # Get Python version
version = ver.major + 0.1*ver.minor
if version < 3.7:
    sys.exit('Python 3.7 or later is required! You are using py{}'.format(version))

# Astropy -----------------------------------------------------
from astropy.table import Table
from astropy.time  import Time
from astropy.io    import fits
from astropy.coordinates import SkyCoord, solar_system, EarthLocation, ICRS
from astropy       import units
# Others  -----------------------------------------------------
from scipy.interpolate import interp1d

import ast
import numpy as np
import os
from os import listdir
from os.path import isfile, join, isdir
import re
import pandas as pd
from itertools import groupby
from operator import itemgetter
import time
import multiprocessing as mp
from functools import partial
import nlopt
from datetime import datetime
# matplotlib  -----------------------------------------------------
import matplotlib
from matplotlib import gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

matplotlib.use('Agg')
import matplotlib.pyplot as plt

font = {'size'   : 6,
        'style'  : 'normal',
        'family' : 'sans-serif'   }

matplotlib.rc('font', **font)

matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['figure.facecolor'] = 'white'

#matplotlib.rcParams['axes.linewidth'] = 2.0
matplotlib.rcParams['axes.labelsize'] = 6

matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

matplotlib.rcParams['xtick.major.top'] = True
matplotlib.rcParams['ytick.major.right'] = True

matplotlib.rcParams['xtick.minor.top'] = True
matplotlib.rcParams['ytick.minor.right'] = True

matplotlib.rcParams['xtick.major.width'] = .6
matplotlib.rcParams['ytick.major.width'] = .6

matplotlib.rcParams['xtick.labelsize'] = 6
matplotlib.rcParams['ytick.labelsize'] = 6

# -------------------------------------------------------------
import warnings
warnings.filterwarnings("ignore")
np.seterr(divide='ignore', invalid='ignore')
