"""
handeling all the argparse setup for all steps
"""

import argparse
import multiprocessing as mp

_epilog = "Contact authors: Asa Stahl (asa.stahl@rice.edu);\
                Shih-Yun Tang (sytang@lowell.edu)"
_igrins_version = '1.5.0 alpha'

def _argparse_step0_s():
    """Take care of all the argparse stuff."""
    parser = argparse.ArgumentParser(
        prog        = 'IGRINS Spectra Radial Velocity Pipeline - Step 0_s',
        description = '''
        The .py do the same as the Step0, but split the AB beam (nodding) \n
        into sets that user assign in the .recipes files. \n
        See more detail on 
        https://github.com/shihyuntang/igrins_rv/wiki/Setup:-PrepData-Files#muti-rvs-per-night-version-of-step0py-night_splitstep0_spy
        This step collects and organizes all relevant information on target \n
        observations and associated telluric standard observations, for ease \n
        of use in later steps. It requires that your observations be listed \n
        in IGRINS_RV_MASTERLOG.csv, which comes with this package in the \n
        /Engine folder. If your target is not listed, you must construct \n
        your own PrepData files. See ReadMe for more details. \n
        ''',
        epilog = _epilog)
    parser.add_argument("targname", action="store",
        help="Enter your *target name", 
        type=str)
    parser.add_argument("-HorK", dest="band", action="store",
        help="Which band to process? H or K?",
        type=str, default='K')

    parser.add_argument("-coord", dest="coord", action="store",
        help="Optional [-XX.xx,-XX.xx] deg, GaiaDR2 coordinates at J2015.5. \
                If give, will calculate BVC base on this info.",
        type=str, default='')
    parser.add_argument("-pm", dest="pm", action="store",
        help="Optional [-XX.xx,-XX.xx] [mas/yr], GaiaDR2 proper motion. \
                If give, will calculate BVC base on this info.",
        type=str, default='')

    parser.add_argument(
        '--version', action='version', 
        version='%(prog)s {}'.format(_igrins_version))

    return parser.parse_args()



def _argparse_step0():
    """Take care of all the argparse stuff."""
    parser = argparse.ArgumentParser(
        prog        = 'IGRINS Spectra Radial Velocity Pipeline - Step 0',
        description = '''
        This step collects and organizes all relevant information on target \n
        observations and associated telluric standard observations, for ease \n
        of use in later steps. It requires that your observations be listed \n
        in IGRINS_RV_MASTERLOG.csv, which comes with this package in the \n
        /Engine folder. If your target is not listed, you must construct \n
        your own PrepData files. See ReadMe for more details. \n
        ''',
        epilog = _epilog)
    parser.add_argument("targname", action="store",
        help="Enter your target name, no space",
        type=str)
    parser.add_argument("-HorK", dest="band", action="store",
        help="Which band to process? H or K?. Default = K",
        type=str, default='K')
    parser.add_argument("-AM", dest="AM_cut", action="store",
        help="AirMass difference allowed between TAR and STD (A0) stars. \
                Default X = 0.3 ",
        type=str, default='0.3')

    parser.add_argument("-coord", dest="coord", action="store",
        help="Optional [-XX.xx,-XX.xx] deg, GaiaDR2 coordinates at J2015.5. \
                If give, will calculate BVC base on this info.",
        type=str, default='')
    parser.add_argument("-pm", dest="pm", action="store",
        help="Optional [-XX.xx,-XX.xx] [mas/yr], GaiaDR2 proper motion. \
                If give, will calculate BVC base on this info.",
        type=str, default='')
    
    parser.add_argument(
        '--version', action='version', 
        version='%(prog)s {}'.format(_igrins_version))
    
    return parser.parse_args()


def _argparse_step1():
    """Take care of all the argparse stuff."""
    parser = argparse.ArgumentParser(
        prog        = 'IGRINS Spectra Radial Velocity Pipeline - Step 1',
        description = '''
        This step 1) defines the wavelength regions to be analyzed based on \n
        user specification 2) generates a synthetic, high-resolution telluric \n
        template for use in later model fits on a night by night basis. \n
        Note that only target star observations will have their fits limited \n
        to the wavelength regions specified. \n
        For A0 observations, only the orders specified will be analyzed, \n
        but each order will be fit as far as there is significant telluric absorption.
        ''',
        epilog = _epilog)
    parser.add_argument("targname", action="store",
        help="Enter your target name, no space",
        type=str)
    parser.add_argument("-HorK", dest="band", action="store",
        help="Which band to process? H or K?. Default = K",
        type=str,   default='K')
    parser.add_argument("-Wr", dest="WRegion", action="store",
        help="Which list of wavelength regions file (./Input/UseWv/WaveRegions_X) \
                to use? Defaults to those chosen by IGRINS RV team, -Wr 1",
        type=int,   default=int(1))

    parser.add_argument("-SN", dest="SN_cut", action="store",
        help="Spectrum S/N quality cut. Spectra with median S/N below this \
                will not be analyzed. Default = 50 ",
        type=str,   default='50')

    parser.add_argument('-c', dest="Nthreads", action="store",
        help="Number of cpu (threads) to use, default is 1/2 of avalible \
                ones (you have %i cpus (threads) avaliable)"%(mp.cpu_count()),
        type=int,   default=int(mp.cpu_count()//2) )
    parser.add_argument('-plot', dest="plotfigs", action="store_true",
        help="If set, will generate plots of A0 fitting results under \
                ./Output/A0Fits/*target/fig/")

    parser.add_argument('-n_use', dest="nights_use", action="store",
        help="If you don't want to process all nights under the \
                ./Input/*target/ folder, specify an array of night you wish \
                to process here. e.g., [20181111,20181112]",
        type=str,   default='')
    parser.add_argument('-DeBug', dest="debug", action="store_true",
        help="If set, DeBug logging will be output, as well as (lots of) \
                extra plots.")
    
    parser.add_argument(
        '--version', action='version', 
        version='%(prog)s {}'.format(_igrins_version))
    
    return parser.parse_args()


def _argparse_step2():
    """Take care of all the argparse stuff."""
    parser = argparse.ArgumentParser(
        prog        = 'IGRINS Spectra Radial Velocity Pipeline - Step 2',
        description = '''
        Required if the average RV of the target star is unknown to >5 km/s precision. \n
        Performs an abbreviated analysis of the target star observations in \n
        order to converge to coarsely accurate RVs, which will be used as \n
        starting points for the more precise analysis in the next step; \n
        simultaneously does the same for target star's vsini. \n
        Only the single most precise wavelength region is used, and all \n
        separate observations for a given exposure are combined into one \n
        higher S/N spectrum before being fit.
        ''',
        epilog = _epilog)
    parser.add_argument("targname", action="store",
        help="Enter your *target name, no space", type=str)
    parser.add_argument("-HorK", dest="band", action="store",
        help="Which band to process? H or K?. Default = K",
        type=str, default='K')
    parser.add_argument("-Wr", dest="WRegion", action="store",
        help="Which list of wavelength regions file (./Input/UseWv/WaveRegions_X) \
                to use? Defaults to those chosen by IGRINS RV team, -Wr 1",
        type=int, default=int(1))
    parser.add_argument("-l_use", dest="label_use", action="store",
        help="Specify ORDER used. Default is the first in WRegion list",
        type=int, default=int(0))
    
    parser.add_argument("-SN", dest="SN_cut", action="store",
        help="Spectrum S/N quality cut. Spectra with median S/N below this \
                will not be analyzed. Default = 50 ",
        type=str, default='50')
    parser.add_argument('-i', dest="initvsini", action="store",
        help="Initial vsini (float, km/s)",
        type=str,   default='' )
    parser.add_argument('-v', dest="vsinivary", action="store",
        help="Range of allowed vsini variation during optimization (float, \
                if set to zero vsini will be held constant), \
                default = 5.0 km/s",
        type=str, default='5' )
    parser.add_argument('-g', dest="guesses", action="store",
        help="Initial guess for RV (int or float, km/s) for all nights. \
                Use -gX instead if you want to reference an Initguesser_results \
                file from a previous run of this step, which will have a \
                different initial guess for each night",
        type=str, default='')
    parser.add_argument('-gX', dest="guessesX", action="store",
        help="The number, X, that refers to the ./*targname/Initguesser_results_X \
                file you wish to use for initial RV guesses",
        type=str, default='')
    
    parser.add_argument('-t', dest="template", action="store",
        help="Stellar template. Pick from 'synthetic', 'PHOENIX', or \
                'livingston'. Default = 'synthetic'",
        type=str, default='synthetic' )
    parser.add_argument('-temp', dest="temperature", action="store",
        help="The synthetic template temperature used, e.g., 5000",
        type=str,   default='' )
    parser.add_argument('-logg', dest="logg", action="store",
        help="The synthetic template logg used, e.g., 4.5",
        type=str, default='' )
    
    parser.add_argument('-B',      dest="B",           action="store",
                        help="The synthetic template B used in kG, e.g., 2.5",
                        type=str,   default='0' )
    parser.add_argument('-c', dest="Nthreads", action="store",
        help="Number of cpu (threads) to use, default is 1/2 of available \
                ones (you have %i cpus (threads) avaliable)"%(mp.cpu_count()),
        type=int,   default=int(mp.cpu_count()//2) )
    parser.add_argument('-plot', dest="plotfigs", action="store_true",
        help="If set, will generate plots of the basic fitting results under \
                ./Output/*targname_*band/figs/main_step2_*band_*runnumber")
    parser.add_argument('-n_use', dest="nights_use", action="store",
        help="If you don't want to process all nights under the ./Input/*target/ \
                folder, specify an array of night you wish to process here. \
                e.g., [20181111,20181112]",
        type=str,   default='')
    parser.add_argument('-DeBug', dest="debug", action="store_true",
        help="If set, DeBug logging will be output, as well as (lots of) \
                extra plots.")
    parser.add_argument('-sk_check', dest="skip", action="store_true",
        help="If set, will skip the input parameters check. Handy when running \
                multiple targets line by line")
    
    parser.add_argument(
        '--version', action='version', 
        version='%(prog)s {}'.format(_igrins_version))
    
    return parser.parse_args()


def _argparse_step3():
    """Take care of all the argparse stuff."""
    parser = argparse.ArgumentParser(
        prog        = 'IGRINS Spectra Radial Velocity Pipeline - Step 3',
        description = '''
        Performs a full analysis of each target star observation to produce \n
        accurate and precise RVs. All the wavelength regions defined in Step 1 \n
        are used, and the code analyzes each observation that is part of a given\n
        exposure separately. \n
        Unless the target vsini is already known to high accuracy, an initial \n
        run of Step 3 in which \vsini is allowed to vary is required. \n
        This provides an estimate of vsini that can then be plugged into the \n
        code as a fixed value in the second run of Step 3. \n
        If the user seeks the best possible RV uncertainty estimates, or if \n
        their target star has a relatively high \vsini ($>$ 10 \kms), they \n
        must run Step 3 once with \vsini held fixed at its estimated value \n
        and once with \vsini held fixed at this value plus or minus one sigma. \n
        The minor differences in the RVs of the two runs (as low as $<$1 \ms \n
        and as high as 7 \ms) can then be incorporated into the final uncertainties. \n
        If \vsini is already well-known, it is not necessary to run Step 3 more \n
        than once, as the code fully converges to the final RVs (within uncertainty) \n
        through just one run.
        ''',
        epilog = _epilog)
    parser.add_argument("targname", action="store",
        help="Enter your *target name", type=str)
    parser.add_argument("-mode", dest="mode", action="store",
        help="RV standard star (STD) or a normal target (TAR)?",
        type=str, default='')
    parser.add_argument("-HorK", dest="band", action="store",
        help="Which band to process? H or K?. Default = K",
        type=str, default='K')
    parser.add_argument("-Wr", dest="WRegion", action="store",
        help="Which list of wavelength regions file (./Input/UseWv/WaveRegions_X) \
                to use? Defaults to those chosen by IGRINS RV team, -Wr 1",
        type=int, default=int(1))
    
    parser.add_argument("-SN", dest="SN_cut", action="store",
        help="Spectrum S/N quality cut. Spectra with median S/N below this \
                will not be analyzed. Default = 50 ",
        type=str, default='50')
    parser.add_argument("-nAB", dest="nAB", action="store",
        help="Minium number of separate A/B exposures within a set for a \
                given observation (ensures accuracy of uncertainly estimates). \
                Default = 2 for STD, 3 for TAR",
        type=str, default='')
    
    parser.add_argument('-i', dest="initvsini", action="store",
        help="Initial vsini (float, km/s). If no literature value known, use \
                the value given by Step 2",
        type=str, default='' )
    parser.add_argument('-v', dest="vsinivary", action="store",
        help="Range of allowed vsini variation during optimization, \
                default = 5.0 km/s. Should be set to 0 for final run.",
        type=str, default='5.0' )
    parser.add_argument('-g', dest="guesses", action="store",
        help="For STD star. Initial RV guess for all nights. Given by \
                Step 2 results (float, km/s)",
        type=str, default='' )
    parser.add_argument('-gS', dest="guesses_source", action="store",
        help="For TAR star. Source for list of initial RV guesses. \
                'init' = Initguesser_results_X = past Step 2 result OR \
                'rvre' = RV_results_X = past Step 3 result",
        type=str, default='')
    parser.add_argument('-gX', dest="guessesX", action="store",
        help="For TAR star. The number, X, under ./*targname/Initguesser_results_X \
                or ./*targname/RV_results_X, that you wish to use. \
                Prefix determined by -gS",
        type=str, default='')
    
    parser.add_argument('-t', dest="template", action="store",
        help="Stellar template. Pick from 'synthetic', 'PHOENIX', or \
                'livingston'. Default = 'synthetic'",
        type=str, default='synthetic' )
    parser.add_argument('-temp', dest="temperature", action="store",
        help="The synthetic template temperature used, e.g., 5000",
        type=str, default='' )
    parser.add_argument('-logg', dest="logg", action="store",
        help="The synthetic template logg used, e.g., 4.5",
        type=str, default='' )
   
    parser.add_argument('-B',      dest="B",           action="store",
                        help="The synthetic template B used in kG, e.g., 2.5",
                        type=str,   default='0' )
    
    parser.add_argument('-c', dest="Nthreads", action="store",
        help="Number of cpu (threads) to use, default is 1/2 of available \
                ones (you have %i cpus (threads) avaliable)"%(mp.cpu_count()),
        type=int,   default=int(mp.cpu_count()//2) )
    parser.add_argument('-plot', dest="plotfigs", action="store_true",
        help="If set, will generate plots of the fitting results under \
                ./Output/*targname_*band/figs/main_step3_*band_*runnumber")
    parser.add_argument('-n_use', dest="nights_use", action="store",
        help="If you don't want to process all nights under the ./Input/*target/ \
                folder, specify an array of night you wish to process here. \
                e.g., [20181111,20181112]",
        type=str, default='')
    parser.add_argument('-DeBug', dest="debug", action="store_true",
        help="If set, DeBug logging will be output, as well as (lots of) extra plots.")
    parser.add_argument('-sk_check', dest="skip", action="store_true",
        help="If set, will skip the input parameters check. Handy when running \
                multiple targets line by line")
    
    parser.add_argument(
        '--version', action='version', 
        version='%(prog)s {}'.format(_igrins_version))
    
    return parser.parse_args()


def _argparse_step4():
    """Take care of all the argparse stuff."""
    parser = argparse.ArgumentParser(
        prog        = 'IGRINS Spectra Radial Velocity Pipeline - Step 4',
        description = '''
        Updates RV uncertainty estimates to take into account uncertainty in vsini. \n
        Takes two runs of Step 3, one with vsini held fixed at the best guess \n
        value and one with vsini held fixed at the best guess value plus or \n
        minus one sigma, and uses the difference between the two to produce \n
        updated RVs and uncertainties.

        For the most part, the uncertainties should change little (~1 m/s), \n
        but for high vsini (>~ 15 km/s) objects, it may increase by ~5-7 m/s or so.
        ''',
        epilog = _epilog)
    parser.add_argument("targname",                          action="store",
        help="Enter your *target name",            type=str)
    
    parser.add_argument("-run1",    dest="run1",             action="store",
        help="First step3 run that will be used, the one with vsini held fixed \
                AT THE BEST GUESS. Takes the string that suffixes \
                'RV_results_', e.g. for 'RV_results_1', you would set this to '1'.",
        type=str,   default='')
    parser.add_argument("-run2",    dest="run2",             action="store",
        help="Second step3 run that will be used, the one with vsini held fixed \
                at the best guess PLUS OR MINUS SIGMA.",
        type=str,   default='')
    
    parser.add_argument("-HorK",    dest="band",             action="store",
        help="Which band to process? H or K?. Default = K",
        type=str,   default='K')
    
    parser.add_argument('-sk_check', dest="skip",           action="store_true",
        help="If set, will skip the input parameters check. Handy when running \
                multiple targets line by line")
    
    
def _argparse_step5():
    """Take care of all the argparse stuff."""
    parser = argparse.ArgumentParser(
        prog        = 'IGRINS Spectra Radial Velocity Pipeline - Step 5',
        description = '''
        Outputs telluric and continuum-corrected stellar residuals based on a \n
        run of Step 3. Combines A's and B's and saves results, but also saves \n
        results on an individual tag-by-tag basis. The latter can be used \n
        to compute bisectors if the corresponding flag is set. If you've \n
        already produced the stellar residuals, you can set a flag to skip \n
        the first part of Step 5 and go straight to calculating bisectors.
        ''',
        epilog = _epilog)
    parser.add_argument("targname", action="store",
        help="Enter your *target name", type=str)
    parser.add_argument("-HorK", dest="band", action="store",
        help="Which band to process? H or K?. Default = K",
        type=str, default='K')
    parser.add_argument("-Wr", dest="WRegion", action="store",
        help="Which list of wavelength regions file (./Input/UseWv/WaveRegions_X) \
                to use? Defaults to those chosen by IGRINS RV team, -Wr 1",
        type=int, default=int(1))

    parser.add_argument("-mode", dest="mode", action="store",
        help="RV standard star (STD) or a normal target (TAR)?",
        type=str, default='')
    parser.add_argument("-run",    dest="run",             action="store",
        help="Takes the string that suffixes \
                'RV_results_', e.g. for 'RV_results_1', you would set this to '1'.",
        type=str,   default='')

    parser.add_argument('-t', dest="template", action="store",
        help="Stellar template. Pick from 'synthetic', 'PHOENIX', or \
                'livingston'. Default = 'synthetic'",
        type=str, default='synthetic' )
    parser.add_argument('-temp', dest="temperature", action="store",
        help="The synthetic template temperature used, e.g., 5000",
        type=str, default='' )
    parser.add_argument('-logg', dest="logg", action="store",
        help="The synthetic template logg used, e.g., 4.5",
        type=str, default='' )

    parser.add_argument('-B',      dest="B",           action="store",
                        help="The synthetic template B used in kG, e.g., 2.5",
                        type=str,   default='0' )

    parser.add_argument('-c', dest="Nthreads", action="store",
        help="Number of cpu (threads) to use, default is 1/2 of available \
                ones (you have %i cpus (threads) avaliable)"%(mp.cpu_count()),
        type=int,   default=int(mp.cpu_count()//2) )
    parser.add_argument('-plot', dest="plotfigs", action="store_true",
        help="If set, will generate plots of the fitting results under \
                ./Output/*targname_*band/figs/main_step3_*band_*runnumber")
    parser.add_argument('-n_use', dest="nights_use", action="store",
        help="If you don't want to process all nights under the ./Input/*target/ \
                folder, specify an array of night you wish to process here. \
                e.g., [20181111,20181112]",
        type=str, default='')

    parser.add_argument('-skipmod',    dest="skipmod",          action="store_true",
        help="If set, will skip computing residuals if you have done so already")
    parser.add_argument('-DeBug', dest="debug", action="store_true",
        help="If set, DeBug logging will be output, as well as (lots of) extra plots.")
    parser.add_argument('-sk_check', dest="skip", action="store_true",
        help="If set, will skip the input parameters check. Handy when running \
                multiple targets line by line")

    parser.add_argument(
        '--version', action='version',
        version='%(prog)s {}'.format(_igrins_version))

    return parser.parse_args()
    
    parser.add_argument(
        '--version', action='version', 
        version='%(prog)s {}'.format(_igrins_version))
    
    return parser.parse_args()
