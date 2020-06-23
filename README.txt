
# Setup and run igrins_rv is easy if you have Conda installed.
# Only tested on Ubuntu machine...

#-------------------------------------------------------------------------------
#-------------------------------------
# Packages installation (part 1)
#-------------------------------------
# Packages installation with conda
# One-click to setup environment with all pkg needed (using the with the environment.yml file) by
--> conda env create

# This will create an environment called "igrins_rv", and you can use
--> conda info --envs
# to check all available environments.

# Use
--> conda activate igrins_rv
# or
--> source activate igrins_rv
# to enter the environment.

# If you success, your command line will look like

(igrins_rv) -->

#-------------------------------------
# If you do not have conda, the basic requirement for running igrins_rv
# is python.3.7, and for the version for individual pkg that tested are listed as follow:
    - python=3.7
    - matplotlib=3.1.3
    - numpy=1.18.1
    - pandas=1.0.3
    - scipy=1.4.1
    - astropy=4.0
    - multiprocess
    - cython=0.29.15
    - requests=2.23.0
    - nlopt=2.6.1
    - pip:
      - pysynphot=0.9.14


#-------------------------------------
# Packages installation (part 2)
#-------------------------------------
# Note the origin Telfit pkg is on https://github.com/kgullikson88/Telluric-Fitter
# by kgullikson88
# One line in the original code in ./src/TelluricFitter.py is modified:
# line 672: resolution=None, vac2air=self.air_wave)
#           --> resolution=None, vac2air=False)

# Install Telfit pkg... and is the tricky part...
# First, check your mechine's gcc version by

--> gcc --version

# Telfit will not work with version later than v9.X.X
# (not sure for v8.X.X)
# v5.5.0 & v7.5.0 work fine

# also, Telfit only recognize "gcc", so please make sure the new installed gcc
# is callable by --> gcc

# If you intend to use more the 6 threads to run this pipeline.
# You are recommended to go in the "setup.py", and modify "num_rundirs = "
# to a number that is THREE times the threads you intend to use.

# Now, enter the igrins_rv environment (because be need to install Telfit in the environment)
# and cd to Telluric-Fitter-master, than do

(igrins_rv) --> python setup.py build
(igrins_rv) --> python setup.py install

# The (igrins_rv) is showing that you are inside the igrins_rv environment.
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#-------------------------------------
# Run igrins_rv
#-------------------------------------
# There are three steps in igrins_rv, you can run them by
(igrins_rv) --> python main_step1.py -h

# The -h command here means --help
# This command will show the basic information of main_step1.py, and tell you what parameters you can put in.
# An example of running the target CITau:
(igrins_rv) --> python main_step1.py CITau -c 10 -plot
