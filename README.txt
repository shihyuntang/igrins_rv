
# Setup and run igrins_rv is easy if you have Conda.

#-------------------------------------------------------------------------------
#-------------------------------------
# Packages installation (part1)
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
# Packages installation (part2)
#-------------------------------------
# Install Telfit pkg...
# Enter the igrins_rv environment and cd to Telluric-Fitter-master, than do
(igrins_rv) --> python setup.py install

# The (igrins_rv) is showing that you are inside the igrins_rv environment.
#-------------------------------------------------------------------------------




#-------------------------------------------------------------------------------
#-------------------------------------
# Running igrins_rv
#-------------------------------------
# There are three steps in igrins_rv, you can run them by
(igrins_rv) --> python main_step1.py -h

# The -h command here means --help
# This command will show the basic information of main_step1.py, and tell you what parameters you can put in.
# An example of running the target CITau:
(igrins_rv) --> python main_step1.py CITau -c 10 -plot
