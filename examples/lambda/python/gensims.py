# Python script to generate Figure 3 in St Pierre, Endy 2008
# RMM, 7 Oct 09

# Include some function we will need
import sys
import os
import getopt
import time
from os import system, stat, putenv
from numpy import linspace

# Set the location of the main configuration files
putenv('SIMULAC_CONFDIR', '../config');

# Default parameters to define the simulation runs
Nvolumes = 8
Ntrials = 8
Trial0 = 0
Volumes = linspace(0.4, 1.8, Nvolumes)
Labels = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']
Subdir = "."

#
# Read command line arguments for additional information
#

# Define a function to print out usage information
def usage():
    print "Usage: python gensims.py"
    return;

# Command line processing using getopt
try:                                
    opts, args = getopt.getopt(sys.argv[1:], "hs:", ["help", "subdir="]) 
except getopt.GetoptError, err:
    print(err); usage();
    sys.exit(2)   

# Now go through the options and process them
for opt, arg in opts:
    if opt in ("-h", "--help"):      
        usage()                     
        sys.exit()                  
    elif opt in ("-s", "--subdir"): 
        Subdir = arg

# Generate a log file indicating what we are doing
logfp = open(Subdir + "/gensims.log", "a");
logfp.write("-------\n");
logfp.write("%s: Nvolumes=%d, Ntrials=%d, Trial0=%d\n" % 
            (time.asctime(), Nvolumes, Ntrials, Trial0));

# Run a bunch of simulations for different cell volumes
for i in range(Volumes.size):
    for trial in range(Ntrials):

        # Create the file name for the simulation
        filename = Subdir + "/lambda-" + Labels[i] + str(trial+Trial0) + ".dat"
        logfile = Subdir + "/lambda-" + Labels[i] + str(trial+Trial0) + ".log"

        # See if the file already exists
        if (os.path.exists(filename)):
            print filename, "already exists; skipping"

        else:
            # Log the run
            logfp.write("%s:   volume=%g, trial=%d -> %s\n" % 
                        (time.asctime(), Volumes[i], trial, filename));
            logfp.flush();

            # Run the simulation
            system("Simulac -v " + str(Volumes[i]) +
                   " --config-file=" + Subdir + "/lambda.cfg" + 
                   " --python-setup=" + Subdir + "/" + "lambda_setup.py" +
                   " -o " + filename + " -d 3 -l " + logfile);

# Close the log file
logfp.close();
