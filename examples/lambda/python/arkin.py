# arkin.py - routines for generating figures from Arkin et al
# RMM, 11 May 2011
#
# This module defines a number of functions for computing and plotting
# quantities out of the paper by Arkin, Ross and McAdams (1997).
#
#  * arkin.fig3         distributions of proteins for lytics and lysogens
#  * arkin.instSubplot  plot individual instance of the simulation
#  * arkin.distSubplot  plot distribution of runs

# Modules that we need
import os                       # operating system calls
import re                       # regular expressions library
import simulac                  # functions to support simulac
import phage                    # lambda utility functions
import numpy as np              # numerical library
import matplotlib as mpl        # plotting library


# Functions to import
from matplotlib.pyplot import * # commands for generating plots
from phage import instSubplot   # plot an instance of a simulation
from phage import distSubplot   # plot a distribution of data

# 
# Main figures
#
# The methods described here generate plots corresponding to some of the
# main figures in the Arkin et al paper
#

def fig3(simset=None, method=0, settitle=None):
    # See what type of simulation data we were passed
    if (type(simset) == tuple):
        # Assume we were passed a list in standard order
        setdata = simset[0];
        setconfig = simset[1];
    elif (type(simset) == dict):
        # Assume we were passed a standard simset object
        setconfig = simset['config'];
        setdata = simset['setdata'];
    elif (simset == None):
        # User wants to query for data
        simset = simulac.readSet();
        setconfig = simset['config'];
        setdata = simset['setdata'];
    else:
        # I'm confused
        print("Can't determine input data; aborting")
        return None
        
    # Make sure that we have some data available
    if (setdata == None or setconfig == None):
        print("arkin.fig3: missing simulation set information");
        return None;

    # Process the set configuration information
    pc = phage.setup(setconfig);

    # Set some plotting parameters
    mpl.rcParams['legend.fontsize'] = 'small'

    # Initilize the plot
    clf(); subplot(221); 
    if (settitle == None):
        settitle = "%s: %s, %s" % (pc.set_name, pc.run_date, pc.run_time);
    title(settitle);
    ylabel("Molecule count"); 

    # Load each of the data files in the directory
    nfiles = 0; runfate = None;
    for run in range(0, setdata.shape[0]):
        simdata = setdata[run,:,:];
        rundata = np.reshape(simdata, (1, simdata.shape[0], simdata.shape[1]));
        ++nfiles;

        # Calculate the fate for this run
        if (runfate == None):
            runfate = [phage.fate(simdata, pc, method=method)]
        else:
            runfate = np.concatenate((runfate, [phage.fate(simdata, pc, 
                                                           method=method)]));

        # Plot the data for this run
        instSubplot(simdata, pc); hold(True);

    # Show the results so far
    instSubplot(None, pc, legendLoc='upper right')
    show()

    # Summarize the data
    print("Read %d files with %d time points, %d variables" % setdata.shape);

    # Generate plots for all runs
    subplot(222); hold(True);
    title('All cells (N = %d)' % setdata.shape[0]);
    (Tmax_all, Cmax_all) = distSubplot(setdata, pc);
    ylabel("Nanomolar"); 

    # Generate plots for lytics
    lytics = np.flatnonzero(np.array(runfate) == 0);
    print("Found %d lytic runs" % lytics.shape[0]);

    subplot(223); hold(True);
    if (lytics.shape[0] != 0):
        (Tmax_lytics, Cmax_lytics) = distSubplot(setdata[lytics,:,:], pc);
    else:
        (Tmax_lytics, Cmax_lytics) = (0, 0);
    title('Lytic subpopulation (N = ' + str(lytics.shape[0]) + ')');
    ylabel("Nanomolar"); 
    xlabel("Time (min)"); 

    # Generate plots for lysogens
    lysogens = np.flatnonzero(np.array(runfate) == 1);
    print("Found %d lysogenic runs" % lysogens.shape[0]);

    subplot(224); hold(True);
    if (lysogens.shape[0] != 0):
        (Tmax_lysogens, Cmax_lysogens) = distSubplot(setdata[lysogens,:,:], pc)
    else:
        (Tmax_lysogens, Cmax_lysogens) = (0, 0);
    title('Lysogenic subpopulation (N = ' + str(lysogens.shape[0]) + ')');
    ylabel("Nanomolar"); 
    xlabel("Time (min)"); 

    # Now reset the axes to be the same
    Cmax = max((Cmax_all, Cmax_lytics, Cmax_lysogens))
    Tmax = max((Tmax_all, Tmax_lytics, Tmax_lysogens))
    subplot(222); axis([0, Tmax, 0, Cmax]);
    subplot(223); axis([0, Tmax, 0, Cmax]);
    subplot(224); axis([0, Tmax, 0, Cmax]);

    # Show the results
    show();

# Plots to look at the loading of the cell
def loading(simset=None, setdata=None, setconfig=None, settitle=None):
    # Load up the simualtion data
    if (simset != None):
        setconfig = simset['config'];
        setdata = simset['setdata'];

    # Make sure that we have some data available
    if (setdata == None or setconfig == None):
        print("arkin.fig3: missing simulation set information");
        return None;

    # Process the set configuration information
    pc = phage.setup(setconfig);

    # Initilize the plot
    clf(); subplot(221); 
    if (settitle == None):
        settitle = "%s: %s, %s" % (pc.set_name, pc.run_date, pc.run_time);
    title(settitle);

    # Compute the number of RNA polymerases
    time = setdata[0,:,0];
    RNAP_count_mean, RNAP_count_std = \
        simulac.computeStats(setdata, setconfig, 'RNAP')
    RNAP_mean, RNAP_std = \
        simulac.computeStats(setdata, setconfig, 'RNAP', conc=True);

    # Plot the data
    plot(time, RNAP_count_mean, 'r');
    ylabel("RNAP count"); 

    subplot(222);
    plot(time, RNAP_mean, 'r');
    ylabel("RNAP concentration"); 
