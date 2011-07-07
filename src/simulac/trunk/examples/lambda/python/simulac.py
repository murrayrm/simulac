# simulac.py - python support for simulac
# RMM, 11 May 2010

# Modules that we need
import os                       # operating system calls
import re                       # regular expressions library
import numpy as np              # numerical library

# Functions to import
from matplotlib.pyplot import plot, errorbar, axis, hold, legend
from ciplot import ciplot       # confidence interval plot

#
# Functions for reading simulac data sets
#

def readSet(path=None):
    # Determine what path to use for data
    if (path == None):
        # Get the path through a dialog box
        import tkFileDialog;
        path = tkFileDialog.askdirectory(mustexist=True);
    
    # Make sure that we got a path
    if (path == None):
        print("No path to data specified");
        return None;

    # Create a dictionary for storing the data
    simset = dict();

    # Import the simulation configuration and store local data
    simset['config'] = readSetConfiguration(path, "lambda_setup.py");
    lastPath = path;            # If we get here, path is valid

    # Load up the data files
    simset['setdata'] = readSetData(path);

    return simset

# Read the simulation parameters from a directory
def readSetConfiguration(path=None, setupFile="simulac_setup.py"):
    # Create a dictionary for this simulation set
    setcfg = dict();
    setcfg['path'] = path;
    setcfg['setupFile'] = setupFile;

    # Open up the simulation setup file
    fp = open(path + "/" + setupFile);

    # Go through the file and parse out the individual lines
    while True:
        # Read a line, check for end of file, strip off the newline
        line = fp.readline();
        if (line == ""): break;
        line.rstrip();

        # Parse the different types of lines, getting rid of comments first
        if (re.match("\s*#.*", line) or re.match("\s*$", line)):
            continue;
        elif (re.match("%.*", line)):
            print("Warning: found line starting with %; MATLAB file?");
            continue;

        # Look for an assignment statement
        m = re.match("([^\s]*)\s*=\s*(.*);", line);
        if (m): 
            setcfg[m.group(1)] = eval(m.group(2));
        else:
            print("Warning: couldn't parse line '%s'" % line);

    return setcfg;

# Load the data from a directory
def readSetData(path=".", datafile_regexp=".*\.dat"):
    # Load each of the data files in the directory
    nfiles = 0; runlength = 0; varcount = 0; setdata = None;
    for file in os.listdir(path):
        # Skip files that are not of the right form
        if (re.search(datafile_regexp, file) == None):
            continue;

        # Load the results of the simulation and store as 3-dim array
        simdata = np.loadtxt(path + "/" + file, comments='%')
        rundata = np.reshape(simdata, (1, simdata.shape[0], simdata.shape[1]));
        ++nfiles;

        # Store the data in an ndarray for later access
        if (runlength == 0):
            # For the initial run, just save the data
            setdata = rundata;

            # Store the length of the run and number of variables
            runlength = simdata.shape[0];
            varcount = simdata.shape[1];

        else:
            # Make sure the data is the same length
            if (simdata.shape[0] != runlength):
                print("%s: run length mismatch; skipping file" % file);
            elif (simdata.shape[1] != varcount):
                print("%s: variable count mismatch; skipping file" % file);
            else:
                # Store the data as ndarrays (n = 3)
                setdata = np.concatenate((setdata, rundata));
                
    return setdata

#
# Functions for reading simulac configuration files
#

# Read a simulac configuration (outline) file
def readSimulacConfiguration(file=None, path="."):
    # Open the main outline file for reading
    outline = open(path + "/" + file);

    # Now go through and parse the contents of the file
    while True:
        line = outline.readline();      # grab a line of input
        if (line == ""): break;         # end of file

        m = re.match("(Cel\.[^\s]*)\s*", line);
        if (m):
            print("Reading cell '%s'" % m.group(1));
            cell = Cell(m.group(1));

    return cell

# Read a simulac parameter file
def readSimulacParameters(file=None):
    return None

#
# Plotting functions
#

# Function to generate a plot of protein concentrations over multiple runs
def distPlot(simset,
                species=("CICI", "CroCro", "CII", "CIII", "N", "Qa"),
                plotdev=True, ploterr=0, colors=None, alpha=0.25,
                Tmax = None, Cmax = None, plotLegend=True, runlist=None
                ):

    # Utility function to parse list-like arguments
    def parseListArg(arg, i):
        if (getattr(arg, '__getitem__', False)): return arg[i % len(arg)]
        return arg;

    # Extract out the data and configuration information
    config = simset['config']
    setdata = simset['setdata']
    if (runlist != None): setdata = setdata[runlist,:,:];
    if (setdata.shape[0] == 0): return (0, 0)

    #
    # Go through the set data and compute means and variances
    #

    # Store the time variable
    time = setdata[0,:,config['time_index']]/60;

    # Keep track of the range for later scaling of the plot
    if (Tmax == None): Tmax = max(time);
    if (Cmax == None): Cmax = 0;

    # Go through each variable and compute the statistics
    for i in range(len(species)):
        var = species[i];

        # Figure out the color to use
        if (colors == None or len(colors) < i):
            # Choose a color from a prespecified list
            color = ('b', 'r', 'g', 'c', 'm', 'k')[i % 6];
        else:
            color = colors[i];

        # Figure out what to plot
        s_plotdev = parseListArg(plotdev, i);
        s_ploterr = parseListArg(ploterr, i);

        s_Cmax = plotConcentration(setdata, config, var, color, 
                                   plotdev=s_plotdev, ploterr=s_ploterr, 
                                   alpha=alpha)
        hold(True);
        if (s_Cmax == None): continue

        Cmax = max((Cmax, s_Cmax))

    # Reset the axes of the plot to show everything
    axis([0, Tmax, 0, Cmax]);
    hold(False);

    # Plot a legend the first time around
    if (plotLegend): legend(species);

    # Return the scale used for this plot
    return (Tmax, Cmax)

# Plot the activity level of a promoter
def activityPlot(simset, list, color='b', 
                 plotdev=True, ploterr=0, alpha=0.25, window=60,
                 runlist=None, scale=1):
    # Extract out the data and configuration information
    config = simset['config']
    setdata = simset['setdata']
    if (runlist != None): setdata = setdata[runlist,:,:];
    if (setdata.shape[0] == 0): return (0, 0)

    # Compute the statistics for this value
    stats = activityStats(setdata, config, list, window=window, scale=scale)
    time = setdata[0,:,config['time_index']]/60;

    # Make sure we found the data
    if (stats == None): return None;

    # Store the statistics for the species
    (s_mean, s_std) = stats;
    s_err = s_std / np.sqrt(setdata.shape[0]);

    return plotDistribution(time, (s_mean, s_std, s_err), color=color, 
                            plotdev=plotdev, ploterr=ploterr, alpha=alpha)


# Plot a species concentration
def plotConcentration(setdata, config, var, color='b', 
                      plotdev=True, ploterr=0, alpha=0.25):
    # Compute the statistics for this value
    stats = computeStats(setdata, config, var, conc=True)
    time = setdata[0,:,config['time_index']]/60;

    # Make sure we found the data
    if (stats == None): return None;

    # Store the statistics for the species
    (s_mean, s_std) = stats;
    s_err = s_std / np.sqrt(setdata.shape[0]);

    return plotDistribution(time, (s_mean, s_std, s_err), color=color, 
                            plotdev=plotdev, ploterr=ploterr, alpha=alpha)

# Plot a distribution
def plotDistribution(time, data, color='b', 
                     plotdev=True, ploterr=0, alpha=0.25):
    s_mean = data[0];
    s_std = data[1];
    s_err = data[2];

    # Plot the desired statistics
    if (plotdev == True):
        # Confidence interval plot for the standard deviation
        ciplot(time, s_mean-s_std, s_mean+s_std, color=color, alpha=alpha)
        hold(True);

        # Compute the range for the axis scaling
        Cmax = max(s_mean) + max(s_std)
    else:
        # Figure out the the max concentration later
        Cmax = 0

    if (ploterr != 0):
        # Save the stride
        s = ploterr;

        # Plot mean with error bars
        errorbar(time[::s], s_mean[::s], yerr = s_err[::s], 
                 ecolor=color, fmt=None)
        hold(True);

        if (Cmax == 0): Cmax = max(s_mean) + max(s_err);

    # Finall plot the mean on top of everything
    plot(time, s_mean, color)
    hold(True);

    # Make sure we pass back the max concentation
    if (Cmax == 0): Cmax = max(s_mean);

    return Cmax

#
# Utility functions
#

# Compute the statistics across a data set
def computeStats(setdata, config, var, conc=False):
    #
    # Figure out what species index to use
    #

    # If we were given a string, try looking it up
    if (type(var) == str):
        # First try as a species
        species = config.get("species_" + var + "_index");

    # Make sure that we got something sensible
    if (species == None or type(species) != int):
        # Couldn't find the species name
        print("computeStats: could not find variable '%s'" % var)
        return None;

    # Get the volume variable
    if (conc):
        volume = config['volume_index'];
    else:
        volume == None;

    # Compute the desired statistics
    if (volume == None):
        data_mean = np.mean(setdata[:,:,species], 0)
        data_std = np.std(setdata[:,:,species], 0)
    else:
        # Normalize by volume
        data_mean = np.mean(np.divide(setdata[:,:,species], 
                                      setdata[:,:,volume]), 0)
        data_std = np.std(np.divide(setdata[:,:,species], 
                                    setdata[:,:,volume]), 0)
    return (data_mean, data_std)

# Calculate activity level
def activityStats(setdata, config, list, window=60, scale=1):
    # Figure out what the time step per iteration is
    timeIndex = config['time_index']
    timeStep = setdata[0, 1, timeIndex] - setdata[0, 0, timeIndex]

    # Decide on the integer window size to use when going through the data
    winsize = np.floor(window/timeStep);

    # Figure out the scale: consists of three factors
    #   (1) the scale factor passed as an argument
    #   (2) window quantization
    #   (3) units (always report in polymerases per minute)
    scale *= (window/timeStep) / winsize * (60.0/window);

    # Sum up the RNAP counts across the list of promoters
    rnapCounts = np.sum(setdata[:,:,list], axis=2);

    # Now figure out how many polymerases per window period
    rnapRate = np.zeros(rnapCounts.shape);

    for i in range(int(winsize), rnapRate.shape[1]):
        rnapRate[:,i] = (rnapCounts[:,i] - rnapCounts[:,i-winsize]) * scale;

    # Finally, compute the statistics for the rates across the runs
    data_mean = np.mean(rnapRate, axis=0)
    data_std = np.std(rnapRate, axis=0)

    return (data_mean, data_std)

# Integrate the concentration of a species
def integrateConcentration(data, species_index, volume_index, time_index=0):
    value = 0;
    for timept in range(0, data.shape[0]-1):
        value += data[timept, species_index]/data[timept, volume_index] * \
            (data[timept+1, time_index] - data[timept, time_index])
    return value
