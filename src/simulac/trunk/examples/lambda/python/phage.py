# phage.py - functions for processing lambda phage simulac data
# RMM, 12 May 2010

# Modules that we need
import re                       # regular expressions library
import numpy as np              # numerical library
import simulac                  # functions to support simulac

# Functions to import
from matplotlib.pyplot import * # commands for generating plots
from ciplot import ciplot       # confidence interval plot
from simulac import computeStats  # compute statistics for a data set

# Define a (empty) class for keeping track of the phage configuration
class PhageConfiguration(object): pass

__doc__ = """
  setup()               set up phage configuration

  Plotting utilities
  ------------------
  instSubplot()         plot the results of a simulation instance
  distSubplot()         plot the distribution of a set of runs
  compareFate()         compare different methods for determining phage
"""

def setup(sc):
    # Start copying elements of the simulac configuration
    pc = PhageConfiguration();
    pc.simulacConfiguration = sc;
    pc.config_file = sc['config_file'];
    pc.time_index = sc['time_index'];
    pc.volume_index = sc['volume_index'];

    # Now augment this with some local definitions for indices
    pc.CI2_index = sc['species_CICI_index'];
    pc.Cro2_index = sc['species_CroCro_index'];
    pc.CII_index = sc['species_CII_index'];
    pc.CIII_index = sc['species_CIII_index'];
    pc.N_index = sc['species_N_index'];

    # Load information that may not be present (depending on system config)
    pc.Q_index = sc.get('species_Q_index');
    pc.Qa_index = sc.get('species_Qa_index');
    pc.VI = sc.get('cell_size_initial');
    pc.V0 = sc.get('cell_size_reference');
    pc.GrowthRate = sc.get('cell_growth_rate');

    # Figure out the configuration information
    m = re.match(".*/(.*)_(.*)/(.*)-(.*)_(.*)$", sc.get('path'));
    if (m == None):
        print("phage: can't process path '%s'" % pc.config_file);
        pc.set_name = 'Unknown';
        pc.set_date = 'Unknown';
        pc.run_name = 'Unknown';
        pc.run_date = 'Unknown';
        pc.run_time = 'Unknown';
    else:
        pc.set_name = m.group(1);
        pc.set_date = m.group(2);
        pc.run_name = m.group(3);
        pc.run_date = m.group(4);
        pc.run_time = m.group(5);

    # Figure out the MOI based on the configuration file
    if (re.search("LambdaQs.*", pc.set_name) != None):
        # Split promoter for PR/PRE + PAQ
        pc.MOI = sc['operator_OperatorPRPRM_0_index'] - \
            sc['operator_OperatorPL_0_index'];
    elif (re.search("LambdaQ.*", pc.set_name) != None):
        # Contains PAQ promoter
        print("phage: MOI computation not implemented for '" +
              pc.set_name + "'; assuming MOI = 8")
        pc.MOI = 8;
    elif (re.search("Lambda.*", pc.set_name) != None):
        # Original version
        pc.MOI = sc['operator_OperatorPRE_0_index'] - \
            sc['operator_OperatorPL_0_index'];
    else:
        print("phage: can't determine MOI for '" +
              pc.set_name + "'; assuming MOI = 8")
        pc.MOI = 8;

    # Return the phage configuration description
    return pc;

# Define a function to get list of PRE promotor states
def PREcols(setup, MOI=None):
    if (MOI == None):
        MOI = setup.MOI;

    # Figure out the increment between promoter states
    sc = setup.simulacConfiguration;
    if (sc.get('operator_OperatorPRE_1_index') == None):
        incr = 1;
    else:
        incr = sc['operator_OperatorPRE_1_index'] - \
            sc['operator_OperatorPRE_0_index'];

    # Information on promotor state for PRE
    return tuple(range(sc['operator_OperatorPRE_0_index'], 
                       sc['operator_OperatorPRE_0_index']+incr*MOI, incr));

# Define a function to get list of PR promotor states
def PRPRMcols(setup, MOI=None):
    if (MOI == None):
        MOI = setup.MOI;

    # Figure out the increment between promoter states
    sc = setup.simulacConfiguration;
    if (sc.get('operator_OperatorPRPRM_1_index') == None):
        incr = 1;
    else:
        incr = sc['operator_OperatorPRPRM_1_index'] - \
            sc['operator_OperatorPRPRM_index']

    # Information on promotor state for PRE
    return tuple(range(sc['operator_OperatorPRPRM_0_index'],
                       sc['operator_OperatorPRPRM_0_index']+incr*MOI, incr));

# Define a function to get list of PRE RNAP counts
def PREcounts(setup, MOI=None):
    if (MOI == None):
        MOI = setup.MOI;

    # Figure out the increment between promoter states
    sc = setup.simulacConfiguration;
    if (sc.get('promoter_PRE_1_index') == None):
        incr = 1;
    else:
        incr = sc['promoter_PRE_1_index'] - sc['promoter_PRE_index']

    # Information on promotor state for PRE
    return tuple(range(sc['promoter_PRE_index'], 
                       sc['promoter_PRE_index'] + incr*MOI, incr));

# Define a function to get list of PRE RNAP counts
def PRcounts(setup, MOI=None):
    if (MOI == None):
        MOI = setup.MOI;

    # Figure out the increment between promoter states
    sc = setup.simulacConfiguration;
    if (sc.get('promoter_PR_1_index') == None):
        incr = 1
    else:
        incr = sc['promoter_PR_1_index'] - \
            sc['promoter_PR_index']

    # Information on promotor state for PRE
    return tuple(range(sc['promoter_PR_index'],
                       sc['promoter_PR_index'] + incr*MOI, incr));

# Define a function to get list of PL RNAP counts
def PLcounts(setup, MOI=None):
    if (MOI == None):
        MOI = setup.MOI;

    # Figure out the increment between promoter states
    sc = setup.simulacConfiguration;
    if (sc.get('promoter_PL_1_index') == None):
        incr = 1
    else:
        incr = sc['promoter_PL_1_index'] - \
            sc['promoter_PL_index']

    # Information on promotor state for PL
    return tuple(range(sc['promoter_PL_index'],
                       sc['promoter_PL_index'] + incr*MOI, incr));

# Define a function to get list of PRM RNAP counts
def PRMcounts(setup, MOI=None):
    if (MOI == None):
        MOI = setup.MOI;

    # Figure out the increment between promoter states
    sc = setup.simulacConfiguration;
    if (sc.get('promoter_PRM_1_index') == None):
        incr = 1
    else:
        incr = sc['promoter_PRM_1_index'] - sc['promoter_PRM_index']

    # Information on promotor state for PRE
    return tuple(range(sc['promoter_PRM_index'],
                       sc['promoter_PRM_index'] + incr*MOI, incr));

def setFate(simset, method=0):
    setdata = simset['setdata'];
    pc = setup(simset['config']);

    runfate = None;
    for run in range(0, setdata.shape[0]):
        simdata = setdata[run,:,:];

        # Calculate the fate for this run
        if (runfate == None):
            runfate = [fate(simdata, pc, method=method)]
        else:
            runfate = np.concatenate((runfate, [fate(simdata, pc, 
                                                     method=method)]));

    return runfate

def fate(data, config=None, method=0, PREmin=1):
    """Determine the cell fate from data run.

    * Method 0 (Arkin et al, 1998): A cell becomes committed to
      lysogeny if there is (i) a sufficient time-integrated
      concentration of CII to activate PRE, and (ii) [CI2] > [Cro2] at
      the end of the 35-min cell cycle. Activation of PRE was defined
      as an average activation level of one open-complex per 2 min
      over a contiguous 4-min period. This level of CII production
      would also activate the other CII-dependent lambda promoters
      that function in execution of the lysogenic pathway (McAdams and
      Shapiro 1995). CI2 concentration greater than that of Cro2 at
      the end of the cell cycle is an additional indication that
      activation of PRE occurred early enough and was productive
      enough to lock on the CI feedback loop.

    * Method 1: Arkin condition (i) - a sufficient time-integrated
      concentration of CII to activate PRE

    * Method 2: Arkin condition (ii) - [CI2] > [Cro2] at the end of
      the 35-min cell cycle"""

    import numpy as np

    # Initialization
    activeTwoMin = False                # PRE active in last two minutes?
    lastActive = data[0, config.time_index]; # Last time we saw an active PRE

    # Figure out what the time step per iteration is
    timeStep = data[1, config.time_index] - data[0, config.time_index];

    # For method 2, just need to look at final concentration of species
    if (method == 2):
        # Just look at the end concentration of CI2 vs Cro2
        if (data[-1, config.CI2_index] > 
              data[-1, config.Cro2_index]):
            return 1;               # Lysogeny
        else:
            return 0;               # Lysis

    # First condition: check if PRE is sufficiently activated
    # Look through the simulation record step by step
    for i in range(np.size(data, axis=0)):
        curTime = data[i, config.time_index];

        # Make sure we have been running long enough
        if (i * timeStep < 120):
            # Not enough time elapsed to average
            continue
        
        #
        # Compute the average number of RNAPs @ PRE over the last 2 min
        #

        # Figure out the beginning time for a two minute window
        startIndex = int(round(i - 120/timeStep));

        if (method == 0 or method == 1):
            # Now count up the number of RNAPs over all copies of PRE
            PRElist = np.sum(data[startIndex:i, PREcounts(config)], axis=1)

            # Now figure out how many RNAPs have launched since start of window
            PREactive = PRElist[-1] - PRElist[0]

        else:
            # Use the old (incorrect) activity count
            PREactive = np.mean(np.sum(data[startIndex:i, PREcols(config)] == 3, 
                                       axis=1));


        # See if we should start a timer for the 4 minute clock
        if ((PREactive > PREmin) and (activeTwoMin == False)):
            startFourMin = curTime;
            activeTwoMin = True;

        # Keep track of whether we have had active complex in last 2 min
        if (not (PREactive > PREmin)):
            activeTwoMin = False;

        # See whether there has been an active complex for four minutes
        if (activeTwoMin and (curTime - startFourMin > 240)):
            # Successful activation of PRE; now check [CI2] > [Cro2]
            if (method == 1 or 
                data[-1, config.CI2_index] > 
                data[-1, config.Cro2_index]):
                return 1;               # Lysogeny
            else:
                return 0;               # Lysis

    # If we get here, we failed to lysogenize
    return 0;

# Define a function to determine the cell fate
def fateArkin(data, config=None, **keywords):
    return fate(data, config, **keywords)

# Determine fate based on a max threshold
def fateMaxThreshold(data, species_index, volume_index, threshold):
    if (threshold > 0):
        return max(data[:, species_index]/data[:, volume_index]) > threshold;
    else:
        return max(data[:, species_index]/data[:, volume_index]) < -threshold;

# Determine fate based on a terminal state
def fateFinalThreshold(data, species_index, volume_index, threshold):
    if (threshold > 0):
        return (data[-1, species_index]/data[-1, volume_index]) > threshold;
    else:
        return (data[-1, species_index]/data[-1, volume_index]) < -threshold;

# Determine the fate based on integral of the value
def fateIntegralThreshold(data, species_index, volume_index, 
                          threshold, time_index=0):
    # Figure out the integral of the concentration
    value = simulac.integrateConcentration(data, species_index, 
                                           volume_index, time_index)

    # Compare to the threshold
    if (threshold > 0):
        return value > threshold;
    else:
        return value < -threshold;

# Determine fate based on a max comparison

#
# Subplots that are used in creating the figures
#

# Function to generate plot of a single run
def instSubplot(rundata, pc, legendLoc=None):
    # Plot the data
    if (rundata != None):
        # Load the raw data for standard variables
        time = rundata[:, pc.time_index]/60;
        volume = rundata[:, pc.volume_index];
        CI2 = rundata[:, pc.CI2_index];
        Cro2 = rundata[:, pc.Cro2_index];
        CII = rundata[:, pc.CII_index];
        CIII = rundata[:, pc.CIII_index];
        N = rundata[:, pc.N_index];
        if (pc.Qa_index != None):
            Q = rundata[:, pc.Qa_index];
        elif (pc.Q_index != None):
            Q = rundata[:, pc.Q_index];

        # Plot the raw data so that we can see it
        plot(time, Cro2, 'm', time, CI2,  'g',
             time, CII,  'k', time, CIII, 'c',
             time, N,    'b');
        if (pc.Q_index != None):
            hold(True);
            plot(time, Q, 'r');

    # Plot the legend
    if (legendLoc != None):
        if (pc.Qa_index != None):
            legend(('Cro2', 'CI2', 'CII', 'CIII', 'N', 'Qa'), loc=legendLoc)
        elif (pc.Q_index != None):
            legend(('Cro2', 'CI2', 'CII', 'CIII', 'N', 'Q'), loc=legendLoc)
        else:
            legend(('Cro2', 'CI2', 'CII', 'CIII', 'N'), loc=legendLoc)

# Function to generate a plot of protein concentrations over multiple runs
def distSubplot(setdata, pc, **keywords):
    return simulac.distPlot(\
      {'setdata':setdata, 'config':pc.simulacConfiguration}, \
      species=("CICI", "CroCro", "CII", "CIII", "N", "Qa"),
      colors=('g', 'm', 'k', 'c', 'b', 'r'),
      **keywords);

# Function to plot multiple histograms next to each other
# Slighted modified from http://www.mail-archive.com/matplotlib-users@
# lists.sourceforge.net/msg05189.html
def multihist(axes, xvals, bins=10, normed=0, bottom=None,
              align='edge', orientation='vertical', width=None,
              log=False, type='overlap', gap=None, patch_kwargs=None,
              **kwargs):
        
    # some integrity checks up front
    if type == 'bi' and len(xvals) != 2:
        raise ValueError, \
            'need exactly two data sets for "bi" multihist: %d given' \
            % len(xvals)
    if patch_kwargs is not None and len(patch_kwargs) != len(xvals):
        raise ValueError, 'need same number of patch kwargs and data sets'
        
    # calculate the common bins, more or less stolen from numpy.histogram
    # xvals = [np.asarray(x).ravel() for x in xvals]
    if not np.iterable(bins):
        mn = float(min([x.min() for x in xvals]))
        mx = float(max([x.max() for x in xvals]))
        if mn == mx:
            mn -= 0.5
            mx += 0.5
        bins = np.linspace(mn, mx, bins, endpoint=False)
        
    # make the histograms using the common bins
    xn = []
    for x in xvals:
        n, bins = np.histogram(x, bins, range=None, normed=normed)
        xn.append(n)
            
    # build the patches parameters depending on type argument
    if width is None: width = 0.9*(bins[1]-bins[0])
    delta = 0
    offset = 0
    paint_width = width
    stay_on_top = True
    if type == 'beside':
        if np.iterable(width):
            raise ValueError, \
                'no sequence of widths allowed for "beside" multihist'
        width /= len(xn)
        delta = width
        if align == 'edge':
            offset = 0
        elif align == 'center':
            offset = ((len(xn) / -2.0 + 0.5) * width)
        else:
            raise ValueError, 'invalid alignment: %s' % align
        if gap is None:
            gap = 0
        paint_width = width - gap
    elif type == 'bi':
        stay_on_top = False
    elif type != 'overlap':
        raise ValueError, 'invalid multihist type: %s' % type
        
    # build the patches
    patch_list = []
    on_top = True
    for n in xn:
        obins = [b + offset for b in bins]
        if on_top:
            rn = n
        else:
            rn = [-v for v in n]
        if orientation == 'horizontal':
            patches = axes.barh(obins, rn, height=paint_width, left=bottom,
                                align=align, log=log)
        elif orientation == 'vertical':
            patches = axes.bar(obins, rn, width=paint_width, bottom=bottom,
                               align=align, log=log)
        else:
            raise ValueError, 'invalid orientation: %s' % orientation
        patch_list.append(cbook.silent_list('Patch', patches))
        offset += delta
        on_top = on_top and stay_on_top
                   
    for i in range(len(patch_list)):
        if patch_kwargs == None:
            kwa = kwargs
        else:
            kwa = patch_kwargs[i]
        for p in patch_list[i]:
            p.update(kwa)

# Arkin-like plot of various fate choices
def compareFate(simset=None, fateExpressions=None, fateLabels=None):    
    # Determine the fate expressions to be evaluated
    if (fateExpressions == None):
        fateExpressions = (
            'fate(simdata, pc, method=0)',
            'fate(simdata, pc, method=1)',
            'fate(simdata, pc, method=2)'
            );
        if (fateLabels == None):
            fateLabels = ('Arkin', 'PRE active', 'CI2 > Cro');

    # Figure out what labels to use if not set already
    if (fateLabels == None): fateLabels = fateExpressions;

    # Save the number of fate methods we will use
    nMethods = len(fateExpressions);

    # Make sure that we have some data available
    if (simset == None):
        print("fateplots.compare: missing simulation set");
        return None;

    # Load up the simulation data
    pc = setup(simset['config']);
    setdata = simset['setdata'];

    # Compute fates for each of the runs in the data set
    fates = None; volumes = None; 
    runfates = np.zeros((nMethods, setdata.shape[0]));
    for run in range(0, setdata.shape[0]):
        simdata = setdata[run,:,:];

        # Plot the data for this run
        instSubplot(simdata, pc); hold(True);

        # Store the run fates
        simfates = np.zeros((1, nMethods));
        for fate in range(nMethods): 
            simfates[0, fate] = eval(fateExpressions[fate]);
            runfates[fate, run] = simfates[0, fate]

        # Now store the data for later plotting
        if (fates == None):
            fates = simfates;
            volumes = np.array([simdata[0, pc.volume_index]]);
        else:
            fates = np.concatenate((fates, simfates));
            volumes = np.concatenate((volumes, [simdata[0, pc.volume_index]]));

    # Summarize the data
    print("Read %d files with %d time points, %d variables" % setdata.shape);

    # Compute bounds for distribution subplots
    (Tmax, Cmax) = distSubplot(setdata, pc);

    # Initialize the plot
    clf();

    # Compute the fates over each possible volume
    vol_data = np.unique(volumes);
    vol_norm = vol_data;
    fate_mean = np.zeros((fates.shape[1], vol_data.shape[0]));
    fate_stderr = np.zeros((fates.shape[1], vol_data.shape[0]));

    # Now iterate over the various fate methods
    for fate in range(0, fates.shape[1]):
        # Calculate lytics and lysogens
        lytics = np.flatnonzero(np.array(runfates[fate]) == 0);
        lysogens = np.flatnonzero(np.array(runfates[fate]) == 1);
        print("Found %d lytic runs" % lytics.shape[0]);

        subplot(3,nMethods,fate+1)
        if (lysogens.shape[0] != 0):
            distSubplot(setdata[lysogens,:,:], pc);
        title("%s\n(%g%% lysogenic)" %
              (fateLabels[fate], 100*lysogens.shape[0]/setdata.shape[0]));
        axis([0, Tmax, 0, Cmax]);

        subplot(3,nMethods,fate+nMethods+1)
        if (lytics.shape[0] != 0):
            distSubplot(setdata[lytics,:,:], pc);
        axis([0, Tmax, 0, Cmax]);
        xlabel("Time (min)"); 

        # Compute out the fate as a function of volume
        for vol in range(0, vol_data.shape[0]):
            # Compute the indices corresponding to this volume
            vol_list = np.flatnonzero(np.array(volumes == vol_data[vol]));

            # Compute  the mean and stderr for each volume/fate pair
            fate_mean[fate, vol] = np.mean(fates[vol_list, fate]);
            fate_stderr[fate, vol] = \
                np.std(fates[vol_list, fate]) / np.sqrt(vol_list.shape[0]);

        # Generate the plots for each fate type
        subplot(3,nMethods,fate+2*nMethods+1)
        ciplot(vol_norm,
              fate_mean[fate,:]-fate_stderr[fate,:],
              fate_mean[fate,:]+fate_stderr[fate,:], 
              alpha=0.25);
        hold(True);

    # Add labels
    subplot(3,nMethods,1); ylabel("Lysogens, nanomolar"); 
    subplot(3,nMethods,nMethods+1); ylabel("Lytics, nanomolar"); 
    subplot(3,nMethods,2*nMethods+1); ylabel("Percent lysogeny");

    # Replot the means on top (with legend)
    for fate in range(0, fates.shape[1]):
        subplot(3,nMethods,fate+2*nMethods+1)
        plot(vol_norm, fate_mean[fate,:], linewidth=2)
        axis([np.min(vol_norm), np.max(vol_norm), 0, 1])
        xlabel("Volume factor");
