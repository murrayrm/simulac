% MATLAB script to run a batch of simulations varying cell volume
% RMM, 11 March 2010
% 
% This script illustrates how to run a batch of simulations and vary
% the cell volume across runs (any other parameter could be used by
% changing the command line arguments).
%
% This script should be run from the directory containing the system
% description files ('Outline.*, Kin.*, Seq.*, etc) for the system.
% You can edit the 'simcmd' string to point to where the Simulac command
% lives on your system.

% Define the parameters for the run
Tmax = 2000;				% running time (in seconds)
Tprint = 10;				% print interval (in seconds)
Ntrials = 10;				% number of trials per parameter value
parlist = linspace(0.4, 2, 10);		% values of the parameters to use

% Some parameters that define where things live
simulac = 'Simulac';			% name/location of simulac executable
datadir = './data-vol';			% directory for storing output
confdir = './config';			% directory with system description
basename = 'lambda';			% base name for simulation files
overwrite = 1;				% overright files if they exist

% Create a command that will be executed for each parameter/trial
% Update the path to Simulac as needed on your system
% '%g' will be replaced by the parameter value
simcmd = [...
  simulac...				% Command name
  ' -v %g'...				% Argument to vary
  ' --matlab-setup=' datadir '/lambda_setup.m'...	% generate setup file
  ' Outline.Lambda'...			% system description
  ' %d %d'...				% time params (set in loop)
  ' > ' datadir '/%s'...		% location to sae output
];

% Set the SIMULAC_CONFDIR environment variable
setenv('SIMULAC_CONFDIR', confdir);

% Make sure the data directory exists
if (exist(datadir, 'file') == 0)
  mkdir('.', datadir);
end
if (exist(datadir, 'dir') ~= 7)
  fprintf(2, '%s is not accessible', datadir);
  return;
end

% Save all of the parmeters that we used so that we know what was run
setupfile = [datadir '/' basename '_setup'];
save(setupfile, 'Tmax', 'Tprint', 'Ntrials', 'parlist', ...
  'simulac', 'datadir', 'basename');

% Iterate over the parameters
for run = 1:length(parlist)
  for trial = 1:Ntrials
    % Figure out the filename to use for output
    outfile = sprintf('%s-%c%d.dat', basename, 'a'+run-1, trial);
    
    % See if the data file already exists
    if (exist([datadir '/' outfile], 'file') && ~overwrite) 
      fprintf(2, 'Skipping %s\n', outfile);
      continue;
    end

    % Generate the command to run
    cmd = sprintf(simcmd, parlist(run), Tmax, Tprint, outfile);
    
    % Run the simulation
    system(cmd);
  end
end
