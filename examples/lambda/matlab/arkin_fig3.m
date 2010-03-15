% MATLAB function to generate Fig 3 in Arkin et al
% RMM, 21 Jan 10
%
% Usage: arkin_fig3(dir)
%

function arkin_fig3(dir)

% Load the parameters that define the simulation
filename = [dir, '/lambda_setup'];
if (exist([filename, '.mat'], 'file'))
  load(filename);			% load MAT file from runbatch
end
if (exist([filename, '.m'], 'file'))
  savedir = pwd;			% save our current location
  cd (dir);				% go into the data directory
  lambda_setup;				% load indices from simulac
  cd (savedir)				% return to original directory
end

% Define a few shorter symbols for indices
sl_time = sl_time_index + 1 ;
sl_vol = sl_volume_index + 1;
sl_CI2 = sl_species_CICI_index + 1;
sl_Cro2 = sl_species_CroCro_index + 1;

% Initilize the plot
clf; subplot(221); 
title('All cells, individual runs');

% Run through all of the data
row = 1;				% keep track of good data
clear runfate;
for run = 1:length(parlist)
  for trial = 1:Ntrials
    % Create the filename
    filename = sprintf('%s/%s-%c%d.dat', dir, basename, 'a'+run-1, trial);
    fprintf(1, 'Loading %s\n', filename);
    
    % Load the results of the simulation
    simulac = load(filename);
    
    % Determine the fate
    runfate(row) = fate(simulac, 10, 1);
    
    % Save the results we want
    if (row == 1)
      % Initialize the data matrices; rows = data, cols = run
      ci2 = (simulac(:, sl_CI2) ./ simulac(:, sl_vol))';
      cro2 = (simulac(:, sl_Cro2) ./ simulac(:, sl_vol))';
      time = simulac(:, sl_time)'/60;
      row = row + 1;
      fprintf(1, '  iniitializing length to %d\n', length(time));

    else
      % Append data for this run
      if (size(simulac(:, sl_time), 1) ~= size(ci2, 2))
	% Incomplete run; just ignore it
	fprintf(2, '  inconsistent length: %d was %d => ignored\n', ...
	  size(simulac(:, sl_time), 1), size(ci2, 2));
      else
	ci2 = vertcat(ci2, (simulac(:, sl_CI2) ./ simulac(:, sl_vol))');
	cro2 = vertcat(cro2, (simulac(:, sl_Cro2) ./ simulac(:, sl_vol))');
	row = row + 1;
      end
    end
    
    % Plot the results
    plot(simulac(:,sl_time)/60, simulac(:,sl_CI2)/simulac(:,sl_vol), 'b', ...
      simulac(:,sl_time)/60, simulac(:,sl_Cro2)/simulac(:,sl_vol), 'g', ...
      simulac(:,sl_time)/60, simulac(:,sl_vol)*10, 'r');
    hold on;
  end
end

%
% Plot theta data in a different way
%
subplot(222); hold on;
title('All cells');

% Compute the mean and variance for the variables
ci2_mean = mean(ci2); ci2_std = std(ci2);
cro2_mean = mean(cro2); cro2_std = std(cro2);

% Plot the means and legend
ph = plot(time, ci2_mean, 'b', time, cro2_mean, 'g', 'LineWidth', 2);
lh = legend('CI2', 'Cro2');
legend(lh, 'Location', 'SouthEast');

% Plot the variance of the CI2 concentration
ciplot(ci2_mean-ci2_std, ci2_mean+ci2_std, time, 'b');
alpha(0.25);

% Plot the variance of the Cro2 concentration
ciplot(cro2_mean-cro2_std, cro2_mean+cro2_std, time, 'g');
alpha(0.25)

% Replot the means of the variables so they are on top
plot(time, ci2_mean, 'b', time, cro2_mean, 'g', 'LineWidth', 2);

%
% Plot data for lysogens
%

% Extract at the data for lysogens
lysogens = find(runfate == 1);
fprintf(1, 'Found %d lysogenic runs\n', length(lysogens));

subplot(224); hold on;
title('Lysogenic subpopulation');

% Compute the mean and variance for the variables
ci2_mean = mean(ci2(lysogens,:), 1); ci2_std = std(ci2(lysogens,:), 1);
cro2_mean = mean(cro2(lysogens,:), 1); cro2_std = std(cro2(lysogens,:), 1);

% Plot the variance of the CI2 concentration
ciplot(ci2_mean-ci2_std, ci2_mean+ci2_std, time, 'b');
alpha(0.25);

% Plot the variance of the Cro2 concentration
ciplot(cro2_mean-cro2_std, cro2_mean+cro2_std, time, 'g');
alpha(0.25)

% Now plot the means of the variables
plot(time, ci2_mean, 'b', time, cro2_mean, 'g', 'LineWidth', 2);

%
% Plot data for lytics
%

% Extract at the data for lytics
lytics = find(runfate == 0);
fprintf(1, 'Found %d lytic runs\n', length(lytics));

subplot(223); hold on;
title('Lytic subpopulation');

% Compute the mean and variance for the variables
ci2_mean = mean(ci2(lytics,:), 1); ci2_std = std(ci2(lytics,:), 1);
cro2_mean = mean(cro2(lytics,:), 1); cro2_std = std(cro2(lytics,:), 1);

% Plot the variance of the CI2 concentration
ciplot(ci2_mean-ci2_std, ci2_mean+ci2_std, time, 'b');
alpha(0.25);

% Plot the variance of the Cro2 concentration
ciplot(cro2_mean-cro2_std, cro2_mean+cro2_std, time, 'g');
alpha(0.25)

% Now plot the means of the variables
plot(time, ci2_mean, 'b', time, cro2_mean, 'g', 'LineWidth', 2);
