% Simple MATLAB script to run a simulation
% RMM, 29 Sep 09
% 
% This script illustrates how to run a simulation from within MATLAB
% and load the results.  The file 'Outline.Lambda' needs to exist in
% the directory from which this script is executed.  Nothing fancy.

% Run a simulation
system('Simulac Outline.Lambda 2000 10 > simulac.dat');

% Load the results of the simulation 
data = load('simulac.dat');


