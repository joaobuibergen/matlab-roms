% Driver for the postprocessing of ROMS model output. 
%
%
% 

close all
clear variables

addpath(genpath('/Users/joao/software/matlab2'))
addpath /Users/joao/software/roms_wilkin/matlab
addpath /Users/joao/software/matlab-roms

% Add necessary folders to PATH env variable
path1 = getenv('PATH');
path1 = [path1 ':/opt/local/bin'];
setenv('PATH', path1);
!echo $PATH

% Get the cwd

cwd = pwd;

% Run identifier

romsID='TEST1';

% log file (we use it to read model and run info)

romsLog = dir([cwd '/Output/log.*' romsID]);

disp([' Reading info from log file: ' romsLog.name])

% Run postprocessing

disp(' Running postprocessing')

romsPostprocessFjords1D


% Run diagnostics

%disp(' Running diagnostics')

%romsDiagnosticsFjords1D
% 
% % Write html report
% 
% disp(' Writing report diagnostics')
% 
% makeReport

