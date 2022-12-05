% Driver for the postprocessing of ROMS model output. 
%
%
% 

close all

% Get the computer name for path costumization
[ret, output] = system('hostname');
output=split(strip(output));
computerName=output{1};

if strcmp(computerName,'fvfh70lhq05p.klientdrift.uib.no')

    addpath(genpath('/Users/joao/software/matlab2'))
    addpath /Users/joao/software/roms_wilkin/matlab
    %addpath /Users/joao/software/matlab-roms

    % Add necessary folders to PATH env variable
    path1 = getenv('PATH');
    path1 = [path1 ':/opt/local/bin'];
    setenv('PATH', path1);
    !echo $PATH

elseif strcmp(computerName,'cyclone.hpc.uib.no')

    addpath(genpath('/home/johor1356/software/matlab_roms'))
    addpath /home/johor1356/software/roms_wilkin/matlab
    %addpath /home/johor1356/software/matlab-roms

    % Add necessary folders to PATH env variable
    path1 = getenv('PATH');
    path1 = [path1 ':/opt/bin'];
    setenv('PATH', path1);
    !echo $PATH

end

% Get the cwd if it wasn't supplied by the script invocation

if ~exist('cwd','var')
    cwd = pwd;
end

% Run identifier

if ~exist('romsID','var')
    error('Please provide a run ID.')
end

disp([' Current working directory: ' cwd])
disp([' ROMS Run ID: ' romsID])

% log file (we use it to read model and run info)

romsLog = dir([cwd 'log.*' romsID]);

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

