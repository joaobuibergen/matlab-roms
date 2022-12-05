% Script to postprocess ROMS 1D setup.


% J. H. Bettencourt, Bergen, August 2022


%  Revision          Details
% ----------         -------
% 22.11.2022         Modified to work together with diagnostics script.

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Definitions 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Application title (from ROMS log file)

  [status,cmdout] = system(strjoin(["gawk 'NR==6 {print $1}' " fullfile(romsLog.folder,romsLog.name)]));
   
  MyAppTITLE = cmdout(1:end-1); % Remove newline character

  disp(['Application title: ' MyAppTITLE])
  
  inputFile = [MyAppTITLE '.in'];

  disp(['Input file: ' inputFile])

  [status,cmdout] = system(strjoin(["gawk '/Header file/ {print $4}' " fullfile(romsLog.folder,romsLog.name)]));

  cppHeaderFile = strip(cmdout);

  disp(['CPP Header file: ' cppHeaderFile])

  str=split(cppHeaderFile,'.');

  MyAppCPP = upper(str{1});

% File to process (his - history file; avg - averages file)

  fileType = 'his';
  
% Define plot variables

  plotVariables={'temp';...
                 'salt';...
                 'rho';...
                 'NO3';...
                 'NH4';...
                 'chlorophyll';...
                 'phytoplankton';...
                 'zooplankton';...
                 'oxygen';...
                 'LdetritusN';...
                 'SdetritusN';...
                 'LdetritusC';...
                 'SdetritusC';...
                 'TIC';...
                 'alkalinity';...
                 'AKv';...
                 'AKt';...
                 'AKs';...
  };

  surfaceForcing={'Pair';...
                  'Tair';...
                  'Uwind';...
                  'Vwind';...
                  'Qair';...
                  'rain';...
                  'swrad';...
                  'lwrad_down';...
                  'SST'
              };
          
  surfaceFluxes={'shflux';...
                 'latent';...
                 'sensible';...
                 'lwrad';...
                 'ssflux';...
                 'EminusP';...
              };
          
  diaBioVariables={'NO3_hadv';...
                   'NO3_vadv';...
                   'NO3_xadv';...
                   'NO3_yadv';...
                   'NO3_vdiff';...
                   'NO3_rate';...
                   'NH4_hadv';...
                   'NH4_vadv';...
                   'NH4_xadv';...
                   'NH4_yadv';...
                   'NH4_vdiff';...
                   'NH4_rate';...
                   'chlorophyll_hadv';...
                   'chlorophyll_vadv';...
                   'chlorophyll_xadv';...
                   'chlorophyll_yadv';...
                   'chlorophyll_vdiff';...
                   'chlorophyll_rate';...
				   'phytoplankton_hadv';...
                   'phytoplankton_vadv';...
                   'phytoplankton_xadv';...
                   'phytoplankton_yadv';...
                   'phytoplankton_vdiff';...
				   'phytoplankton_rate';...
                   'zooplankton_hadv';...
                   'zooplankton_vadv';...
                   'zooplankton_xadv';...
                   'zooplankton_yadv';...
                   'zooplankton_vdiff';...
				   'zooplankton_rate';...
                   'LdetritusN_hadv';...
                   'LdetritusN_vadv';...
                   'LdetritusN_xadv';...
                   'LdetritusN_yadv';...
                   'LdetritusN_vdiff';...
				   'LdetritusN_rate';...
				   'SdetritusN_hadv';...
                   'SdetritusN_vadv';...
                   'SdetritusN_xadv';...
                   'SdetritusN_yadv';...
                   'SdetritusN_vdiff';...
				   'SdetritusN_rate';...
                   'LdetritusC_hadv';...
                   'LdetritusC_vadv';...
                   'LdetritusC_xadv';...
                   'LdetritusC_yadv';...
                   'LdetritusC_vdiff';...
				   'LdetritusC_rate';...
				   'SdetritusC_hadv';...
                   'SdetritusC_vadv';...
                   'SdetritusC_xadv';...
                   'SdetritusC_yadv';...
                   'SdetritusC_vdiff';...
				   'SdetritusC_rate';...
				   'TIC_hadv';...
                   'TIC_vadv';...
                   'TIC_xadv';...
                   'TIC_yadv';...
                   'TIC_vdiff';...
				   'TIC_rate';...
				   'alkalinity_hadv';...
                   'alkalinity_vadv';...
                   'alkalinity_xadv';...
                   'alkalinity_yadv';...
                   'alkalinity_vdiff';...
				   'alkalinity_rate';...
				   'oxygen_hadv';...
                   'oxygen_vadv';...
                   'oxygen_xadv';...
                   'oxygen_yadv';...
                   'oxygen_vdiff';...
				   'oxygen_rate';...
                   'denitrification';...
                   'CO2_airsea';...
                   'pCO2';...
                   'O2_airsea';...
                   'P_Production';...
                   'NO3_uptake';...
                   'SO2C';...
                   'DO2C';...
                   'PAR';...
                   'Nitrif_flux';...
                   'Zoo_respiration';...
                   'Bac_respiration';...
                   'OM_sinking'...
              };

  nPlotVariables = length(plotVariables);
  nSurfaceForcing = length(surfaceForcing);
  nSurfaceFluxes = length(surfaceFluxes);
  nDiaBioVariables = length(diaBioVariables);

  
% Plot types

  plotInitial = 1;                    % Plot profiles of initial conditions
  plotProfiles = 0;                   % Plot profiles at specified time
  plotMaps = 0;                       % Plot time-depth maps
  
% Plot on undisturbed depths?

  plotZeta0 = 1;
  
% Times to plot

  if plotProfiles
      
     profilePlotTimes = datenum(2020,03,05);
     
  end
  
  if plotMaps
      
     mapPlotTimes = -1; % If -1, plot all simulation
     
  end
  
% Sill depth (m)

  depthSill = 180;
  
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Get grid and other metadata 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Output folder

[status,cmdout] = system(strjoin(["gawk '/Output Restart File:/ {print $4}' " fullfile(romsLog.folder,romsLog.name)]));

str=split(cmdout,'//');

romsOutputFolder = [str{1} '/output']

% Data files
  
  romsDataFile = strip(split(ncreadatt([fullfile(romsOutputFolder,MyAppTITLE) '_his.nc'],'/',...
      [ fileType '_file']),','));

  nRomsData = length(romsDataFile);

% Get the names of initial, forcing and output files. If multiple files,
% get their names as a list.

  romsIniFile = ncreadatt([fullfile(romsOutputFolder,MyAppTITLE) '_his.nc'],'/','ini_file');
  disp(['Roms IC File: ' romsIniFile])
  romsFrcFile = strip(split(ncreadatt([fullfile(romsOutputFolder,MyAppTITLE) '_his.nc'],'/','frc_file_01'),','));
  disp(['Roms Forcing Files: ' strjoin(romsFrcFile)])
  nRomsFrc = length(romsFrcFile);
  romsGridFile = ncreadatt([fullfile(romsOutputFolder,MyAppTITLE) '_his.nc'],'/','grd_file');
  disp(['Roms Grid File: ' romsGridFile])
  romsAvgFile = strip(split(ncreadatt([fullfile(romsOutpuFolder,MyAppTITLE) '_his.nc'],'/','avg_file'),','));
  disp(['Roms Averages Files: ' romsAvgFile])
  nRomsAvg = length(romsAvgFile);
  romsDiaFile = strip(split(ncreadatt([fullfile(romsOutputFolder,MyAppTITLE) '_his.nc'],'/','dia_file'),','));
  disp(['Roms Diagnostic Files: ' romsDiaFile])
  nRomsDia = length(romsDiaFile);

% Get ROMS grid (depths for t=0)

  romsGrid = get_roms_grid(romsGridFile, romsDataFile{1});
  
% Get CPP options, model version, run dir, run date

  %romsCPP = split(ncreadatt(romsDataFile,'/','CPP_options'),',');%roms_cpplist([lower(MyAppCPP) '.h']);
  romsCPP = string(roms_cpplist([fullfile(romsLog.folder,MyAppTITLE) '_his.nc']));
  getStr = split(ncreadatt([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'/','history'),',');
  romsVersion = [getStr{1} ' ' getStr{2} ' rev. '...
      ncreadatt([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'/','svn_rev')];
  
  runDateTime = [getStr{3} ', ' getStr{4}];
  
  romsRunDir = ncreadatt([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'/','header_dir');
  romsSourceDir = ncreadatt([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'/','code_dir');
  
% Get baroclinic time step, barotropic time step and number of time steps between output
  
  romsDt = ncread([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'dt');
  romsDtFast = ncread([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'dtfast');
  
  %if strcmp(fileType,'his')
  
      romsNHisOut = ncread([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'nHIS');
  
  %elseif strcmp(fileType,'avg')
      
      romsNAvgOut = ncread([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'nAVG');
      
  %end
  
  romsDOut = romsDt*romsNHisOut; % Interval in secs between output fields
  




% Get units of plot variables (if units attribute is not present, assumed 
% nondimensional variable)

  plotVariablesUnits={};

  for i=1:nPlotVariables
      
     % Get units attribute
     try
         tmpStr = ncreadatt(romsDataFile{1},plotVariables{i},'units');
     catch
         tmpStr = '-';
     end
     tmpStr = replace(tmpStr,'_',' ');
     tmpStr = replace(tmpStr,'meter','m');
     tmpStr = replace(tmpStr,'kilogram','kg');
     tmpStr = replace(tmpStr,'millimole','mmole');
     tmpStr = replace(tmpStr,'milliequivalents','meqvs');
     tmpStr = replace(tmpStr,'milligrams','mg');
     tmpStr = replace(tmpStr,'nitrogen','N');
     tmpStr = replace(tmpStr,'oxygen','O');
     tmpStr = replace(tmpStr,'carbon','C');
     tmpStr = replace(tmpStr,'chlorophyll','Chl');
     
     
     
     plotVariablesUnits{i}=tmpStr;
      
      
  end
  
  surfaceForcingUnits={};

  for i=1:nSurfaceForcing
      
     % Get units attribute
     try
         tmpStr = ncreadatt(romsFrcFile{1},surfaceForcing{i},'units');
     catch
         tmpStr = '-';
     end
     tmpStr = replace(tmpStr,'second','s');
     tmpStr = replace(tmpStr,'watt','W');
     tmpStr = replace(tmpStr,'meter','m');
     surfaceForcingUnits{i}=tmpStr;
      
      
  end
  
  surfaceFluxesUnits={};

  for i=1:nSurfaceFluxes
      
     % Get units attribute
     try
         tmpStr = ncreadatt(romsDataFile{1},surfaceFluxes{i},'units');
     catch
         tmpStr = '-';
     end
     tmpStr = replace(tmpStr,'second','s');
     tmpStr = replace(tmpStr,'watt','W');
     tmpStr = replace(tmpStr,'meter','m');
     surfaceFluxesUnits{i}=tmpStr;
      
      
  end
  
  diaBioVariablesUnits={};

  for i=1:nDiaBioVariables
      
     % Get units attribute
     try
         tmpStr = ncreadatt(romsDiaFile{1},diaBioVariables{i},'units');
     catch
         tmpStr = '-';
     end
     tmpStr = replace(tmpStr,'_',' ');
     tmpStr = replace(tmpStr,'meter','m');
     tmpStr = replace(tmpStr,'kilogram','kg');
     tmpStr = replace(tmpStr,'second','s');
     tmpStr = replace(tmpStr,'millimole','mmole');
     tmpStr = replace(tmpStr,'milliequivalents','meqvs');
     tmpStr = replace(tmpStr,'milligrams','mg');
     tmpStr = replace(tmpStr,'nitrogen','N');
     tmpStr = replace(tmpStr,'oxygen','O');
     tmpStr = replace(tmpStr,'carbon','C');
     tmpStr = replace(tmpStr,'chlorophyll','Chl');
     diaBioVariablesUnits{i}=tmpStr;
      
      
  end
  
% Get forcing time in matlab format

  frcTimeNative=[];

  for i=1:nRomsFrc

      frcTimeNative = [frcTimeNative ncread(romsFrcFile{i},'frc_time')'];

  end

  frcTimeUnits = ncreadatt(romsFrcFile{1},'frc_time','units');
  
  getStr = split(frcTimeUnits);
  
  frcRefTime = datetime([getStr{3} ' ' getStr{4}]);
  
  if strcmp(getStr{1},'days')
  
      frcTime = frcRefTime + days(frcTimeNative);
  
  elseif strcmp(getStr{1},'seconds')
      
      frcTime = frcRefTime + seconds(frcTimeNative);
      
  end
  
  nFrcTimes = length(frcTimeNative);
  
% Get solution time in matlab format

  romsTimeNative=[];

  for i=1:nRomsData

      romsTimeNative = [romsTimeNative ncread(romsDataFile{i},'ocean_time')'];

  end

  romsTimeUnits = ncreadatt(romsDataFile{1},'ocean_time','units');
  
  getStr = split(romsTimeUnits);
  
  romsRefTime = datetime([getStr{3} ' ' getStr{4}]);
  
  romsTime = romsRefTime + seconds(romsTimeNative);
  
  nTimes = length(romsTimeNative);
  
% Get diagnostics time in matlab format

  diaTimeNative=[];

  for i=1:nRomsDia

      diaTimeNative = [diaTimeNative ncread(romsDiaFile{i},'ocean_time')'];

  end

  diaTimeUnits = ncreadatt(romsDiaFile{1},'ocean_time','units');
  
  getStr = split(diaTimeUnits);
  
  diaRefTime = datetime([getStr{3} ' ' getStr{4}]);
  
  diaTime = diaRefTime + seconds(diaTimeNative);
  
  nDiaTimes = length(diaTimeNative);


% ROMS start date
%strjoin(["gawk '/Header file/ {print $4}' " fullfile(romsLog.folder,romsLog.name)])
  [status,cmdout] = system(strjoin(["gawk '/DSTART =/ {print $3}'  " inputFile]));
  
  tmp=split(cmdout(1:end-1),'d');
  day0=str2num(tmp{1});
  ex0=str2num(tmp{2});
  startDay=day0*10^ex0;
  romsStartTime = romsRefTime + days(startDay);

% Print info



% TO BE DONE

% % Write relevant data to output file
% 
%   fprintf(outFID,'ROMSGridFile = %s\n',romsGridFile);
%   fprintf(outFID,'ROMSIniFile = %s\n',romsIniFile);
% 
%   testOut = romsDataFile';
%   fprintf(outFID,'ROMSDataFile = %s\n',testOut{:});%romsDataFile);
%   
%   testOut = romsFrcFile';
%   fprintf(outFID,'ROMSFrcFile = %s\n',testOut{:});%romsFrcFile);
% 
%   fprintf(outFID,'ROMSVersion = %s\n',romsVersion);
%   
%   tmpStr=[];
%   spx = {' '};
%   
%   for i=1:length(romsCPP)
%       
%       tmpStr=strcat(tmpStr,spx,strip(romsCPP(i)));
%       
%   end
%   
%   fprintf(outFID,'CPPOptions = %s\n',strip(tmpStr));
%   fprintf(outFID,'RUNDir = %s\n',romsRunDir);
%   fprintf(outFID,'SOURCEDir = %s\n',romsSourceDir);
%   fprintf(outFID,'RUNDate = %s\n',runDateTime);
%   fprintf(outFID,'STARTDate = %s\n',datestr(romsTime(1),'yyyy/mm/dd HH:MM:SS'));
%   fprintf(outFID,'ENDDate = %s\n',datestr(romsTime(end),'yyyy/mm/dd HH:MM:SS'));
%   fprintf(outFID,'OUTStep = %s\n',num2str(romsDOut/3600));
%   fprintf(outFID,'DT = %s\n',num2str(romsDt));
%   fprintf(outFID,'DTFAST = %s\n',num2str(romsDtFast));
%   fprintf(outFID, 'NVERT = %s\n',num2str(romsGrid.N));
%   fprintf(outFID, 'VTRANSF = %s\n',num2str(romsGrid.Vtransform));
%   fprintf(outFID, 'VSTRETCH = %s\n',num2str(romsGrid.Vstretching));
%   fprintf(outFID, 'THETAS = %s\n',num2str(romsGrid.theta_s));
%   fprintf(outFID, 'THETAB = %s\n',num2str(romsGrid.theta_b));
%   fprintf(outFID, 'TCLINE = %s\n',num2str(romsGrid.Tcline));
%   fprintf(outFID, 'HC = %s\n',num2str(romsGrid.hc));
  
  
  
  
% Get biogeochemical model parameters

bioParVariables={'AttSW';...
                 'AttChl';...
                 'PARfrac'
                 'Vp0';...
                 'I_thNH4';...
                 'D_p5NH4';...
                 'NitriR';...
                 'K_NO3';...
                 'K_NH4';...
                 'K_PO4';...
                 'K_Phy';...
                 'Chl2C_m';...
                 'ChlMin';...
                 'PhyCN';...
                 'R_P2N';...
                 'PhyIP';...
                 'PhyIS';...
                 'PhyMin';...
                 'PhyMR';...
                 'ZooAE_N';...
                 'ZooBM';...
                 'ZooCN';...
                 'ZooER';...
                 'ZooGR';...
                 'ZooMin';...
                 'ZooMR';...
                 'LDeRRN';...
                 'LDeRRC';...
                 'CoagR';...
                 'SDeRRN';...
                 'SDeRRC';...
                 'RDeRRN';...
                 'wPhy';...
                 'wLDet';...
                 'wSDet';...
                 'pCO2air';...
                 };
  bioParValues={};
  bioParUnits={};
  bioParLongNames={};
  
% Write header in bio paramater file

%   fprintf(bioParFID,'Parameter, Symbol, Value, Units\n');
  
  for k=1:length(bioParVariables)
% 	double AttSW ;
% 		AttSW:long_name = "light attenuation by seawater" ;
% 		AttSW:units = "meter-1" ;
    bioParValues{k} = ncread(romsDataFile{1},bioParVariables{k});
    % Have to do some string replacements to shorten variable units 
    bioParLongNames{k} = replace(ncreadatt(romsDataFile{1},bioParVariables{k},'long_name'),',','');
    bioParLongNames{k} = replace(bioParLongNames{k},'phytoplankton','phyto');
    bioParLongNames{k} = replace(bioParLongNames{k},'zooplankton','zoo');
    bioParLongNames{k} = replace(bioParLongNames{k},'carbon','C');
    bioParLongNames{k} = replace(bioParLongNames{k},'Carbon','C');
    bioParLongNames{k} = replace(bioParLongNames{k},'nitrogen','N');
    bioParLongNames{k} = replace(bioParLongNames{k},'Nitrogen','N');
    try 
      bioParUnits{k} = replace(ncreadatt(romsDataFile{1},bioParVariables{k},'units'),'_',' ');
      bioParUnits{k} = replace(bioParUnits{k},'meter','m');
      bioParUnits{k} = replace(bioParUnits{k},'milligram','mg');
      bioParUnits{k} = replace(bioParUnits{k},'parts per million','ppm');
      bioParUnits{k} = replace(bioParUnits{k},'milli','m');
      bioParUnits{k} = replace(bioParUnits{k},'watt','W');     
    catch
      bioParUnits{k} = '-';
    end
    
%     % Write to output file
%     outStr=[bioParLongNames{k} ',' ...
%         replace(bioParVariables{k},'_','\_') ',' num2str(bioParValues{k}) ',' bioParUnits{k}];
%     fprintf(bioParFID,'%s\n',outStr);
    
  end

  Values=cell2mat(bioParValues)';
  Description=char(bioParLongNames);
  Units=char(bioParUnits);
  Params=char(bioParVariables);
  bioParTable = table(Values,Units,Description,'RowNames',bioParVariables)
  
  
%   % Close output files
% 
%   fclose(outFID)
%   fclose(bioParFID);

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Aux stuff
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sec2day = 3600*24; % Seconds to days

  cellArea = 1./((romsGrid.pm).*(romsGrid.pn));

  nPlot = romsGrid.N;
  zPlot = squeeze(romsGrid.z_r(3,3,:));
  wPlot = squeeze(romsGrid.z_w(3,3,:));
  kBasin = find(zPlot<-1*depthSill,1,'last');
  k100 = find(zPlot>=-1*100.,1,'first');
  k40 = find(zPlot>=-1*40.,1,'first');

  cellThickness = diff(romsGrid.z_w,1,3);
  cellVolume = repmat(cellArea,[1 1 nPlot]).*cellThickness;
  cellVolume = squeeze(cellVolume(3,3,:));


  basinVolume = sum(cellVolume(1:kBasin));
  vol100 = sum(cellVolume(k100:nPlot));

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Plot parameters
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  posMap3 = [680 240 730 740];
  posMap4 = [680 190 730 790];
  posMap5 = [680  30 730 950];

% Figure counter

  kf = 0;

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Plot model evolution
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% T,S contour plots (upper 100m)

  nPlot = romsGrid.N;
  zMap = repmat(squeeze(romsGrid.z_r(1,1,:)),[1 nTimes]);
  zMapDia = repmat(squeeze(romsGrid.z_r(1,1,:)),[1 nDiaTimes]);
  tMap = repmat(romsTime,[nPlot 1]);
  tMapDia = repmat(diaTime,[nPlot 1]);
  
  kf=kf+1;
  figure(kf) % Temp, salt, rho
  subplot 311
  plotVar = 'temp';
  index = find(strcmp(plotVariables, plotVar));
  plotVarUnits = plotVariablesUnits{index};
  plotMap=[];
  for i=1:nRomsData
      plotMap = [plotMap squeeze(ncread(romsDataFile{i},plotVar,[3 3 k100 1],[1 1 Inf Inf]))];
  end
  v=5:14;
  [C,H] = contour(datenum(tMap(k100:end,:)),zMap(k100:end,:),plotMap,v,...
      '-k','LineWidth',1.5);
  clabel(C,H,'FontSize',12,'Color','red')
  xlabel('Time')
  datetick('x','mm/yyyy','keeplimits')
  ylabel('z (m)')
  set(gca,'FontSize',14)
  title([plotVar ' (' plotVarUnits ')'])
  
  subplot 312
  plotVar = 'salt';
  index = find(strcmp(plotVariables, plotVar));
  plotVarUnits = plotVariablesUnits{index};
  plotMap=[];
  for i=1:nRomsData
      plotMap = [plotMap squeeze(ncread(romsDataFile{i},plotVar,[3 3 k100 1],[1 1 Inf Inf]))];
  end
  v=21.6:0.8:35.2;
  [C,H] = contour(datenum(tMap(k100:end,:)),zMap(k100:end,:),plotMap,v,...
      '-k','LineWidth',1.5);
  clabel(C,H,'FontSize',12,'Color','red')
  xlabel('Time')
  datetick('x','mm/yyyy','keeplimits')
  ylabel('z (m)')
  set(gca,'FontSize',14)
  title([plotVar ' (' plotVarUnits ')'])

  subplot 313
  plotVar = 'rho';
  index = find(strcmp(plotVariables, plotVar));
  plotVarUnits = plotVariablesUnits{index};
  plotMap=[];
  for i=1:nRomsData
      plotMap = [plotMap squeeze(ncread(romsDataFile{i},plotVar,[3 3 k100 1],[1 1 Inf Inf]))];
  end
  v=22:0.5:27;
  [C,H] = contour(datenum(tMap(k100:end,:)),zMap(k100:end,:),plotMap,v,...
      '-k','LineWidth',1.5);
  clabel(C,H,'FontSize',12,'Color','red')
  xlabel('Time')
  datetick('x','mm/yyyy','keeplimits')
  ylabel('z (m)')
  set(gca,'FontSize',14)
  title([plotVar ' (' plotVarUnits ')'])
  
  %sgtitle(TITLE)
  set(gcf,'Position',[497   206   934   1012])
  print('TSUpper100m.png','-r200', '-dpng')
  close

% Stratitification and Vertical mixing profiles on January, April, July,
% October

  nProfiles = 4; % January, April, July, October
 
  year0=year(romsTime(1));
  year1=year(romsTime(end));

  profileMonths=[1 4 7 10];

  dateProfiles = [];
  idxProfiles = [];

  for y=year0:year1
      for m=1:3:10

          dateProfiles=[dateProfiles datetime(y,m,15,12,0,0)];

          deltaT = abs(romsTime-datetime(y,m,15,12,0,0));

          idxProfiles=[idxProfiles find(deltaT==min(deltaT),1,'first')];

      end
  end

  idxProfiles=idxProfiles(idxProfiles>1 & idxProfiles<nTimes);

  idxProfiles=reshape(idxProfiles,nProfiles,year1-year0);

  profileRho=zeros(nPlot,3,nProfiles);
  profileRhoZ=zeros(nPlot,nProfiles);
  profileAkt=zeros(nPlot+1,3,nProfiles);

  dZr=zPlot(3:end)-zPlot(1:end-2);
  i=1;
  for j=1:nProfiles

      % Rho
      plotVar='rho';
      map=squeeze(ncread(romsDataFile{i},plotVar,[3 3 1 1],[1 1 Inf Inf]));
      profileRho(:,1,j)=min(map(:,idxProfiles(j,:)),[],2);
      profileRho(:,2,j)=mean(map(:,idxProfiles(j,:)),2);
      profileRho(:,3,j)=max(map(:,idxProfiles(j,:)),[],2);
      tmp=profileRho(3:end,2,j)-profileRho(1:end-2,2,j);
      
      profileRhoZ(2:end-1,j)=tmp./dZr;

      % Akt
      plotVar='AKt';
      map=squeeze(ncread(romsDataFile{i},plotVar,[3 3 1 1],[1 1 Inf Inf]));
      profileAkt(:,1,j)=min(map(:,idxProfiles(j,:)),[],2);
      profileAkt(:,2,j)=mean(map(:,idxProfiles(j,:)),2);
      profileAkt(:,3,j)=max(map(:,idxProfiles(j,:)),[],2);

  end

  kf=kf+1;
  figure(kf) % Rho profiles
  plotVar = 'rho';
  index = find(strcmp(plotVariables, plotVar));
  plotVarUnits = plotVariablesUnits{index};
  titles={'January','April','July','October'};
  for i=1:nProfiles
      subplot (2,2,i)
      ha=plot(profileRho(:,2,i),zPlot,'-r','LineWidth',2);
      % Plot max/min range
      xP=[profileRho(:,1,i);flipud(profileRho(:,3,i))];
      yP=[zPlot;flipud(zPlot)];
      hb=patch(xP,yP,'r','FaceAlpha',.3);
      xlabel([plotVar ' (' plotVarUnits ')'])
      ylabel('z (m)')
      legend('\sigma(z)','d\sigma(z)/dz')
      grid on
      %title(titles{i})
      axis([22 30 -480 0])
      hAx(1)=gca;
      set(hAx(1),'FontSize',14)
      hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','right','color','none');
      hold(hAx(2),'on')
      %plot(hAx(2),x2,y1)
      hc=plot(hAx(2),profileRhoZ(:,i),zPlot,'-c','LineWidth',2);
      set(hAx(2),'FontSize',14)
      hAx(2).YLim=[-480 0];
      hAx(2).YTick=[];
      hAx(2).XLabel.String='d\sigma(z)/dz (kg m^{-4})';
      legend([ha,hb,hc], '\sigma(z)','min/max \sigma(z)','d\sigma(z)/dz',...
          'Location','southwest');
  end
  set(gcf,'Position',[ 680   251   607   726])
  print('rhoProfiles.png','-r200', '-dpng')
  close

  kf=kf+1;
  figure(kf) % Rho profiles
  plotVar = 'AKt';
  index = find(strcmp(plotVariables, plotVar));
  plotVarUnits = plotVariablesUnits{index};
  titles={'January','April','July','October'};
  for i=1:nProfiles
      subplot (2,2,i)
      plot(profileAkt(:,2,i)*10^4,wPlot,'-r','LineWidth',2)
      xP=[profileAkt(:,1,i)*10^4;flipud(profileAkt(:,3,i)*10^4)];
      yP=[wPlot;flipud(wPlot)];
      patch(xP,yP,'r','FaceAlpha',.3)
      xlabel(['10^4 ' plotVar ' (' plotVarUnits ')'])
      ylabel('z (m)')
      grid on
      title(titles{i})
      axis([0 4 -480 0])
      set(gca,'FontSize',14)
  end
  set(gcf,'Position',[ 680   251   607   726])
  print('aktProfiles.png','-r200', '-dpng')
  close

% Nutrient profiles in on January, April, July,
% October
  
  profileNO3=zeros(nPlot,3,nProfiles);
  profileNH4=zeros(nPlot,3,nProfiles);

  i=1;
  for j=1:nProfiles

      % NO3
      plotVar='NO3';
      map=squeeze(ncread(romsDataFile{i},plotVar,[3 3 1 1],[1 1 Inf Inf]));
      profileNO3(:,1,j)=min(map(:,idxProfiles(j,:)),[],2);
      profileNO3(:,2,j)=mean(map(:,idxProfiles(j,:)),2);
      profileNO3(:,3,j)=max(map(:,idxProfiles(j,:)),[],2);

      % NH4
      plotVar='NH4';
      map=squeeze(ncread(romsDataFile{i},plotVar,[3 3 1 1],[1 1 Inf Inf]));
      profileNH4(:,1,j)=min(map(:,idxProfiles(j,:)),[],2);
      profileNH4(:,2,j)=mean(map(:,idxProfiles(j,:)),2);
      profileNH4(:,3,j)=max(map(:,idxProfiles(j,:)),[],2);

  end

  kf=kf+1;
  figure(kf) % NO3,NH4 profiles
  plotVar = 'NO3';
  index = find(strcmp(plotVariables, plotVar));
  plotVarUnits = plotVariablesUnits{index};
  titles={'Jan.','Apr.','Jul.','Oct.'};
  for i=1:nProfiles
      subplot (2,2,i)
      ha=plot(profileNO3(:,2,i),zPlot,'-r','LineWidth',2);
      % Plot max/min range
      xP=[profileNO3(:,1,i);flipud(profileNO3(:,3,i))];
      yP=[zPlot;flipud(zPlot)];
      hb=patch(xP,yP,'r','FaceAlpha',.3);
      xlabel([plotVar ' ' titles{i}  ' (' plotVarUnits ')'])
      ylabel('z (m)')
      legend('NO3','min/max NO3')
      grid on
      %title(titles{i})
      axis([0 140 -480 0])
      hAx(1)=gca;
      set(hAx(1),'FontSize',14)
      hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','right','color','none');
      hold(hAx(2),'on')
      %plot(hAx(2),x2,y1)
      hc=plot(hAx(2),profileNH4(:,2,i),zPlot,'-c','LineWidth',2);
      hold on
      % Plot max/min range
      xP=[profileNH4(:,1,i);flipud(profileNH4(:,3,i))];
      yP=[zPlot;flipud(zPlot)];
      hd=patch(xP,yP,'c','FaceAlpha',.3);
      set(hAx(2),'FontSize',14)
      hAx(2).XLim=[0 5];
      hAx(2).YLim=[-480 0];
      hAx(2).YTick=[];
      hAx(2).XLabel.String=['NH4 ' titles{i}   ' (' plotVarUnits ')'];
      legend([ha,hb,hc,hd], 'NO3','min/max NO3','NH4','min/max NH4',...
          'Location','southeast');
  end
  set(gcf,'Position',[ 680   251   607   726])
  print('no3Profiles.png','-r200', '-dpng')
  close

  

