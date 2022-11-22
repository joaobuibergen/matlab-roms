% Computes diagnostics of the roms model implementation for fjords. 
% This version of the script applied to 1D water column implementations.
% The diagnostics are discussed in the notebook Simulations/Diagnostics of
% Fjord models.

% J. H. Bettencourt, Bergen, November 2022


%  Revision     Details
% ----------    -------
% 09.11.2022    Adapted from post_process_ROMS_1D.m
% 17.11.2022    Modified to be called from d_postprocessing.m
  
% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % Definitions 
% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 
% % Application title (from ROMS log file)
% 
%   [status,cmdout] = system(strjoin(["gawk 'NR==6 {print $1}' " fullfile(romsLog.folder,romsLog.name)]));
%    
%   MyAppTITLE = cmdout(1:end-1); % Remove newline character
% 
%   disp(['Application title: ' MyAppTITLE])
%   
%   inputFile = [MyAppTITLE '.in'];
% 
%   disp(['Input file: ' inputFile])
% 
%   [status,cmdout] = system(strjoin(["gawk '/Header file/ {print $4}' " fullfile(romsLog.folder,romsLog.name)]));
% 
%   cppHeaderFile = strip(cmdout);
% 
%   disp(['CPP Header file: ' cppHeaderFile])
% 
%   str=split(cppHeaderFile,'.');
% 
%   MyAppCPP = upper(str{1});
% 
% % File to process (his - history file; avg - averages file)
% 
%   fileType = 'his';
%   
% % Define plot variables
% 
%   plotVariables={'temp';...
%                  'salt';...
%                  'rho';...
%                  'NO3';...
%                  'NH4';...
%                  'chlorophyll';...
%                  'phytoplankton';...
%                  'zooplankton';...
%                  'oxygen';...
%                  'LdetritusN';...
%                  'SdetritusN';...
%                  'LdetritusC';...
%                  'SdetritusC';...
%                  'TIC';...
%                  'alkalinity';...
%                  'AKv';...
%                  'AKt';...
%                  'AKs';...
%   };
% 
%   surfaceForcing={'Pair';...
%                   'Tair';...
%                   'Uwind';...
%                   'Vwind';...
%                   'Qair';...
%                   'rain';...
%                   'swrad';...
%                   'lwrad_down';...
%                   'SST'
%               };
%           
%   surfaceFluxes={'shflux';...
%                  'latent';...
%                  'sensible';...
%                  'lwrad';...
%                  'ssflux';...
%                  'EminusP';...
%               };
%           
%   diaBioVariables={'NO3_hadv';...
%                    'NO3_vadv';...
%                    'NO3_xadv';...
%                    'NO3_yadv';...
%                    'NO3_vdiff';...
%                    'NO3_rate';...
%                    'NH4_hadv';...
%                    'NH4_vadv';...
%                    'NH4_xadv';...
%                    'NH4_yadv';...
%                    'NH4_vdiff';...
%                    'NH4_rate';...
%                    'chlorophyll_hadv';...
%                    'chlorophyll_vadv';...
%                    'chlorophyll_xadv';...
%                    'chlorophyll_yadv';...
%                    'chlorophyll_vdiff';...
%                    'chlorophyll_rate';...
% 				   'phytoplankton_hadv';...
%                    'phytoplankton_vadv';...
%                    'phytoplankton_xadv';...
%                    'phytoplankton_yadv';...
%                    'phytoplankton_vdiff';...
% 				   'phytoplankton_rate';...
%                    'zooplankton_hadv';...
%                    'zooplankton_vadv';...
%                    'zooplankton_xadv';...
%                    'zooplankton_yadv';...
%                    'zooplankton_vdiff';...
% 				   'zooplankton_rate';...
%                    'LdetritusN_hadv';...
%                    'LdetritusN_vadv';...
%                    'LdetritusN_xadv';...
%                    'LdetritusN_yadv';...
%                    'LdetritusN_vdiff';...
% 				   'LdetritusN_rate';...
% 				   'SdetritusN_hadv';...
%                    'SdetritusN_vadv';...
%                    'SdetritusN_xadv';...
%                    'SdetritusN_yadv';...
%                    'SdetritusN_vdiff';...
% 				   'SdetritusN_rate';...
%                    'LdetritusC_hadv';...
%                    'LdetritusC_vadv';...
%                    'LdetritusC_xadv';...
%                    'LdetritusC_yadv';...
%                    'LdetritusC_vdiff';...
% 				   'LdetritusC_rate';...
% 				   'SdetritusC_hadv';...
%                    'SdetritusC_vadv';...
%                    'SdetritusC_xadv';...
%                    'SdetritusC_yadv';...
%                    'SdetritusC_vdiff';...
% 				   'SdetritusC_rate';...
% 				   'TIC_hadv';...
%                    'TIC_vadv';...
%                    'TIC_xadv';...
%                    'TIC_yadv';...
%                    'TIC_vdiff';...
% 				   'TIC_rate';...
% 				   'alkalinity_hadv';...
%                    'alkalinity_vadv';...
%                    'alkalinity_xadv';...
%                    'alkalinity_yadv';...
%                    'alkalinity_vdiff';...
% 				   'alkalinity_rate';...
% 				   'oxygen_hadv';...
%                    'oxygen_vadv';...
%                    'oxygen_xadv';...
%                    'oxygen_yadv';...
%                    'oxygen_vdiff';...
% 				   'oxygen_rate';...
%                    'denitrification';...
%                    'CO2_airsea';...
%                    'pCO2';...
%                    'O2_airsea';...
%                    'P_Production';...
%                    'NO3_uptake';...
%                    'SO2C';...
%                    'DO2C';...
%                    'PAR';...
%                    'Nitrif_flux';...
%                    'Zoo_respiration';...
%                    'Bac_respiration';...
%                    'OM_sinking'...
%               };
% 
%   nPlotVariables = length(plotVariables);
%   nSurfaceForcing = length(surfaceForcing);
%   nSurfaceFluxes = length(surfaceFluxes);
%   nDiaBioVariables = length(diaBioVariables);
% 
%   
% % Plot types
% 
%   plotInitial = 1;                    % Plot profiles of initial conditions
%   plotProfiles = 0;                   % Plot profiles at specified time
%   plotMaps = 0;                       % Plot time-depth maps
%   
% % Plot on undisturbed depths?
% 
%   plotZeta0 = 1;
%   
% % Times to plot
% 
%   if plotProfiles
%       
%      profilePlotTimes = datenum(2020,03,05);
%      
%   end
%   
%   if plotMaps
%       
%      mapPlotTimes = -1; % If -1, plot all simulation
%      
%   end
%   
% % Sill depth (m)
% 
%   depthSill = 180;
%   
% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % Get grid and other metadata 
% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%  
% % Data files
%   
%   romsDataFile = strip(split(ncreadatt([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'/',...
%       [ fileType '_file']),','));
% 
%   nRomsData = length(romsDataFile);
% 
% % Get the names of initial, forcing and output files. If multiple files,
% % get their names as a list.
% 
%   romsIniFile = ncreadatt([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'/','ini_file');
%   romsFrcFile = strip(split(ncreadatt([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'/','frc_file_01'),','));
%   nRomsFrc = length(romsFrcFile);
%   romsGridFile = ncreadatt([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'/','grd_file');
%   romsAvgFile = strip(split(ncreadatt([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'/','avg_file'),','));
%   nRomsAvg = length(romsAvgFile);
%   romsDiaFile = strip(split(ncreadatt([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'/','dia_file'),','));
%   nRomsDia = length(romsDiaFile);
% 
% % Get ROMS grid (depths for t=0)
% 
%   romsGrid = get_roms_grid(romsGridFile, romsDataFile{1});
%   
% % Get CPP options, model version, run dir, run date
% 
%   %romsCPP = split(ncreadatt(romsDataFile,'/','CPP_options'),',');%roms_cpplist([lower(MyAppCPP) '.h']);
%   romsCPP = string(roms_cpplist([fullfile(romsLog.folder,MyAppTITLE) '_his.nc']));
%   getStr = split(ncreadatt([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'/','history'),',');
%   romsVersion = [getStr{1} ' ' getStr{2} ' rev. '...
%       ncreadatt([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'/','svn_rev')];
%   
%   runDateTime = [getStr{3} ', ' getStr{4}];
%   
%   romsRunDir = ncreadatt([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'/','header_dir');
%   romsSourceDir = ncreadatt([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'/','code_dir');
%   
% % Get baroclinic time step, barotropic time step and number of time steps between output
%   
%   romsDt = ncread([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'dt');
%   romsDtFast = ncread([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'dtfast');
%   
%   %if strcmp(fileType,'his')
%   
%       romsNHisOut = ncread([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'nHIS');
%   
%   %elseif strcmp(fileType,'avg')
%       
%       romsNAvgOut = ncread([fullfile(romsLog.folder,MyAppTITLE) '_his.nc'],'nAVG');
%       
%   %end
%   
%   romsDOut = romsDt*romsNHisOut; % Interval in secs between output fields
%   
% 
% 
% 
% 
% % Get units of plot variables (if units attribute is not present, assumed 
% % nondimensional variable)
% 
%   plotVariablesUnits={};
% 
%   for i=1:nPlotVariables
%       
%      % Get units attribute
%      try
%          tmpStr = ncreadatt(romsDataFile{1},plotVariables{i},'units');
%      catch
%          tmpStr = '-';
%      end
%      tmpStr = replace(tmpStr,'_',' ');
%      tmpStr = replace(tmpStr,'meter','m');
%      tmpStr = replace(tmpStr,'kilogram','kg');
%      tmpStr = replace(tmpStr,'millimole','mmole');
%      tmpStr = replace(tmpStr,'milliequivalents','meqvs');
%      tmpStr = replace(tmpStr,'milligrams','mg');
%      tmpStr = replace(tmpStr,'nitrogen','N');
%      tmpStr = replace(tmpStr,'oxygen','O');
%      tmpStr = replace(tmpStr,'carbon','C');
%      tmpStr = replace(tmpStr,'chlorophyll','Chl');
%      
%      
%      
%      plotVariablesUnits{i}=tmpStr;
%       
%       
%   end
%   
%   surfaceForcingUnits={};
% 
%   for i=1:nSurfaceForcing
%       
%      % Get units attribute
%      try
%          tmpStr = ncreadatt(romsFrcFile{1},surfaceForcing{i},'units');
%      catch
%          tmpStr = '-';
%      end
%      tmpStr = replace(tmpStr,'second','s');
%      tmpStr = replace(tmpStr,'watt','W');
%      tmpStr = replace(tmpStr,'meter','m');
%      surfaceForcingUnits{i}=tmpStr;
%       
%       
%   end
%   
%   surfaceFluxesUnits={};
% 
%   for i=1:nSurfaceFluxes
%       
%      % Get units attribute
%      try
%          tmpStr = ncreadatt(romsDataFile{1},surfaceFluxes{i},'units');
%      catch
%          tmpStr = '-';
%      end
%      tmpStr = replace(tmpStr,'second','s');
%      tmpStr = replace(tmpStr,'watt','W');
%      tmpStr = replace(tmpStr,'meter','m');
%      surfaceFluxesUnits{i}=tmpStr;
%       
%       
%   end
%   
%   diaBioVariablesUnits={};
% 
%   for i=1:nDiaBioVariables
%       
%      % Get units attribute
%      try
%          tmpStr = ncreadatt(romsDiaFile{1},diaBioVariables{i},'units');
%      catch
%          tmpStr = '-';
%      end
%      tmpStr = replace(tmpStr,'_',' ');
%      tmpStr = replace(tmpStr,'meter','m');
%      tmpStr = replace(tmpStr,'kilogram','kg');
%      tmpStr = replace(tmpStr,'second','s');
%      tmpStr = replace(tmpStr,'millimole','mmole');
%      tmpStr = replace(tmpStr,'milliequivalents','meqvs');
%      tmpStr = replace(tmpStr,'milligrams','mg');
%      tmpStr = replace(tmpStr,'nitrogen','N');
%      tmpStr = replace(tmpStr,'oxygen','O');
%      tmpStr = replace(tmpStr,'carbon','C');
%      tmpStr = replace(tmpStr,'chlorophyll','Chl');
%      diaBioVariablesUnits{i}=tmpStr;
%       
%       
%   end
%   
% % Get forcing time in matlab format
% 
%   frcTimeNative=[];
% 
%   for i=1:nRomsFrc
% 
%       frcTimeNative = [frcTimeNative ncread(romsFrcFile{i},'frc_time')'];
% 
%   end
% 
%   frcTimeUnits = ncreadatt(romsFrcFile{1},'frc_time','units');
%   
%   getStr = split(frcTimeUnits);
%   
%   frcRefTime = datetime([getStr{3} ' ' getStr{4}]);
%   
%   if strcmp(getStr{1},'days')
%   
%       frcTime = frcRefTime + days(frcTimeNative);
%   
%   elseif strcmp(getStr{1},'seconds')
%       
%       frcTime = frcRefTime + seconds(frcTimeNative);
%       
%   end
%   
%   nFrcTimes = length(frcTimeNative);
%   
% % Get solution time in matlab format
% 
%   romsTimeNative=[];
% 
%   for i=1:nRomsData
% 
%       romsTimeNative = [romsTimeNative ncread(romsDataFile{i},'ocean_time')'];
% 
%   end
% 
%   romsTimeUnits = ncreadatt(romsDataFile{1},'ocean_time','units');
%   
%   getStr = split(romsTimeUnits);
%   
%   romsRefTime = datetime([getStr{3} ' ' getStr{4}]);
%   
%   romsTime = romsRefTime + seconds(romsTimeNative);
%   
%   nTimes = length(romsTimeNative);
%   
% % Get diagnostics time in matlab format
% 
%   diaTimeNative=[];
% 
%   for i=1:nRomsDia
% 
%       diaTimeNative = [diaTimeNative ncread(romsDiaFile{i},'ocean_time')'];
% 
%   end
% 
%   diaTimeUnits = ncreadatt(romsDiaFile{1},'ocean_time','units');
%   
%   getStr = split(diaTimeUnits);
%   
%   diaRefTime = datetime([getStr{3} ' ' getStr{4}]);
%   
%   diaTime = diaRefTime + seconds(diaTimeNative);
%   
%   nDiaTimes = length(diaTimeNative);
% 
% 
% % ROMS start date
% %strjoin(["gawk '/Header file/ {print $4}' " fullfile(romsLog.folder,romsLog.name)])
%   [status,cmdout] = system(strjoin(["gawk '/DSTART =/ {print $3}'  " inputFile]));
%   
%   tmp=split(cmdout(1:end-1),'d');
%   day0=str2num(tmp{1});
%   ex0=str2num(tmp{2});
%   startDay=day0*10^ex0;
%   romsStartTime = romsRefTime + days(startDay);
% 
% % Print info
% 
% 
% 
% % TO BE DONE
% 
% % % Write relevant data to output file
% % 
% %   fprintf(outFID,'ROMSGridFile = %s\n',romsGridFile);
% %   fprintf(outFID,'ROMSIniFile = %s\n',romsIniFile);
% % 
% %   testOut = romsDataFile';
% %   fprintf(outFID,'ROMSDataFile = %s\n',testOut{:});%romsDataFile);
% %   
% %   testOut = romsFrcFile';
% %   fprintf(outFID,'ROMSFrcFile = %s\n',testOut{:});%romsFrcFile);
% % 
% %   fprintf(outFID,'ROMSVersion = %s\n',romsVersion);
% %   
% %   tmpStr=[];
% %   spx = {' '};
% %   
% %   for i=1:length(romsCPP)
% %       
% %       tmpStr=strcat(tmpStr,spx,strip(romsCPP(i)));
% %       
% %   end
% %   
% %   fprintf(outFID,'CPPOptions = %s\n',strip(tmpStr));
% %   fprintf(outFID,'RUNDir = %s\n',romsRunDir);
% %   fprintf(outFID,'SOURCEDir = %s\n',romsSourceDir);
% %   fprintf(outFID,'RUNDate = %s\n',runDateTime);
% %   fprintf(outFID,'STARTDate = %s\n',datestr(romsTime(1),'yyyy/mm/dd HH:MM:SS'));
% %   fprintf(outFID,'ENDDate = %s\n',datestr(romsTime(end),'yyyy/mm/dd HH:MM:SS'));
% %   fprintf(outFID,'OUTStep = %s\n',num2str(romsDOut/3600));
% %   fprintf(outFID,'DT = %s\n',num2str(romsDt));
% %   fprintf(outFID,'DTFAST = %s\n',num2str(romsDtFast));
% %   fprintf(outFID, 'NVERT = %s\n',num2str(romsGrid.N));
% %   fprintf(outFID, 'VTRANSF = %s\n',num2str(romsGrid.Vtransform));
% %   fprintf(outFID, 'VSTRETCH = %s\n',num2str(romsGrid.Vstretching));
% %   fprintf(outFID, 'THETAS = %s\n',num2str(romsGrid.theta_s));
% %   fprintf(outFID, 'THETAB = %s\n',num2str(romsGrid.theta_b));
% %   fprintf(outFID, 'TCLINE = %s\n',num2str(romsGrid.Tcline));
% %   fprintf(outFID, 'HC = %s\n',num2str(romsGrid.hc));
%   
%   
%   
%   
% % Get biogeochemical model parameters
% 
% bioParVariables={'AttSW';...
%                  'AttChl';...
%                  'PARfrac'
%                  'Vp0';...
%                  'I_thNH4';...
%                  'D_p5NH4';...
%                  'NitriR';...
%                  'K_NO3';...
%                  'K_NH4';...
%                  'K_PO4';...
%                  'K_Phy';...
%                  'Chl2C_m';...
%                  'ChlMin';...
%                  'PhyCN';...
%                  'R_P2N';...
%                  'PhyIP';...
%                  'PhyIS';...
%                  'PhyMin';...
%                  'PhyMR';...
%                  'ZooAE_N';...
%                  'ZooBM';...
%                  'ZooCN';...
%                  'ZooER';...
%                  'ZooGR';...
%                  'ZooMin';...
%                  'ZooMR';...
%                  'LDeRRN';...
%                  'LDeRRC';...
%                  'CoagR';...
%                  'SDeRRN';...
%                  'SDeRRC';...
%                  'RDeRRN';...
%                  'wPhy';...
%                  'wLDet';...
%                  'wSDet';...
%                  'pCO2air';...
%                  };
%   bioParValues={};
%   bioParUnits={};
%   bioParLongNames={};
%   
% % Write header in bio paramater file
% 
% %   fprintf(bioParFID,'Parameter, Symbol, Value, Units\n');
%   
%   for k=1:length(bioParVariables)
% % 	double AttSW ;
% % 		AttSW:long_name = "light attenuation by seawater" ;
% % 		AttSW:units = "meter-1" ;
%     bioParValues{k} = ncread(romsDataFile{1},bioParVariables{k});
%     % Have to do some string replacements to shorten variable units 
%     bioParLongNames{k} = replace(ncreadatt(romsDataFile{1},bioParVariables{k},'long_name'),',','');
%     bioParLongNames{k} = replace(bioParLongNames{k},'phytoplankton','phyto');
%     bioParLongNames{k} = replace(bioParLongNames{k},'zooplankton','zoo');
%     bioParLongNames{k} = replace(bioParLongNames{k},'carbon','C');
%     bioParLongNames{k} = replace(bioParLongNames{k},'Carbon','C');
%     bioParLongNames{k} = replace(bioParLongNames{k},'nitrogen','N');
%     bioParLongNames{k} = replace(bioParLongNames{k},'Nitrogen','N');
%     try 
%       bioParUnits{k} = replace(ncreadatt(romsDataFile{1},bioParVariables{k},'units'),'_',' ');
%       bioParUnits{k} = replace(bioParUnits{k},'meter','m');
%       bioParUnits{k} = replace(bioParUnits{k},'milligram','mg');
%       bioParUnits{k} = replace(bioParUnits{k},'parts per million','ppm');
%       bioParUnits{k} = replace(bioParUnits{k},'milli','m');
%       bioParUnits{k} = replace(bioParUnits{k},'watt','W');     
%     catch
%       bioParUnits{k} = '-';
%     end
%     
% %     % Write to output file
% %     outStr=[bioParLongNames{k} ',' ...
% %         replace(bioParVariables{k},'_','\_') ',' num2str(bioParValues{k}) ',' bioParUnits{k}];
% %     fprintf(bioParFID,'%s\n',outStr);
%     
%   end
% 
%   Values=cell2mat(bioParValues)';
%   Description=char(bioParLongNames);
%   Units=char(bioParUnits);
%   Params=char(bioParVariables);
%   bioParTable = table(Values,Units,Description,'RowNames',bioParVariables)
%   
%   
% %   % Close output files
% % 
% %   fclose(outFID)
% %   fclose(bioParFID);
% 
% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % Aux stuff
% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 
%   sec2day = 3600*24; % Seconds to days
% 
%   cellArea = 1./((romsGrid.pm).*(romsGrid.pn));
% 
%   nPlot = romsGrid.N;
%   zPlot = squeeze(romsGrid.z_r(3,3,:));
%   kBasin = find(zPlot<-1*depthSill,1,'last');
% 
%   cellThickness = diff(romsGrid.z_w,1,3);
%   cellVolume = repmat(cellArea,[1 1 romsGrid.N]).*cellThickness;
% 
%   basinVolume = sum(sum(sum(cellVolume(:,:,1:kBasin))));
% 
% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % Plot parameters
% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 
%   posMap3 = [680 240 730 740];
%   posMap4 = [680 190 730 790];
%   posMap5 = [680  30 730 950];
% 
% % Figure counter
% 
%   kf = 0;

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Compute and plot diagnostics
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%
% 1. Depth of euphotic layer
%

  plotVar = 'PAR';
  index = find(strcmp(diaBioVariables, plotVar));
  plotVarUnits = diaBioVariablesUnits{index};
  varData = [];
  for i=1:nRomsDia
      varData = [varData squeeze(ncread(romsDiaFile{i},plotVar,[3 3 1 1],[1 1 Inf Inf]))'];
  end
  
  ke=[]; % Index where I(ke) is 1% of I0
  de=[]; % Depth at k = ke
  eVolume=[]; % Volume of the euphotic layer
  slVolume=[]; % Volume of the surface layer below the euphotic layer
  
  for i=1:length(varData)

      ke=[ke find(varData(i,:)<=0.01*varData(i,end),1,'last')];
      de=[de romsGrid.z_r(3,3,ke(end))];

      eVolume=[eVolume sum(sum(sum(cellVolume(:,:,ke:romsGrid.N))))];
      slVolume=[slVolume sum(sum(sum(cellVolume(:,:,kBasin:ke))))];

  end

%
% 1. Basin Mean O2
%
  
  elO2 = [];
  slO2 = [];
  baO2 = [];

  for i=1:nRomsData

      dataInfo = ncinfo(romsDataFile{i});

      for j=1:dataInfo.Dimensions(14).Length
      
      %o2Basin = ncread(romsDataFile{i},'oxygen',[1 1 1 1],[Inf Inf kBasin Inf]);
      
      %o2BasinMean = [o2BasinMean squeeze(sum(sum(sum(o2Basin.*repmat(cellVolume(:,:,1:kBasin),[1 1 1 length(o2Basin)])))))/basinVolume];

          tmp = ncread(romsDataFile{i},'oxygen',[1 1 1 j],[Inf Inf Inf 1]);
    
          elO2 = [elO2 squeeze(sum(sum(sum(tmp(:,:,ke(min(j,length(ke))):romsGrid.N,:).*...
              cellVolume(:,:,ke(min(j,length(ke))):romsGrid.N)))))/eVolume(min(j,length(ke)))];
    
          slO2 = [slO2 squeeze(sum(sum(sum(tmp(:,:,kBasin:ke(min(j,length(ke))),:).*...
              cellVolume(:,:,kBasin:ke(min(j,length(ke))))))))/slVolume(min(j,length(ke)))];
    
          baO2 = [baO2 squeeze(sum(sum(sum(tmp(:,:,1:kBasin,:).*...
              cellVolume(:,:,1:kBasin)))))/basinVolume];

      end


  end

%
% 2. Diffusive flux of O2 into basin
%

  o2BasinVDiffMap = zeros([size(romsGrid.pm) length(diaTime)]);

  l0=1;
  for i=1:nRomsDia

      diaInfo = ncinfo(romsDiaFile{i});

      hZ = repmat(romsGrid.Hz(:,:,kBasin),[1 1 diaInfo.Dimensions(13).Length]);

      l1=l0+diaInfo.Dimensions(13).Length-1;
      o2BasinVDiffMap(:,:,l0:l1) = squeeze(ncread(romsDiaFile{i},'oxygen_vdiff',[1 1 kBasin 1],...
          [Inf Inf 1 Inf])).*hZ;

      l0=l1+1;

  end

  o2BasinVDiff=squeeze(mean(mean(o2BasinVDiffMap,1),2))*sec2day;


%
% 3. Vertical DO flux in to the basin from the sediment layer
%

  sedO2ConMap = zeros([size(romsGrid.pm) length(diaTime)]);

  l0=1;
  for i=1:nRomsDia

      diaInfo = ncinfo(romsDiaFile{i});

      l1=l0+diaInfo.Dimensions(13).Length-1;
      sedO2ConMap(:,:,l0:l1) = ncread(romsDiaFile{i},'SO2C');

      l0=l1+1;

  end

  sedO2Con=squeeze(mean(mean(sedO2ConMap,1),2));

%
% 4. Dissolved oxygen consumption rate in the basin
%

%   o2BasinConMean = [];
% 
%   for i=1:nRomsDia
%       
%       o2BasinCon = ncread(romsDiaFile{i},'DO2C',[1 1 1 1],[Inf Inf kBasin Inf]);
%       
%       o2BasinConMean = [o2BasinConMean squeeze(sum(sum(sum(o2BasinCon.*repmat(cellVolume(:,:,1:kBasin),[1 1 1 length(o2BasinCon)])))))/basinVolume];
% 
%   end

  elO2C = [];
  slO2C = [];
  baO2C = [];

  for i=1:nRomsDia

      diaInfo = ncinfo(romsDiaFile{i});

      for j=1:diaInfo.Dimensions(13).Length

          tmp = ncread(romsDiaFile{i},'DO2C',[1 1 1 j],[Inf Inf Inf 1]);

          elO2C = [elO2C squeeze(sum(sum(sum(tmp(:,:,ke(min(j,length(ke))):romsGrid.N,:).*...
              cellVolume(:,:,ke(min(j,length(ke))):romsGrid.N)))))/eVolume(min(j,length(ke)))];

          slO2C = [slO2C squeeze(sum(sum(sum(tmp(:,:,kBasin:ke(min(j,length(ke))),:).*...
              cellVolume(:,:,kBasin:ke(min(j,length(ke))))))))/slVolume(min(j,length(ke)))];

          baO2C = [baO2C squeeze(sum(sum(sum(tmp(:,:,1:kBasin,:).*...
              cellVolume(:,:,1:kBasin)))))/basinVolume];

      end

  end

%
% 5. Primary production (EL)
%

  primaryProdRate = [];

  for i=1:nRomsDia

      diaInfo = ncinfo(romsDiaFile{i});

      for j=1:diaInfo.Dimensions(13).Length

          primaryProd = ncread(romsDiaFile{i},'P_Production',[1 1 ke(j) j],[Inf Inf Inf 1]);

          primaryProdRate = [primaryProdRate squeeze(sum(sum(sum(primaryProd.*...
              cellVolume(:,:,ke(j):romsGrid.N)))))/eVolume(j)];

      end

  end

%
% 6. Nutrient level (EL/SL/BA)
%

  elNO3 = [];
  slNO3 = [];
  baNO3 = [];

  elNH4 = [];
  slNH4 = [];
  baNH4 = [];
     
  for i=1:nRomsData

      dataInfo = ncinfo(romsDataFile{i});

      for j=1:dataInfo.Dimensions(14).Length

          tmp = ncread(romsDataFile{i},'NO3',[1 1 1 j],[Inf Inf Inf 1]);

          elNO3 = [elNO3 squeeze(sum(sum(sum(tmp(:,:,ke(min(j,length(ke))):romsGrid.N,:).*...
              cellVolume(:,:,ke(min(j,length(ke))):romsGrid.N)))))/eVolume(min(j,length(ke)))];

          slNO3 = [slNO3 squeeze(sum(sum(sum(tmp(:,:,kBasin:ke(min(j,length(ke))),:).*...
              cellVolume(:,:,kBasin:ke(min(j,length(ke))))))))/eVolume(min(j,length(ke)))];

          baNO3 = [baNO3 squeeze(sum(sum(sum(tmp(:,:,1:kBasin,:).*...
              cellVolume(:,:,1:kBasin)))))/eVolume(min(j,length(ke)))];

          tmp = ncread(romsDataFile{i},'NH4',[1 1 1 j],[Inf Inf Inf 1]);

          elNH4 = [elNH4 squeeze(sum(sum(sum(tmp(:,:,ke(min(j,length(ke))):romsGrid.N,:).*...
              cellVolume(:,:,ke(min(j,length(ke))):romsGrid.N)))))/eVolume(min(j,length(ke)))];
          
          slNH4 = [slNH4 squeeze(sum(sum(sum(tmp(:,:,kBasin:ke(min(j,length(ke))),:).*...
              cellVolume(:,:,kBasin:ke(min(j,length(ke))))))))/eVolume(min(j,length(ke)))];

          baNH4 = [baNH4 squeeze(sum(sum(sum(tmp(:,:,1:kBasin,:).*...
              cellVolume(:,:,1:kBasin)))))/eVolume(min(j,length(ke)))];

      end

  end

%
% 7. Freshwater nutrient fluxes
%

   freshwNO3Map = zeros([size(romsGrid.pm) length(romsTime)]);

  l0=1;
  for i=1:nRomsDia

      diaInfo = ncinfo(romsDataFile{i});

      l1=l0+diaInfo.Dimensions(14).Length-1;
      freshwNO3Map(:,:,l0:l1) = ncread(romsDataFile{i},'NO3_sflux');

      l0=l1+1;

  end

  freshwNO3=squeeze(mean(mean(freshwNO3Map,1),2))*sec2day;

%
% 8. Zooplankton respiration (EL/SL/BA)
% 9. Heterotrophic bacteria respiration (EL/SL/BA)
% 10. Nitrification (EL/SL/BA)
% 11. Horizontal flux of DO (EL/SL) (NOT YET IMPLEMENTED!!!)
% 12. Vertical flux of organic matter (EL/SL/BA)
%

  elZoR = [];
  slZoR = [];
  baZoR = [];

  elBaR = [];
  slBaR = [];
  baBaR = [];

  elNit = [];
  slNit = [];
  baNit = [];

  slVOM = [];
  baVOM = [];
  sdVOM = [];

     
  for i=1:nRomsDia

      diaInfo = ncinfo(romsDiaFile{i});

      for j=1:diaInfo.Dimensions(13).Length

          tmp = ncread(romsDiaFile{i},'Zoo_respiration',[1 1 1 j],[Inf Inf Inf 1]);

          elZoR = [elZoR squeeze(sum(sum(sum(tmp(:,:,ke(min(j,length(ke))):romsGrid.N,:).*...
              cellVolume(:,:,ke(min(j,length(ke))):romsGrid.N)))))/eVolume(min(j,length(ke)))];

          slZoR = [slZoR squeeze(sum(sum(sum(tmp(:,:,kBasin:ke(min(j,length(ke))),:).*...
              cellVolume(:,:,kBasin:ke(min(j,length(ke))))))))/slVolume(min(j,length(ke)))];

          baZoR = [baZoR squeeze(sum(sum(sum(tmp(:,:,1:kBasin,:).*...
              cellVolume(:,:,1:kBasin)))))/basinVolume];

          tmp = ncread(romsDiaFile{i},'Bac_respiration',[1 1 1 j],[Inf Inf Inf 1]);

          elBaR = [elBaR squeeze(sum(sum(sum(tmp(:,:,ke(min(j,length(ke))):romsGrid.N,:).*...
              cellVolume(:,:,ke(min(j,length(ke))):romsGrid.N)))))/eVolume(min(j,length(ke)))];

          slBaR = [slBaR squeeze(sum(sum(sum(tmp(:,:,kBasin:ke(min(j,length(ke))),:).*...
              cellVolume(:,:,kBasin:ke(min(j,length(ke))))))))/slVolume(min(j,length(ke)))];

          baBaR = [baBaR squeeze(sum(sum(sum(tmp(:,:,1:kBasin,:).*...
              cellVolume(:,:,1:kBasin)))))/basinVolume];

          tmp = ncread(romsDiaFile{i},'Nitrif_flux',[1 1 1 j],[Inf Inf Inf 1]);

          elNit = [elNit 2.0*squeeze(sum(sum(sum(tmp(:,:,ke(min(j,length(ke))):romsGrid.N,:).*...
              cellVolume(:,:,ke(min(j,length(ke))):romsGrid.N)))))/eVolume(min(j,length(ke)))];

          slNit = [slNit 2.0*squeeze(sum(sum(sum(tmp(:,:,kBasin:ke(min(j,length(ke))),:).*...
              cellVolume(:,:,kBasin:ke(min(j,length(ke))))))))/slVolume(min(j,length(ke)))];

          baNit = [baNit 2.0*squeeze(sum(sum(sum(tmp(:,:,1:kBasin,:).*...
              cellVolume(:,:,1:kBasin)))))/basinVolume];

          tmp = ncread(romsDiaFile{i},'OM_sinking',[1 1 1 j],[Inf Inf Inf 1]);

          slVOM = [slVOM squeeze(sum(sum(sum(tmp(:,:,ke(min(j,length(ke))),:)))))];

          baVOM = [baVOM squeeze(sum(sum(sum(tmp(:,:,kBasin,:)))))];

          sdVOM = [sdVOM squeeze(sum(sum(sum(tmp(:,:,1,:)))))];
         
      end

  end

  kf=kf+1;
  figure(kf)

  %subplot 311
  TITLE = {'Basin vol. average of O_2 concentration and O_2 consumption rate'};
  
  yyaxis left
  plotVar1 = 'O_2 Concentration';
  plotVar1Units = 'mmol O_2 m^{-3}';
  plot(romsTime, baO2,'LineWidth',2)
  xlabel('Time')
  ylabel({plotVar1;['(' plotVar1Units ')']})
  grid on
  
  yyaxis right
  plotVar2 = 'O_2 Consumption';
  plotVar2Units = 'mmol O_2 m^{-3} day^{-1}';
  plot(diaTime,baO2C,'LineWidth',2)
  hold on
  plot(diaTime,sedO2Con,'LineWidth',2)
  ylabel({plotVar2;['(' plotVar2Units ')']})
  set(gca,'FontSize',14,'FontName','Times')
  legend('O_2','O_2 CONS BASIN','O_2 CONS SED')
  title(TITLE,'FontSize',18,'FontName','Times')
  set(gcf,'Units','centimeters','Position', [23.95 19.61 25.75 14.81])
  print('basinMeanO2.png','-r200', '-dpng')
  close
  
  %subplot 312
  kf=kf+1;
  figure(kf)
  TITLE = {'Basin vol. average of Zoop. and Bact. respiration and nitrification'};
  
  yyaxis left
  plotVar1 = 'Respiration';
  plotVar1Units = 'mmol O_2 m^{-3} day^{-1}';
  plot(diaTime,baZoR,'LineWidth',2)
  hold on
  plot(diaTime,baBaR,'LineWidth',2)
  xlabel('Time')
  ylabel({plotVar1;['(' plotVar1Units ')']})
  grid on
  yyaxis right
  plotVar2 = 'Nitrification';
  plotVar2Units = 'mmol O_2 m^{-3} day^{-1}';
  plot(diaTime,baNit,'LineWidth',2)
  ylabel({plotVar2;['(' plotVar2Units ')']})

  legend('Zoop. Res','Het. Bac. Resp.','Nitrification')
  set(gca,'FontSize',14,'FontName','Times')
  title(TITLE,'FontSize',18,'FontName','Times')
  set(gcf,'Units','centimeters','Position', [23.95 19.61 25.75 14.81])
  print('basinMeanRespNit.png','-r200', '-dpng')
  close
  
  subplot 313
  TITLE = {'Vertical flux of O_2 at sill depth and sediment O_2 flux'};
  
  yyaxis left
  plotVar1 = 'Vertical O_2 flux';
  plotVar1Units = 'mmol O_2 m^{-2} day^{-1}';
  plot(diaTime,o2BasinVDiff,'LineWidth',2)
  xlabel('Time')
  ylabel({plotVar1;['(' plotVar1Units ')']})
  grid on
  yyaxis right
  plotVar2 = 'Vertical OM flux';
  plotVar2Units = 'mmol N m^{-2} day^{-1}';
  plot(diaTime,baVOM,'LineWidth',2)
  hold on
  plot(diaTime,100*sdVOM,'LineWidth',2)
  ylabel({plotVar2;['('      plotVar2Units ')']})

  legend('Vertical O_2 flux','Vert. OM flux','Vert. OM flux Sed. Interface x 100')
  title(TITLE)
  
  set(gcf,'Position', posMap3)

  print('basinMeanO2Budget.eps','-r200', '-depsc2')
  close


  kf=kf+1;
  figure(kf)

  subplot 311
  TITLE = {'Nutrient input and levels in the Euphotic layer'};
  
  yyaxis left
  plotVar1 = 'Nutrient input';
  plotVar1Units = 'mmol N m^{-2} day^{-1} ';
  plot(romsTime, freshwNO3,'LineWidth',2)
  xlabel('Time')
  ylabel({plotVar1;['(' plotVar1Units ')']})
  grid on
  
  yyaxis right
  plotVar2 = 'Nutrient concentration';
  plotVar2Units = 'mmol N m^{-3}';
  plot(romsTime,elNO3,'LineWidth',2)
  hold on
  plot(romsTime,10*elNH4,'LineWidth',2)
  ylabel({plotVar2;['(' plotVar2Units ')']})

  legend('Flux NO3','NO3','10*NH4')
  title(TITLE)
  
  subplot 312
  TITLE = {'Euphotic layer vol. averages of O_2 concentration and O_2 consumption rate'};
 
  yyaxis left
  plotVar1 = 'O_2 Concentration';
  plotVar1Units = 'mmol O_2 m^{-3}';
  plot(romsTime, elO2,'LineWidth',2)
  xlabel('Time')
  ylabel({plotVar1;['(' plotVar1Units ')']})
  grid on
  
  yyaxis right
  plotVar2 = 'O_2 Consumption';
  plotVar2Units = 'mmol O_2 m^{-3} day^{-1}';
  plot(diaTime,elO2C,'LineWidth',2)
  ylabel({plotVar2;['(' plotVar2Units ')']})

  legend('O_2','O_2 CONS')
  title(TITLE)
  
  subplot 313
  TITLE = {'Euphotic layer vol. averages of Zoop. and Bact. respiration and nitrification'};
  
  yyaxis left
  plotVar1 = 'Respiration';
  plotVar1Units = 'mmol O_2 m^{-3} day^{-1}';
  plot(diaTime,elZoR,'LineWidth',2)
  hold on
  plot(diaTime,elBaR,'LineWidth',2)
  xlabel('Time')
  ylabel({plotVar1;['(' plotVar1Units ')']})
  grid on
  yyaxis right
  plotVar2 = 'Nitrification';
  plotVar2Units = 'mmol O_2 m^{-3} day^{-1}';
  plot(diaTime,elNit,'LineWidth',2)
  ylabel({plotVar2;['(' plotVar2Units ')']})

  legend('Zoop. Res','Het. Bac. Resp.','Nitrification')
  title(TITLE)
  
  set(gcf,'Position', posMap3)

  print('euphoticLayerO2Budget.eps','-r200', '-depsc2')
  close

  
    
