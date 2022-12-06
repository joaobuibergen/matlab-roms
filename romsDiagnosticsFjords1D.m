% Computes diagnostics of the roms model implementation for fjords. 
% This version of the script applied to 1D water column implementations.
% The diagnostics are discussed in the notebook Simulations/Diagnostics of
% Fjord models.

% J. H. Bettencourt, Bergen, November 2022


%  Revision     Details
% ----------    -------
% 09.11.2022    Adapted from post_process_ROMS_1D.m
% 17.11.2022    Modified to be called from d_postprocessing.m
% 24.11.2022    Modified to be run after romsPostprocessFjords1D.m
  


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

      eVolume=[eVolume sum(cellVolume(ke:romsGrid.N))];
      slVolume=[slVolume sum(cellVolume(kBasin:ke))];

  end

%
% 1. Basin Mean O2
%
  
  elO2 = [];
  slO2 = [];
  baO2 = [];
  o2Basin = zeros(kBasin,nTimes);

  l0=1;

  for i=1:nRomsData

      dataInfo = ncinfo(romsDataFile{i});

      l1=l0+dataInfo.Dimensions(14).Length-1;
      
      tmp = squeeze(ncread(romsDataFile{i},'oxygen',[3 3 1 1],[1 1 Inf Inf]));

      o2Basin(:,l0:l1)=tmp(1:kBasin,:);

      l0=l1;

      baO2 = [baO2 squeeze(sum(tmp(1:kBasin,:).*...
              cellVolume(1:kBasin)))/basinVolume];

      for j=1:dataInfo.Dimensions(14).Length

          elO2 = [elO2 squeeze(sum(tmp(ke(min(j,length(ke))):romsGrid.N,j).*...
              cellVolume(ke(min(j,length(ke))):romsGrid.N)))/eVolume(min(j,length(ke)))];
    
          slO2 = [slO2 squeeze(sum(tmp(kBasin:ke(min(j,length(ke))),j).*...
              cellVolume(kBasin:ke(min(j,length(ke))))))/slVolume(min(j,length(ke)))];
    
          

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

  elO2C = [];
  slO2C = [];
  baO2C = [];

  o2CBasin = zeros(kBasin,nTimes);

  l0=1;

  for i=1:nRomsDia

      diaInfo = ncinfo(romsDiaFile{i});

      l1=l0+diaInfo.Dimensions(13).Length-1;
      
      tmp = squeeze(ncread(romsDiaFile{i},'DO2C',[3 3 1 1],[1 1 Inf Inf]));

      o2CBasin(:,l0:l1)=tmp(1:kBasin,:);

      l0=l1;

      baO2C = [baO2C squeeze(sum(tmp(1:kBasin,:).*...
              cellVolume(1:kBasin)))/basinVolume];

      for j=1:diaInfo.Dimensions(13).Length

          elO2C = [elO2C squeeze(sum(tmp(ke(min(j,length(ke))):romsGrid.N,j).*...
              cellVolume(ke(min(j,length(ke))):romsGrid.N)))/eVolume(min(j,length(ke)))];
    
          slO2C = [slO2C squeeze(sum(tmp(kBasin:ke(min(j,length(ke))),j).*...
              cellVolume(kBasin:ke(min(j,length(ke))))))/slVolume(min(j,length(ke)))];

      end

  end

%
% O2 and O2 consumption plots in the basin.
%

  kf=kf+1;
  figure(kf)

  subplot 311
  plotVar = 'oxygen';
  index = find(strcmp(plotVariables, plotVar));
  plotVarUnits = plotVariablesUnits{index};
  v=0:300;
  [C,H] = contour(datenum(tMap(1:kBasin,:)),zMap(1:kBasin,:),o2Basin,...v,...
      '-k','LineWidth',1.5);
  clabel(C,H,'FontSize',12,'Color','red')
  xlabel('Time')
  datetick('x','mm/yyyy','keeplimits')
  ylabel('z (m)')
  set(gca,'FontSize',14)
  title([plotVar ' (' plotVarUnits ')'])

  subplot 312
  plotVar = 'DO2C';
  index = find(strcmp(diaBioVariables, plotVar));
  plotVarUnits = diaBioVariablesUnits{index};
  v=0:300;
  [C,H] = contour(datenum(tMap(1:kBasin,:)),zMap(1:kBasin,:),o2CBasin,...v,...
      '-k','LineWidth',1.5);
  clabel(C,H,'FontSize',12,'Color','red')
  xlabel('Time')
  datetick('x','mm/yyyy','keeplimits')
  ylabel('z (m)')
  set(gca,'FontSize',14)
  ax0=gca;
  xTicks=ax0.XTick;
  xTicksLabels=ax0.XTickLabel;
  xLims=ax0.XLim;
  title([plotVar ' (' plotVarUnits ')'])

  subplot 313
  TITLE = {'Basin vol. average of O_2 concentration and O_2 consumption rate'};
  
  yyaxis left
  plotVar1 = 'O_2 Concentration';
  plotVar1Units = 'mmol O_2 m^{-3}';
  plot(datenum(romsTime), baO2,'-r','LineWidth',2)
  xlabel('Time')
  ylabel({plotVar1;['(' plotVar1Units ')']})
  grid on
  ax=gca;
  ax.YColor = 'red';
  ax.XTick=xTicks;
  ax.XTickLabel=xTicksLabels;
  ax.XLim = xLims;
  
  yyaxis right
  plotVar2 = 'O_2 Consumption';
  plotVar2Units = 'mmol O_2 m^{-3} day^{-1}';
  plot(datenum(diaTime),10*baO2C,'-c','LineWidth',2)
  hold on
  plot(datenum(diaTime),sedO2Con,'--c','LineWidth',2)
  ylabel({plotVar2;['(' plotVar2Units ')']})
  ax=gca;
  ax.YColor = 'cyan';
  ax.YLim = [0 0.5];
  set(gca,'FontSize',14)%,'FontName','Times')
  legend('O_2','O_2 CONS BASIN * 10','O_2 CONS SED','Location','southwest')
  title(TITLE,'FontSize',18,'FontName','Times')
  set(gcf,'Units','centimeters','Position', [23.95 2.99 31.39 31.43])
  print('basinMeanO2.png','-r200', '-dpng')
%   close
%   
%   %subplot 312
%   kf=kf+1;
%   figure(kf)
%   TITLE = {'Basin vol. average of Zoop. and Bact. respiration and nitrification'};
%   
%   yyaxis left
%   plotVar1 = 'Respiration';
%   plotVar1Units = 'mmol O_2 m^{-3} day^{-1}';
%   plot(diaTime,baZoR,'LineWidth',2)
%   hold on
%   plot(diaTime,baBaR,'LineWidth',2)
%   xlabel('Time')
%   ylabel({plotVar1;['(' plotVar1Units ')']})
%   grid on
%   yyaxis right
%   plotVar2 = 'Nitrification';
%   plotVar2Units = 'mmol O_2 m^{-3} day^{-1}';
%   plot(diaTime,baNit,'LineWidth',2)
%   ylabel({plotVar2;['(' plotVar2Units ')']})
% 
%   legend('Zoop. Res','Het. Bac. Resp.','Nitrification')
%   set(gca,'FontSize',14,'FontName','Times')
%   title(TITLE,'FontSize',18,'FontName','Times')
%   set(gcf,'Units','centimeters','Position', [23.95 19.61 25.75 14.81])
%   print('basinMeanRespNit.png','-r200', '-dpng')
%   close
%   
%   subplot 313
%   TITLE = {'Vertical flux of O_2 at sill depth and sediment O_2 flux'};
%   
%   yyaxis left
%   plotVar1 = 'Vertical O_2 flux';
%   plotVar1Units = 'mmol O_2 m^{-2} day^{-1}';
%   plot(diaTime,o2BasinVDiff,'LineWidth',2)
%   xlabel('Time')
%   ylabel({plotVar1;['(' plotVar1Units ')']})
%   grid on
%   yyaxis right
%   plotVar2 = 'Vertical OM flux';
%   plotVar2Units = 'mmol N m^{-2} day^{-1}';
%   plot(diaTime,baVOM,'LineWidth',2)
%   hold on
%   plot(diaTime,100*sdVOM,'LineWidth',2)
%   ylabel({plotVar2;['('      plotVar2Units ')']})
% 
%   legend('Vertical O_2 flux','Vert. OM flux','Vert. OM flux Sed. Interface x 100')
%   title(TITLE)
%   
%   set(gcf,'Position', posMap3)
% 
%   print('basinMeanO2Budget.eps','-r200', '-depsc2')
%   close

%
% 5. Chlorophyll, Primary production and POC
%

  chlA = zeros(romsGrid.N,nTimes);

  l0=1;

  for i=1:nRomsData

      dataInfo = ncinfo(romsDataFile{i});

      l1=l0+dataInfo.Dimensions(14).Length-1;
      
      tmp = squeeze(ncread(romsDataFile{i},'chlorophyll',[3 3 1 1],[1 1 Inf Inf]));

      chlA(:,l0:l1)=tmp;

      l0=l1;

  end

  primaryProdRate = [];



  for i=1:nRomsDia

      diaInfo = ncinfo(romsDiaFile{i});

      primaryProd = squeeze(ncread(romsDiaFile{i},'P_Production',[3 3 1 1],[1 1 Inf Inf]));

      primaryProdRate = squeeze(sum(primaryProd.*cellVolume))/vol100;

  end

  pOC = zeros(romsGrid.N,nTimes);

  l0=1;

  for i=1:nRomsData

      dataInfo = ncinfo(romsDataFile{i});

      l1=l0+dataInfo.Dimensions(14).Length-1;
      
      tmpL = squeeze(ncread(romsDataFile{i},'LdetritusC',[3 3 1 1],[1 1 Inf Inf]));
      tmpS = squeeze(ncread(romsDataFile{i},'SdetritusC',[3 3 1 1],[1 1 Inf Inf]));

      pOC(:,l0:l1)=tmpL+tmpS;

      l0=l1;

  end




  kf=kf+1;
  figure(kf)

  subplot 311
  plotVar = 'chlorophyll';
  index = find(strcmp(plotVariables, plotVar));
  plotVarUnits = plotVariablesUnits{index};
  v1=2:2:10;
  v2=0:0.2:2;
  [C1,H1] = contour(datenum(tMap(k40:nPlot,:)),zMap(k40:nPlot,:),chlA(k40:nPlot,:),v1,...
      '-k','LineWidth',1.5);
  hold on
  [C2,H2] = contour(datenum(tMap(k40:nPlot,:)),zMap(k40:nPlot,:),chlA(k40:nPlot,:),v2,...
      '--k','LineWidth',1);
  clabel(C1,H1,'FontSize',10,'Color','red')
  xlabel('Time')
  datetick('x','mm/yyyy','keeplimits')
  ylabel('z (m)')
  set(gca,'FontSize',14)
  title([plotVar ' (' plotVarUnits ')'])

  subplot 312
  plotVar = 'P_Production';
  index = find(strcmp(diaBioVariables, plotVar));
  plotVarUnits = diaBioVariablesUnits{index};
  v1=0.5;
  v2=[0.001,0.01,0.1];
  [C1,H1] = contour(datenum(tMapDia(k40:nPlot,:)),zMapDia(k40:nPlot,:),primaryProd(k40:nPlot,:),[v1 v1],...
      '-k','LineWidth',1.5);
  hold on
  [C2,H2] = contour(datenum(tMapDia(k40:nPlot,:)),zMapDia(k40:nPlot,:),primaryProd(k40:nPlot,:),v2,...
      '--c','LineWidth',0.5);
  clabel(C1,H1,'FontSize',10,'Color','red')
  %clabel(C2,H2,'FontSize',10,'Color','red')
  xlabel('Time')
  datetick('x','mm/yyyy','keeplimits')
  ylabel('z (m)')
  set(gca,'FontSize',14)
  title([plotVar ' (' plotVarUnits ')'])

  subplot 313
  plotVar = 'LdetritusC';
  index = find(strcmp(plotVariables, plotVar));
  plotVarUnits = plotVariablesUnits{index};
  v1=0:5:40;
  v2=[0.01,0.1,1];
  [C1,H1] = contour(datenum(tMap),zMap,pOC,v1,...
      '-k','LineWidth',1.5);
  hold on
  [C2,H2] = contour(datenum(tMap),zMap,pOC,v2,...
      '--k','LineWidth',0.5);
  clabel(C1,H1,'FontSize',10,'Color','red')
  clabel(C2,H2,'FontSize',10,'Color','red')
  xlabel('Time')
  datetick('x','mm/yyyy','keeplimits')
  ylabel('z (m)')
  set(gca,'FontSize',14)
  title(['POC (' plotVarUnits ')'])


  set(gcf,'Units','centimeters','Position', [23.95 2.99 31.39 31.43])
%
% 6. Nutrient levels (EL/SL/BA)
%

  elNO3 = [];
  slNO3 = [];
  baNO3 = [];

  elNH4 = [];
  slNH4 = [];
  baNH4 = [];
     
  for i=1:nRomsData

      dataInfo = ncinfo(romsDataFile{i});

      tmp = squeeze(ncread(romsDataFile{i},'NO3',[3 3 1 1],[1 1 Inf Inf]));

      for j=1:dataInfo.Dimensions(14).Length

          elNO3 = [elNO3 squeeze(sum(tmp(ke(min(j,length(ke))):romsGrid.N,j).*...
              cellVolume(ke(min(j,length(ke))):romsGrid.N)))/eVolume(min(j,length(ke)))];

          slNO3 = [slNO3 squeeze(sum(tmp(kBasin:ke(min(j,length(ke))),j).*...
              cellVolume(kBasin:ke(min(j,length(ke))))))/slVolume(min(j,length(ke)))];

          baNO3 = [baNO3 squeeze(sum(tmp(1:kBasin)),j).*...
              cellVolume(1:kBasin)/basinVolume];

      end
      
      tmp = squeeze(ncread(romsDataFile{i},'NH4',[3 3 1 1],[1 1 Inf Inf]));

      for j=1:dataInfo.Dimensions(14).Length

          elNH4 = [elNH4 squeeze(sum(tmp(ke(min(j,length(ke))):romsGrid.N,j).*...
              cellVolume(ke(min(j,length(ke))):romsGrid.N)))/eVolume(min(j,length(ke)))];
          
          slNH4 = [slNH4 squeeze(sum(tmp(kBasin:ke(min(j,length(ke))),j).*...
              cellVolume(kBasin:ke(min(j,length(ke))))))/slVolume(min(j,length(ke)))];

          baNH4 = [baNH4 squeeze(sum(tmp(1:kBasin)),j).*...
              cellVolume(1:kBasin)/basinVolume];

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

      tmp = squeeze(ncread(romsDiaFile{i},'Zoo_respiration',[3 3 1 1],[1 1 Inf Inf]));

      for j=1:diaInfo.Dimensions(13).Length

          elZoR = [elZoR squeeze(sum(tmp(ke(min(j,length(ke))):romsGrid.N,j).*...
              cellVolume(ke(min(j,length(ke))):romsGrid.N)))/eVolume(min(j,length(ke)))];

          slZoR = [slZoR squeeze(sum(tmp(kBasin:ke(min(j,length(ke))),j).*...
              cellVolume(kBasin:ke(min(j,length(ke))))))/slVolume(min(j,length(ke)))];

          baZoR = [baZoR squeeze(sum(tmp(1:kBasin)),j).*...
              cellVolume(1:kBasin)/basinVolume];

      end

      tmp = squeeze(ncread(romsDiaFile{i},'Bac_respiration',[3 3 1 1],[1 1 Inf Inf]));

      for j=1:diaInfo.Dimensions(13).Length

          elBaR = [elBaR squeeze(sum(tmp(ke(min(j,length(ke))):romsGrid.N,j).*...
              cellVolume(ke(min(j,length(ke))):romsGrid.N)))/eVolume(min(j,length(ke)))];

          slBaR = [slBaR squeeze(sum(tmp(kBasin:ke(min(j,length(ke))),j).*...
              cellVolume(kBasin:ke(min(j,length(ke))))))/slVolume(min(j,length(ke)))];

          baBaR = [baBaR squeeze(sum(tmp(1:kBasin)),j).*...
              cellVolume(1:kBasin)/basinVolume];

      end

      tmp = squeeze(ncread(romsDiaFile{i},'Nitrif_flux',[3 3 1 1],[1 1 Inf Inf]));

      for j=1:diaInfo.Dimensions(13).Length

          elNit = [elNit 2.0*squeeze(sum(tmp(ke(min(j,length(ke))):romsGrid.N,j).*...
              cellVolume(ke(min(j,length(ke))):romsGrid.N)))/eVolume(min(j,length(ke)))];

          slNit = [slNit 2.0*squeeze(sum(tmp(kBasin:ke(min(j,length(ke))),j).*...
              cellVolume(kBasin:ke(min(j,length(ke))))))/slVolume(min(j,length(ke)))];

          baNit = [baNit 2.0*squeeze(sum(tmp(1:kBasin)),j).*...
              cellVolume(1:kBasin)/basinVolume];

      end

      tmp = squeeze(ncread(romsDiaFile{i},'OM_sinking',[3 3 1 1],[1 1 Inf Inf]));

      for j=1:diaInfo.Dimensions(13).Length

          slVOM = [slVOM squeeze(sum(sum(sum(tmp(ke(min(j,length(ke))),:)))))];

          baVOM = [baVOM squeeze(sum(sum(sum(tmp(kBasin,:)))))];

          sdVOM = [sdVOM squeeze(sum(sum(sum(tmp(1,:)))))];
         
      end

  end

  

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

  
    
