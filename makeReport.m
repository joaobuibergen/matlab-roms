% This m-file generates an html report with the results
%
%

% Base classes
import mlreportgen.report.* 
import mlreportgen.dom.* 

% Create the report object
rpt = Report(MyAppTITLE,'html-file'); 

% Add a title page
tp = TitlePage; 
tp.Title = MyAppTITLE; 
tp.Subtitle = 'Diagnostics'; 
tp.Author = ''; 
append(rpt,tp);

% Add a table of contents.
append(rpt,TableOfContents); 

% Add report contents
ch1 = Chapter; 
ch1.Title = 'Model configuration'; 
% Run info
sec0 = Section; 
sec0.Title = 'Run info';
para = Paragraph(['Date: ' runDateTime]);
append(sec0,para)
para = Paragraph(['Run directory: ' romsRunDir]);
append(sec0,para)
para = Paragraph(['Source directory: ' romsSourceDir]);
append(sec0,para)
append(ch1,sec0)
% CPP Options
sec1 = Section; 
sec1.Title = 'CPP Options';
para = Paragraph(['CPP Header file: ' [lower(MyAppCPP) '.h']]); 
append(sec1,para)
para = Paragraph(strjoin(['CPP Options: ' strjoin(romsCPP,',')])); 
append(sec1,para)
append(ch1,sec1)
% Files
sec2 = Section; 
sec2.Title = 'Input and output files';
para = Paragraph(['Input: ' inputFile ' ' ...
    'Grid: ' romsGridFile ' '...
    'Initial: ' romsIniFile ' '...
    'Forcing: ' strjoin(romsFrcFile,',') ' '...
    'History: ' strjoin(romsDataFile,',') ' '...
    'Averages: ' strjoin(romsAvgFile,',') ' '...
    'Diagnostics: ' strjoin(romsDiaFile,',') ' '...
    '']);
append(sec2,para)
append(ch1,sec2)
% Grid
sec3 = Section; 
sec3.Title = 'Grid Information';
para = Paragraph(['Interior xi rho-points Lm=' num2str(romsGrid.Lm)]);
append(sec3,para)
para = Paragraph(['Interior eta rho-points Mm=' num2str(romsGrid.Mm)]);
append(sec3,para)
para = Paragraph(['Vertical levels N=' num2str(romsGrid.N)]);
append(sec3,para)
para = Paragraph(['Vtransform=' num2str(romsGrid.Vtransform)]);
append(sec3,para)
para = Paragraph(['Vstretching=' num2str(romsGrid.Vstretching)]);
append(sec3,para)
para = Paragraph(['Surface stretching theta s=' num2str(romsGrid.theta_s)]);
append(sec3,para)
para = Paragraph(['Bottom stretching theta b=' num2str(romsGrid.theta_b)]);
append(sec3,para)
para = Paragraph(['Thermocline depth Tcline=' num2str(romsGrid.Tcline) ' m']);
append(sec3,para)
para = Paragraph(['hc=' num2str(romsGrid.hc) ' m']);
append(sec3,para)
append(ch1,sec3)
% Time steps
sec4 = Section; 
sec4.Title = 'Time steps';
para = Paragraph(['Baroclinic Dt=' num2str(romsDt) ' s']);
append(sec4,para)
para = Paragraph(['Barotropic Dtfast=' num2str(romsDt/romsDtFast) ' s']);
append(sec4,para)
para = Paragraph(['History file output=' num2str(romsNHisOut*romsDt) ' s']);
append(sec4,para)
para = Paragraph(['Averages file output=' num2str(romsNAvgOut*romsDt) ' s']);
append(sec4,para)

append(ch1,sec4)
% Solution history
sec5 = Section; 
sec5.Title = 'Solution history';
para = Paragraph(['Reference time: ' datestr(romsRefTime,'yyyy.mm.dd HH:MM:SS')]);
append(sec5,para)
para = Paragraph(['Solution start date: ' datestr(romsStartTime,'yyyy.mm.dd HH:MM:SS')]);
append(sec5,para)

append(ch1,sec5)
append(rpt,ch1); 

% Fennel model
ch2 = Chapter; 
ch2.Title = 'Fennel model'; 
% 
append(rpt,ch2); 
table1 = BaseTable(magic(5));
table1.Title = "Parameters of Fennel model";
table1.Content=bioParTable;
add(rpt,table1);

% Results
ch3 = Chapter; 
ch3.Title = 'Results'; 
% Solution history
sec6 = Section; 
sec6.Title = 'Hydrography';

% Temperature, Salinity and Density Anomaly
image1 = FormalImage();
image1.Image = which('TSUpper100m.png');
image1.Caption = ['Hydrography in the upper 100 m. Top: temperature; '...
    'Middle: salinity; Bottom: density anomaly.'];
image1.Height = '24.75cm';
image1.Width = '25.5cm';
add(sec6,image1);

% Average January, April, July, October profiles of density anomaly and
% density gradient
image2 = FormalImage();
image2.Image = which('rhoProfiles.png');
image2.Caption = ['Average profiles of density anomaly and density anomaly gradient.'...
    ' Top left: January; Top right: April; Bottom left: July; Bottom right: October'];
image2.Height = '25.6117cm';
image2.Width = '21.4136cm';
add(sec6,image2);

% Average January, April, July, October profiles of vertical mixing
image2 = FormalImage();
image2.Image = which('aktProfiles.png');
image2.Caption = ['Average profiles of vertical diffusivity.'...
    ' Top left: January; Top right: April; Bottom left: July; Bottom right: October'];
image2.Height = '25.6117cm';
image2.Width = '21.4136cm';
add(sec6,image2);

% Average January, April, July, October profiles of NO3 and NH4
image2 = FormalImage();
image2.Image = which('no3Profiles.png');
image2.Caption = ['Average profiles of NO3 and NH4.'...
    ' Top left: January; Top right: April; Bottom left: July; Bottom right: October'];
image2.Height = '25.6117cm';
image2.Width = '21.4136cm';
add(sec6,image2);


append(ch3,sec6)
append(rpt,ch3); 

% Diagnostics
ch4 = Chapter; 
ch4.Title = 'Diagnostics'; 
% Solution history
sec7 = Section; 
sec7.Title = 'Basin Dissolved Oxygen Concentration and Consumption Rate';

para = Paragraph(['Oxygen consumption has two components: nitrification ' ...
    'and respiration. These oxygen sinks are described by the following terms in eq. (1) of Fennel et al. (2013):']);
append(sec7,para)
add(sec7,Equation('-2\hat{n}NH4-R_{O2:NH4}(lZoo+\hat{r}_{SD}SDet+\hat{r}_{LD}LDet)'));% - R_{O2:NH4} \( lZoo+ r_{SD} SDet + r_{LD} LDet \)'));
para = Paragraph(['The first term is the nitrification sink with \hat{n} '...
    'the nitrification flux. The second term is the oxygen due to zooplankton '...
    'and heterotrohic bacteria, with $l$ the zooplankton excretion rate '...
    ' and \hat{r}_{SD} and \hat{r}_{LD} the remineralization rates of '...
    'small and large detritus.']);
append(sec7,para)
% Oxygen in the basin
image2 = FormalImage();
image2.Image = which('basinMeanO2.png');
image2.Caption = ['Basin oxygen. Top: Oxygen concentration; Middle: oxygen '...
'consumption rate; Bottom: Volume averaged oxygen concentration and consumption '...
 'rate.'];
image2.Height = '23.55cm'; 
image2.Width = '23.55cm';
add(sec7,image2);


append(ch4,sec7)
append(rpt,ch4); 

% Finalize report
close(rpt)
%rptview(rpt)

