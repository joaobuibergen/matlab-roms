% This m-file generates an html report with the results
%
%

% Base classes
import mlreportgen.report.* 
import mlreportgen.dom.* 

% Create the report object
rpt = Report(MyAppTITLE,'html'); 

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

% Diagnostics
ch3 = Chapter; 
ch3.Title = 'Diagnostics'; 
% Solution history
sec6 = Section; 
sec6.Title = 'Basin Dissoved Oxygen';

% Time series of dissolved oxygen and dissolved oxygen consumption rates
image1 = FormalImage();
image1.Image = which('basinMeanO2.png');
para = Paragraph(['Figure 3.1 Time series of basin volume-weighted averages of O2, O2 consumption rate' ...
    'and O2 sediment consumption rate. O2 consumption rate is the sum of zooplankton respiration ' ...
    ' (negligible in the basin), heterotrophic bacterial respiration and nitrification. The volume-weighted ' ...
    'average is computed between bottom and sill depth.']);
image1.Caption = para;
image1.Height = '14.81cm';
image1.Width = '25.75cm';
add(sec6,image1);

% Time series of dissolved oxygen and dissolved oxygen consumption rates
image2 = FormalImage();
image2.Image = which('basinMeanRespNit.png');
para = Paragraph(['Figure 3.2 Time series of basin volume-weighted averages of Zooplankton respiration, ' ...
    'heteretrophic bacteria respiration and nitrification. The volume-weighted ' ...
    'average is computed between bottom and sill depth.']);
image2.Caption = para;
image2.Height = '14.81cm';
image2.Width = '25.75cm';
add(sec6,image2);




append(ch3,sec6)
append(rpt,ch3); 

% Finalize report
close(rpt)
rptview(rpt)

