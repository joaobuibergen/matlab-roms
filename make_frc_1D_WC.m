% Script to write the forcing file for the 1d simulations of
% coupled physical-biogeochemistry in Masfjorden, at station M26, with the 
% Fennel (2006) model . The script uses the existing forcing
% file of the ROMS test case BIO TOY, updating it with the data from the 
% norwegian Met Institute forecast model archive, using the opendat
% protocol.
%
%
% Joao Bettencourt, Bergen, August 2022
%
% Revision          Details
% --------          -------
%

  clear variables
  %startup
  addpath /Users/joao/software/matlab2/netcdf
  addpath /Users/joao/software/matlab2/utility
  
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Definitions (KB2020603 M26 Station (from mat file with CTD data))
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% M26 position

%   stationLongitude = 5.41666666666667; 
%   stationLatitude = 60.8721666666667;
  stationLatitude = 60+52.19047/60;
  stationLongitude = 5+21.9086/60;
  
% Station date: '04-Mar-2020 09:10:48'. Rounding to nearest hour in the past

  stationObservationTime = datetime(2019,10,1,0,0,0);
  %stationObservationTime = datetime(2020,12,1,12,0,0);
  
% Advance/delay in hours relative to station date/ run end date

  advanceDate = 0;

% Simulation length in years (assume 1 year = 365 days);

  simulationLength = 3;
  
% Dates for data retrieval

  dataStartDate = stationObservationTime - advanceDate/24;
  yearStartDate = datetime(year(dataStartDate),1,1);
  
  endSimulationDate = dataStartDate + years(simulationLength);
  yearEndDate = datetime(year(endSimulationDate),12,31);

  simulationLengthHours = yearEndDate - yearStartDate;
  
  dataRetrievalDates = yearStartDate + hours(0:1:hours(simulationLengthHours));
  
  nRetrievalDates = length(dataRetrievalDates);

% Build date lists per year

  years = year(yearStartDate):1:year(yearEndDate);
  months = 1:12;
  
  getDates = {};
  forcYear = [];
  forcMonth = [];

  k=1;
  for i=1:length(years)
      for j=1:length(months)

          dateSelect = year(dataRetrievalDates)==years(i) & month(dataRetrievalDates)==months(j);
          getDates{k}=dataRetrievalDates(dateSelect);
          forcYear = [forcYear years(i)];
          forcMonth = [forcMonth months(j)];
          k=k+1;

      end
  end

  spherical = false;
  nctype    = 'nc_double'; 
  Unlimited = true;  

  mode = netcdf.getConstant('CLOBBER');                    % overwrite!!!
  mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET'));

% Output file name

  romsForcing = 'mf_m24_frc_';
  refDateForcing = datetime('2001-1-1 00:00:00');
  
% OpenDAP URL base
% https://thredds.met.no/thredds/dodsC/metpparchive/2020/03/04/met_analysis_1_0km_nordic_20200304T00Z.nc

  opendapUrlBase = 'https://thredds.met.no/thredds/dodsC/metpparchive/';
  productName = '/met_analysis_1_0km_nordic_';
  
% Variables to retrieve

  varNamesROMS={'Pair';...
                'Uwind';...
                'Vwind';...
                'Tair';...
                'Qair';...
                'rain';...
                'swrad';...
                'lwrad_down';... % supplying downward lw flux instead of net
                'cloud';...
                'SST';...
                };
            
 
  varNamesMETNO={'air_pressure_at_sea_level';... % Pa ---> mbar
                'wind_direction_10m';...         % Degrees from ()
                'wind_speed_10m';...             % m/s
                'air_temperature_2m';...         % Kelvin ---> Celsius 
                'relative_humidity_2m';...       % [0-1] --> [0-100]%
                'precipitation_amount';...       % kg/m2 --> kg/m2/s
 'integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time';... % W s/m2 ---> W/m2
                '';...
                'cloud_area_fraction';...        % [0-1]
                '';...
                'x_wind_10m';...
                'y_wind_10m';...
                };
            
  nVariables = length(varNamesROMS);
  nVariablesMet = length(varNamesMETNO);
  
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Get MET.NO data 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  dataUrl=[opendapUrlBase datestr(dataRetrievalDates(1),'yyyy/mm/dd') ...
      productName datestr(dataRetrievalDates(1),'yyyymmddTHHZ') '.nc'];
%
% Get lon/lat grids to select grid point nearest to station
%

  lonMet = ncread(dataUrl,'longitude');   
  latMet = ncread(dataUrl,'latitude');
  
  lonSta = ones(size(lonMet))*stationLongitude;
  latSta = ones(size(latMet))*stationLatitude;
  
% Distances in degrees between grid points and station

  [Dist,Az]=distance(latMet,lonMet,latSta,lonSta);
  
  [Imin,Jmin] = find(Dist==min(min(Dist))); % Indices of grid point closest to station
  
  clear lonMet latMet lonSta latSta

%  
%% Cycle time  
%  

  for k=13:length(getDates)

      tmpData = zeros(nVariablesMet,length(getDates{k}));

      naData = [];

      retrievalDates=getDates{k};

      for i=1:length(retrievalDates)

          disp([' ** Data retrieval: ' datestr(retrievalDates(i),'yyyymmddTHHZ')])

          dataUrl=[opendapUrlBase datestr(retrievalDates(i),'yyyy/mm/dd') ...
              productName datestr(retrievalDates(i),'yyyymmddTHHZ') '.nc'];

          disp([' **   URL: ' dataUrl])

          try
              ncid = netcdf.open(dataUrl);
          catch
              disp('Data not available - File error')
              naData=[naData i];
              tmpData(:,i) = tmpData(:,max(1,i-1));
              continue
          end

          for j=1:nVariablesMet

              if ~isempty(varNamesMETNO{j})

                  try
                      varid = netcdf.inqVarID(ncid,varNamesMETNO{j});
                  catch
                      disp(['Data not available -  Variable error: ' ...
                          varNamesMETNO{j}])
                      tmpData(j,i) = realmax;
                      naData=[naData i];
                      continue
                  end

                  try
                      tmpData(j,i) = netcdf.getVar(ncid,varid,[Imin-1 Jmin-1 0],[1 1 1]);
                      %theData(j) = netcdf.getVar(ncid,varid,[Imin-1 Jmin-1 0],[1 1 1]);
                  catch
                      disp('Data not available - Read error')
                      naData=[naData i];
                      tmpData(j,i) = tmpData(j,max(1,i-1));
                      continue
                  end

              else

                  tmpData(j,i) = realmax('single');


              end

          end % End variable loop

          netcdf.close(ncid);

      end

      % Apply transformations/conversions to MET data


      % Air pressure: Pascal to milibar

      tmpData(1,:) = 0.01*tmpData(1,:);

      % Wind speed and direction to northward, eastward wind velocities

      windF = tmpData(2,:)==realmax;

      tmpU = - tmpData(3,~windF).*sin(pi/180*tmpData(2,~windF));
      tmpV = - tmpData(3,~windF).*cos(pi/180*tmpData(2,~windF));

      tmpData(2,:) = tmpData(11,:); % store u wind in wind direction line
      tmpData(3,:) = tmpData(12,:); % store v wind in wind speed line

      tmpData(2,~windF) = tmpU; % store u wind in wind direction line
      tmpData(3,~windF) = tmpV; % store v wind in wind speed line

      % Air temp: Kelvin to Celsius

      tmpData(4,:) = tmpData(4,:) - 273.15;

      % Relative humidity: [0-1] to [0-100%]

      tmpData(5,:) = tmpData(5,:) * 100;

      % Precipitation amount: kg/m2 --> kg/m2/s
      % I'm just going to divide by 3600s

      tmpData(6,:) = tmpData(6,:) / 3600;

      % Shortwave flux: W s/m2 --> W/m2
      % I'm just going to divide by 3600s

      tmpData(7,:) = tmpData(7,:) / 3600;

      % Compute downward longwave flux from Idso and Jackson (1969) and
      % Konzelmann et al. (1994), taken from Duarte et al (2006)

      sigma = 5.67051*10^-8; % Steffan-Boltzmann constant
      T0 = 273;
      T4 = (tmpData(4,:)+273.15).^4;
      C4 = tmpData(9,:).^4;

      eac = 1-0.261*exp(-7.77*10^-4*(T0-tmpData(4,:)).^2);

      lwDownClearSky = sigma*eac.*T4;

      cloudCorrection = (1-C4);
      lwDown = lwDownClearSky .* cloudCorrection + 0.952*sigma*C4.*T4;


      

      % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      % Create forcing file
      % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      %
      % Most of the code in this section is inspired by d_ecmwf2roms.m

%       forcYear = year(retrievalDates(1));
%       forcMonth = month(retrievalDates(1));

      romsForcingNc = [romsForcing num2str(forcYear(k)) '_' num2str(forcMonth(k),'%02d') ...
          '.nc'];

      S.Filename = romsForcingNc;

      S.Attributes(1).Name      = 'type';
      S.Attributes(1).Value     = 'FORCING file';

      S.Attributes(2).Name      = 'title';
      S.Attributes(2).Value     = ['Met.no Reanalysis 1.0 km, ',           ...
          'Masfjorden 1D '];

      S.Attributes(3).Name      = 'history';
      S.Attributes(3).Value     = ['Forcing file created with ',            ...
          which(mfilename) ' on ' date_stamp];

      %   S.Dimensions(1).Name      = 'lon';
      %   S.Dimensions(1).Length    = Im;
      %   S.Dimensions(1).Unlimited = false;
      %
      %   S.Dimensions(2).Name      = 'lat';
      %   S.Dimensions(2).Length    = Jm;
      %   S.Dimensions(2).Unlimited = false;

      S.Dimensions(1).Name      = 'frc_time';
      S.Dimensions(1).Length    = nc_constant('nc_unlimited');
      S.Dimensions(1).Unlimited = true;

      S.Variables(1) = roms_metadata('frc_time', [], [], Unlimited);

      % Edit the time variable 'units" attribute for the correct reference
      % time and add calendar attribute.

      natts = length(S.Variables(1).Attributes);
      iatt  = strcmp({S.Variables(1).Attributes.Name}, 'units');
      S.Variables(1).Attributes(iatt).Value = ['days since ' datestr(refDateForcing,31)];

      S.Variables(1).Attributes(natts+1).Name  = 'calendar';
      S.Variables(1).Attributes(natts+1).Value = 'gregorian';

      for i=1:nVariables

          S.Variables(1+i) = roms_metadata(varNamesROMS{i}, spherical, nctype, Unlimited);

          % Since this is an 1-D application, the forcing fields only have 1
          % dimension (time), so we need to remove the spatial dimensions
          % created by roms_metadata()

          S.Variables(1+i).Dimensions(1)=[];
          S.Variables(1+i).Dimensions(1)=[];

          % And we change the time dimension name to 'frc_time'

          S.Variables(1+i).Dimensions(1).Name = 'frc_time';

      end

      % Check ROMS metadata structure.  Fill unassigned fields.

      S = check_metadata(S);

      % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      % Write data to forcing file
      % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ncid = nc_create(S.Filename, mode, S); % create a new NetCDF file

      % Write times

      %ncwrite(romsForcingNc,'frc_time',days(dataRetrievalDates-refDateForcing));
      status = nc_write(S.Filename, 'frc_time',  days(retrievalDates-refDateForcing));
      % Loop roms variables list and write fields

      for i=1:nVariables

          if i==8
              field = lwDown;
          else
              field = tmpData(i,:);
          end

          status = nc_write(S.Filename, S.Variables(i+1).Name,  field);

      end


  end


