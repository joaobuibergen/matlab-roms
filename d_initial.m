%
%  D_INITIAL:  Driver script to create a ROMS initial conditions file.
%
%  This a user modifiable script that can be used to prepare ROMS initial
%  conditions NetCDF file.  It sets-up all the necessary parameters and
%  variables. USERS can use this as a prototype for their application.
%

% svn $Id: d_initial.m 1102 2022-01-06 21:50:44Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%  Set input/output NetCDF files.

 my_root = '/Users/joao/SEAS/Models/Masfjorden/GS2019114_M24';

 GRDname = fullfile(my_root, 'mf_1d_grd.nc');
 OAname  = fullfile(my_root, 'Data/OA',      'oa4_lev94_feb.nc');
%INIname = fullfile(my_root, 'Data/netcdf3', 'damee4_levfeb_b.nc');

 INIname = 'm24_ini_gs2019114.nc';

%  Set local variables.

 OA_INTERPOLATE = 0;               % switch to interpolate from OA fields
%method='linear';                  % linear interpolation
 method='spline';                  % spline interpolation

%--------------------------------------------------------------------------
%  Set application parameters in structure array, S.
%--------------------------------------------------------------------------

%  Initial conditions output file name.

S.ncname = INIname;

%  Spherical grid switch.
%
%            [0] Cartesian grid
%            [1] Spherical grid

S.spherical = 1;

%  Number of interior RHO-points (Lm, Mm) in the X- and Y-directions. These
%  are the values specified in ROMS input script.

S.Lm = 4;
S.Mm = 4;

%  Number of vertical levels (N) at RHO-points (center grid cell).

S.N = 120;

%  Number of total active and passive tracers. Usually, NT=2 for potential
%  temperature and salinity.

S.NT = 16; %(Temp, Salt and 14 bio tracers)

%  Set vertical stretching parameters. These values are specified in ROMS
%  input script.

%  Vertical transfomation equation (Vtransform).
%    (Check  https://www.myroms.org/wiki/index.php/Vertical_S-coordinate)
%
%            [1] Original tranformation equation
%            [2] New transformation equation

S.Vtransform = 2;

%  Vertical stretching function (Vstretching).
%    (Check  https://www.myroms.org/wiki/index.php/Vertical_S-coordinate)
%
%            [1] Original function (Song and Haidvogel, 1994)
%            [2] A. Shchepetkin function
%            [3] R. Geyer function

S.Vstretching = 4;

%  Vertical stretching parameters:
%
%            tetha_s:  S-coordinate surface control parameter.
%            tetha_b:  S-coordinate bottom control parameter.
%            Tcline:   S-coordinate surface/bottom stretching width (m)

S.theta_s = 7.0;
S.theta_b = 6.0;
S.Tcline  = 200.0;

S.hc = S.Tcline;

%--------------------------------------------------------------------------
%  Create initial condition Netcdf file.
%--------------------------------------------------------------------------

[~]=c_initial(S);

[~]=c_initial_bio_fennel(S);

%  Set attributes for "ocean_time".

avalue='days since 2000-01-01 00:00:00';
[~]=nc_attadd(INIname,'units',avalue,'ocean_time');
  
avalue='360.0 days in every year';
%[~]=nc_attadd(INIname,'calendar',avalue,'ocean_time');

%--------------------------------------------------------------------------
%  Set grid variables.
%--------------------------------------------------------------------------

V=nc_vnames(GRDname);
nvars=length(V.Variables);

%  Horizontal grid variables. Read in for input GRID NetCDF file.

if (S.spherical),
  S.lon_rho = nc_read(GRDname, 'lon_rho');
  S.lat_rho = nc_read(GRDname, 'lat_rho');
  
  S.lon_u   = nc_read(GRDname, 'lon_u');
  S.lat_u   = nc_read(GRDname, 'lat_u');
  
  S.lon_v   = nc_read(GRDname, 'lon_v');
  S.lat_v   = nc_read(GRDname, 'lat_v');
else  
  S.x_rho   = nc_read(GRDname, 'x_rho');
  S.y_rho   = nc_read(GRDname, 'y_rho');
  
  S.x_u     = nc_read(GRDname, 'x_u');
  S.y_u     = nc_read(GRDname, 'y_u');
  
  S.x_v     = nc_read(GRDname, 'x_v');
  S.y_v     = nc_read(GRDname, 'y_v');  
end  

%  Read in Land/Sea mask, if appropriate.

for n=1:nvars,
  name=char(V.Variables(n).Name);
  switch (name),
    case 'mask_rho'
      S.mask_rho = nc_read(GRDname, 'mask_rho');
    case 'mask_u'
      S.mask_u   = nc_read(GRDname, 'mask_u');
    case 'mask_v'
      S.mask_v   = nc_read(GRDname, 'mask_v');
  end,
end,


%  Bathymetry.

S.h = nc_read(GRDname, 'h');

%  Set vertical grid variables.

[S.s_rho, S.Cs_r]=stretching(S.Vstretching, ...
                             S.theta_s, S.theta_b, S.hc, S.N,           ...
			     0, 1);

[S.s_w,   S.Cs_w]=stretching(S.Vstretching, ...
                             S.theta_s, S.theta_b, S.hc, S.N,           ...
			     1, 1);

%--------------------------------------------------------------------------
%  Set zero initial conditions.
%--------------------------------------------------------------------------

Lr = S.Lm+2;   Lu = Lr-1;   Lv = Lr;
Mr = S.Mm+2;   Mu = Mr;     Mv = Mr-1;

S.zeta = zeros([Lr Mr]);
S.ubar = zeros([Lu Mu]);
S.vbar = zeros([Lv Mv]);
S.u    = zeros([Lu Mu S.N]);
S.v    = zeros([Lv Mv S.N]);
S.temp = zeros([Lr Mr S.N]);
S.salt = zeros([Lr Mr S.N]);

%  If Land/Sea masking arrays are not found, initialize them to unity.

if (~isfield(S, 'mask_rho')),  S.mask_rho = ones([Lr Mr]);  end,
if (~isfield(S, 'mask_u'  )),  S.mask_u   = ones([Lu Mu]);  end,
if (~isfield(S, 'mask_v'  )),  S.mask_v   = ones([Lv Mv]);  end,

%--------------------------------------------------------------------------
%  Write out grid variables.
%--------------------------------------------------------------------------
			 
[~]=nc_write(INIname,   'spherical',   S.spherical);

[~]=nc_write(INIname,   'Vtransform',  S.Vtransform);
[~]=nc_write(INIname,   'Vstretching', S.Vstretching);
[~]=nc_write(INIname,   'theta_s',     S.theta_s);
[~]=nc_write(INIname,   'theta_b',     S.theta_b);
[~]=nc_write(INIname,   'Tcline',      S.Tcline);
[~]=nc_write(INIname,   'hc',          S.hc);

[~]=nc_write(INIname,   's_rho',       S.s_rho);
[~]=nc_write(INIname,   's_w',         S.s_w);
[~]=nc_write(INIname,   'Cs_r',        S.Cs_r);
[~]=nc_write(INIname,   'Cs_w',        S.Cs_w);

[~]=nc_write(INIname,   'h',           S.h);

if (S.spherical),
  [~]=nc_write(INIname, 'lon_rho',     S.lon_rho);
  [~]=nc_write(INIname, 'lat_rho',     S.lat_rho);
  [~]=nc_write(INIname, 'lon_u',       S.lon_u);
  [~]=nc_write(INIname, 'lat_u',       S.lat_u);
  [~]=nc_write(INIname, 'lon_v',       S.lon_v);
  [~]=nc_write(INIname, 'lat_v',       S.lat_v);
else
  [~]=nc_write(INIname, 'x_rho',       S.x_rho);
  [~]=nc_write(INIname, 'y_rho',       S.y_rho);
  [~]=nc_write(INIname, 'x_u',         S.x_u);
  [~]=nc_write(INIname, 'y_u',         S.y_u);
  [~]=nc_write(INIname, 'x_v',         S.x_v);
  [~]=nc_write(INIname, 'y_v',         S.y_v);
end

%--------------------------------------------------------------------------
%  Compute depths at horizontal and vertical RHO-points.
%--------------------------------------------------------------------------

igrid = 1;

[z_r] = set_depth(S.Vtransform, S.Vstretching,                          ...
                  S.theta_s, S.theta_b, S.hc, S.N,                      ...
                  igrid, S.h, S.zeta);

figure
plot(1:S.N,squeeze(z_r(1,1,:)))
hold on
plot(1:S.N,squeeze(z_r(1,1,:)),'ro')
xlabel('Sigma level')
ylabel('Depth [m]')


%--------------------------------------------------------------------------
%  Interpolate OA of temperature and salinity from standard levels to
%  model depths.
%--------------------------------------------------------------------------

if (OA_INTERPOLATE),

  disp(' ')
  disp('Interpolating from OA fields, please wait ...');
  
  InpRec = 1;

  Zoa=nc_read(OAname, 'zout');

  oa_temp=nc_read(OAname, 'temp', InpRec);
  oa_salt=nc_read(OAname, 'salt', InpRec);

  for j=1:Mr,
    for i=1:Lr,
      Zroms = squeeze(z_r(i,j,:));
      Toa   = squeeze(oa_temp(i,j,:));
      Soa   = squeeze(oa_salt(i,j,:));

      S.temp(i,j,:) = interp1(Zoa, Toa, Zroms, method);
      S.salt(i,j,:) = interp1(Zoa, Soa, Zroms, method);
    end,
  end,
  
end,

%---------------------------------------------------------------------------
%  Compute initial conditions for biological variables.
%---------------------------------------------------------------------------

stationLat = 60+52.19047/60;
stationLon = 5+21.9086/60;

bottleDataFile = 'GS2019114_M24.csv';

T = readtable(bottleDataFile);

depthBottle = flipud(T.Depth);

tic = flipud(T.DIC__mol_kg_);
alk = flipud(T.Alk__mol_kg_);
no3 = flipud(T.NO3__mol_l_);

% Read temperature, salinity and density from calibrated CTD data  

T2=readtable('/Volumes/Ctd_Data_Vestlandsfjordar/CTD/Calibrated/CTD_2019_GS2019114.csv');

indexM24=T2.station==382;

stationCTDTemperature=T2.temperature_conservative_C_(indexM24);
stationCTDPressure=T2.pressure_dbar_(indexM24);
stationCTDSalinity=T2.salinity_practical_(indexM24);
stationCTDDensityAnomaly=T2.sigma0_kg_m3_(indexM24);
stationCTDOxygen=T2.oxygen_mL_L_(indexM24)-1000.;
  
% Get depth from pressure to later interpolate to ROMS depths
  
stationCTDDepth = sw_dpth(stationCTDPressure,stationLat);

% Read temperature, salinity and density from calibrated CTD data  

load('/Users/joao/SEAS/Data/Masfjorden/CTD-data/GEOSARS_2019/CTD/GS2019_Calibrated.mat');

stationCTDIndex = 3;


stationCTDPressure2=CTD(stationCTDIndex).P(1:end-1)';
stationCTDFluorescence=CTD(stationCTDIndex).flC(1:end-1)';

% Get depth from pressure to later interpolate to ROMS depths

stationCTDDepth2 = sw_dpth(stationCTDPressure2,stationLat);

zroms=-1*squeeze(z_r(1,1,:));

% Interpolate Temperature, salinity and oxygen

iniTemp = interp1(stationCTDDepth,stationCTDTemperature,zroms,'linear','extrap');
iniSalt = interp1(stationCTDDepth,stationCTDSalinity,zroms,'linear','extrap');
iniOxygen = interp1(stationCTDDepth,stationCTDOxygen,zroms,'linear','extrap');

figure
subplot 131
plot(iniTemp,zroms)
hold on 
plot(stationCTDTemperature, stationCTDDepth,'.r')
xlabel('Temperature (deg C)')
ylabel('depth (m)')
set(gca,'YDir','reverse')
legend('ROMS I.C.','CTD')
subplot 132
plot(iniSalt,zroms)
hold on 
plot(stationCTDSalinity, stationCTDDepth,'.r')
xlabel('Salinity (PSU)')
set(gca,'YDir','reverse')
legend('ROMS I.C.','CTD')
subplot 133
plot(iniOxygen,zroms)
hold on 
plot(stationCTDOxygen, stationCTDDepth,'.r')
xlabel('Dissolved Oxygen (ml l^{-1})')
set(gca,'YDir','reverse')
legend('ROMS I.C.','CTD')

% Interpolate NO3, TIC, ALK

iniNO3 = interp1(depthBottle,no3,zroms,'linear','extrap');
iniNH4 = 0.05*iniNO3;

iniTIC = interp1(depthBottle,tic,zroms,'linear','extrap');
iniALK = interp1(depthBottle,alk,zroms,'linear','extrap');

figure
subplot 131
plot(iniNO3,zroms)
hold on 
plot(no3, depthBottle,'or')
xlabel('NO3 (mmol l^{-1})')
ylabel('depth (m)')
set(gca,'YDir','reverse')
legend('ROMS I.C.','Sample')
subplot 132
plot(iniTIC,zroms)
hold on 
plot(tic, depthBottle,'or')
xlabel('TIC (mmol kg^{-1})')
set(gca,'YDir','reverse')
legend('ROMS I.C.','Sample')
subplot 133
plot(iniALK,zroms)
hold on 
plot(alk, depthBottle,'or')
xlabel('Alkalinity (mol kg^{-1})')
set(gca,'YDir','reverse')
legend('ROMS I.C.','Sample')

% Estimate remaining biological variables phytoplankton, chlorophyll,
% zooplankton and detritral species

iniChl = interp1(stationCTDDepth2,stationCTDFluorescence,zroms,'linear','extrap');

% Phytoplankton (using the Chl-a/P ration of 1.59; Gutknecht et al. 2013)

iniPhy = iniChl/1.59;

% Zooplankton (using 0.15 ratio of Zooplankton to Phytoplankton biomass)

iniZoo = 0.15 * iniPhy;

iniLdN = 0.0025*ones(size(iniPhy));
iniLdC = iniLdN;

iniSdN = 0.05*ones(size(iniPhy));
iniSdC = iniSdN;

%---------------------------------------------------------------------------
%  Write out initial conditions.
%---------------------------------------------------------------------------

IniRec = 1;                               % NetCDF time record

S.ocean_time = 7214;                % initial conditions time (s)

[~]=nc_write(INIname, 'ocean_time', S.ocean_time, IniRec);

[~]=nc_write(INIname, 'zeta', S.zeta, IniRec);
[~]=nc_write(INIname, 'ubar', S.ubar, IniRec);
[~]=nc_write(INIname, 'vbar', S.vbar, IniRec);
[~]=nc_write(INIname, 'u',    S.u,    IniRec);
[~]=nc_write(INIname, 'v',    S.v,    IniRec);
% [~]=nc_write(INIname, 'temp', S.temp, IniRec);
% [~]=nc_write(INIname, 'salt', S.salt, IniRec);

% Write temp and salt variables

tmp=permute(repmat(iniTemp,1,1,Lr,Mr),[4 3 1 2]);

[~]=nc_write(INIname, 'temp', tmp, IniRec);

tmp=permute(repmat(iniSalt,1,1,Lr,Mr),[4 3 1 2]);

[~]=nc_write(INIname, 'salt', tmp, IniRec);

% Write bio vars

tmp=permute(repmat(iniOxygen,1,1,Lr,Mr),[4 3 1 2]);

[~]=nc_write(INIname, 'oxygen', tmp, IniRec);

tmp=permute(repmat(iniNO3,1,1,Lr,Mr),[4 3 1 2]);

[~]=nc_write(INIname, 'NO3', tmp, IniRec);

tmp=permute(repmat(iniNH4,1,1,Lr,Mr),[4 3 1 2]);

[~]=nc_write(INIname, 'NH4', tmp, IniRec);

tmp=permute(repmat(iniPhy,1,1,Lr,Mr),[4 3 1 2]);

[~]=nc_write(INIname, 'phytoplankton', tmp, IniRec);

tmp=permute(repmat(iniZoo,1,1,Lr,Mr),[4 3 1 2]);

[~]=nc_write(INIname, 'zooplankton', tmp, IniRec);

tmp=permute(repmat(iniChl,1,1,Lr,Mr),[4 3 1 2]);

[~]=nc_write(INIname, 'chlorophyll', tmp, IniRec);

tmp=permute(repmat(iniTIC,1,1,Lr,Mr),[4 3 1 2]);

[~]=nc_write(INIname, 'TIC', tmp, IniRec);

tmp=permute(repmat(iniALK,1,1,Lr,Mr),[4 3 1 2]);

[~]=nc_write(INIname, 'alkalinity', tmp, IniRec);

tmp=permute(repmat(iniLdN,1,1,Lr,Mr),[4 3 1 2]);

[~]=nc_write(INIname, 'LdetritusN', tmp, IniRec);

tmp=permute(repmat(iniLdC,1,1,Lr,Mr),[4 3 1 2]);

[~]=nc_write(INIname, 'LdetritusC', tmp, IniRec);

tmp=permute(repmat(iniSdN,1,1,Lr,Mr),[4 3 1 2]);

[~]=nc_write(INIname, 'SdetritusN', tmp, IniRec);

tmp=permute(repmat(iniSdC,1,1,Lr,Mr),[4 3 1 2]);

[~]=nc_write(INIname, 'SdetritusC', tmp, IniRec);


