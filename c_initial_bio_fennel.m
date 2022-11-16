function [status]=c_initial_bio_fennel(S)

%
% C_INITIAL_BIO_FENNEL:  Add Fennel model's variables to a ROMS initial conditions NetCDF file
%
% [status]=c_initial_bio_fennel(S)
%
% This function adds the variables of the Fennel model to a ROMS initial conditions NetCDF file using
% specified parameters in structure array, S.
%
% On Input:
%
%    S           Initial condidions creation parameters (structure array):
%
%                  S.ncname           NetCDF file name
%                  S.spherical        Spherical grid switch
%                  S.Vtransform       Vertical transformation equation
%                  S.Lm               Number of interior RHO-points in X
%                  S.Mm               Number of interior RHO-points in Y
%                  S.N                Number of vertical levels
%                  S.NT               Number of active and passive tracers
%
% On Output:
%
%    status      Error flag.
%

% svn $Id: c_initial.m 1105 2022-01-19 02:26:38Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
%  Get initial condition creation parameters.
%--------------------------------------------------------------------------

if (isfield(S,'ncname')),
  ncname=S.ncname;
else
  error(['C_INITIAL - Cannot find dimension parameter: ncname, ',       ...
         'in structure array S']);
end

if (isfield(S,'spherical')),
  spherical=S.spherical;
else
  spherical=0;
end

if (isfield(S,'Vtransform')),
  Vtransform=S.Vtransform;
else
  error(['C_INITIAL - Cannot find dimension parameter: Vtransform, ',   ...
         'in structure array S']);
end

if (isfield(S,'Lm')),
  Lp=S.Lm+2;
else
  error(['C_INITIAL - Cannot find dimension parameter: Lm, ',           ...
         'in structure array S']);
end

if (isfield(S,'Mm')),
  Mp=S.Mm+2;
else
  error(['C_INITIAL - Cannot find dimension parameter: Mm, ',           ...
         'in structure array S']);
end

if (isfield(S,'N')),
  N=S.N;
else
  error(['C_INITIAL - Cannot find dimension parameter: N, ',            ...
         'in structure array S']);
end

if (isfield(S,'NT')),
  NT=S.NT;
else
  error(['C_INITIAL - Cannot find dimension parameter: NT, ',           ...
         'in structure S']);
end

% %--------------------------------------------------------------------------
% %  Create dimensions.
% %--------------------------------------------------------------------------
% 
Dname = {'xi_rho',  'xi_u',  'xi_v',  'xi_psi',                       ...
         'eta_rho', 'eta_u', 'eta_v', 'eta_psi',                      ...
         's_rho', 's_w', 'tracer', 'ocean_time'};
% 
% for value = Dname
%   dim = char(value);
%   switch (dim)
%     case {'xi_rho', 'xi_v'}
%       Dsize.(dim) = Lp;
%     case {'xi_psi', 'xi_u'}
%       Dsize.(dim) = Lp-1;
%     case {'eta_rho', 'eta_u'}
%       Dsize.(dim) = Mp;
%     case {'eta_psi', 'eta_v'}
%       Dsize.(dim) = Mp-1;
%     case {'s_rho'}
%       Dsize.(dim) = N;
%     case {'s_w'}
%       Dsize.(dim) = N+1;
%     case {'tracer'}
%       Dsize.(dim) = NT;
%     case {'ocean_time'}
%       Dsize.(dim) = netcdf.getConstant('UNLIMITED');
%   end
% end
% 
%--------------------------------------------------------------------------
%  Set grid Variables.
%--------------------------------------------------------------------------

grd_vars = {'spherical', 'Vtransform', 'Vstretching',                   ...
            'theta_s', 'theta_b', 'Tcline', 'hc',                       ...
            's_rho', 's_w', 'Cs_r', 'Cs_w', 'h'};
	    
if (spherical)
  grd_vars = [grd_vars, 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi',     ...
                        'lon_u',   'lat_u',   'lon_v',   'lat_v'];
else
  grd_vars = [grd_vars, 'x_rho', 'y_rho', 'x_psi', 'y_psi',             ...
                        'x_u',   'y_u',   'x_v',   'y_v'];
end

grd_vars = [grd_vars, 'mask_rho', 'mask_psi', 'mask_u', 'mask_v'];

%--------------------------------------------------------------------------
%  Set Fennel model variables.
%--------------------------------------------------------------------------

bio_vars = {'ocean_time','NO3', 'NH4', 'phytoplankton', 'chlorophyll', 'zooplankton', ... 
              'LdetritusN','LdetritusC', 'SdetritusN', 'SdetritusC', ...
              'TIC', 'alkalinity', 'oxygen'};

%--------------------------------------------------------------------------
%  Open initial conditions NetCDF file.
%--------------------------------------------------------------------------

%mode = netcdf.getConstant('NOCLOBBER');
%mode = bitor(mode, netcdf.getConstant('64BIT_OFFSET'));

ncid = netcdf.open(ncname,'WRITE');
netcdf.reDef(ncid);

%--------------------------------------------------------------------------
%  Get dimensions.
%--------------------------------------------------------------------------

for value = Dname,
  field = char(value);
  did.(field) = netcdf.inqDimID(ncid, field);
end

%--------------------------------------------------------------------------
%  Define NetCDF variables.
%--------------------------------------------------------------------------

vars = [grd_vars, bio_vars];

for value = vars
  field = char(value);
  Vname.(field) = field;
end

%--------------------------------------------------------------------------
%  Define initial conditions biological variables.
%--------------------------------------------------------------------------

Var.name            = Vname.NO3;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time did.s_rho did.eta_rho did.xi_rho];
Var.long_name       = 'nitrate concentration';
Var.units           = 'millimole_nitrogen meter-3';
Var.time            = Vname.ocean_time;
if (spherical)
  Var.coordinates   = strcat([Vname.lon_rho,' ',Vname.lat_rho,' ',      ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
else
  Var.coordinates   = strcat([Vname.x_rho,' ',Vname.y_rho,' ',          ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.NH4;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time did.s_rho did.eta_rho did.xi_rho];
Var.long_name       = 'ammonium concentration';
Var.units           = 'millimole_nitrogen meter-3';
Var.time            = Vname.ocean_time;
if (spherical)
  Var.coordinates   = strcat([Vname.lon_rho,' ',Vname.lat_rho,' ',      ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
else
  Var.coordinates   = strcat([Vname.x_rho,' ',Vname.y_rho,' ',          ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.phytoplankton;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time did.s_rho did.eta_rho did.xi_rho];
Var.long_name       = 'phytoplankton concentration';
Var.units           = 'millimole_nitrogen meter-3';
Var.time            = Vname.ocean_time;
if (spherical)
  Var.coordinates   = strcat([Vname.lon_rho,' ',Vname.lat_rho,' ',      ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
else
  Var.coordinates   = strcat([Vname.x_rho,' ',Vname.y_rho,' ',          ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.chlorophyll;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time did.s_rho did.eta_rho did.xi_rho];
Var.long_name       = 'chlorophyll concentration';
Var.units           = 'milligrams_chlorophyll meter-3';
Var.time            = Vname.ocean_time;
if (spherical)
  Var.coordinates   = strcat([Vname.lon_rho,' ',Vname.lat_rho,' ',      ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
else
  Var.coordinates   = strcat([Vname.x_rho,' ',Vname.y_rho,' ',          ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.zooplankton;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time did.s_rho did.eta_rho did.xi_rho];
Var.long_name       = 'zooplankton concentration';
Var.units           = 'millimole_nitrogen meter-3';
Var.time            = Vname.ocean_time;
if (spherical)
  Var.coordinates   = strcat([Vname.lon_rho,' ',Vname.lat_rho,' ',      ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
else
  Var.coordinates   = strcat([Vname.x_rho,' ',Vname.y_rho,' ',          ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.LdetritusN;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time did.s_rho did.eta_rho did.xi_rho];
Var.long_name       = 'large nitrogen-detritus concentration';
Var.units           = 'millimole_nitrogen meter-3';
Var.time            = Vname.ocean_time;
if (spherical)
  Var.coordinates   = strcat([Vname.lon_rho,' ',Vname.lat_rho,' ',      ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
else
  Var.coordinates   = strcat([Vname.x_rho,' ',Vname.y_rho,' ',          ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.SdetritusN;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time did.s_rho did.eta_rho did.xi_rho];
Var.long_name       = 'small nitrogen-detritus concentration';
Var.units           = 'millimole_nitrogen meter-3';
Var.time            = Vname.ocean_time;
if (spherical)
  Var.coordinates   = strcat([Vname.lon_rho,' ',Vname.lat_rho,' ',      ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
else
  Var.coordinates   = strcat([Vname.x_rho,' ',Vname.y_rho,' ',          ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.LdetritusC;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time did.s_rho did.eta_rho did.xi_rho];
Var.long_name       = 'large carbon-detritus concentration';
Var.units           = 'millimole_carbon meter-3';
Var.time            = Vname.ocean_time;
if (spherical)
  Var.coordinates   = strcat([Vname.lon_rho,' ',Vname.lat_rho,' ',      ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
else
  Var.coordinates   = strcat([Vname.x_rho,' ',Vname.y_rho,' ',          ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.SdetritusC;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time did.s_rho did.eta_rho did.xi_rho];
Var.long_name       = 'small carbon-detritus concentration';
Var.units           = 'millimole_carbon meter-3';
Var.time            = Vname.ocean_time;
if (spherical)
  Var.coordinates   = strcat([Vname.lon_rho,' ',Vname.lat_rho,' ',      ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
else
  Var.coordinates   = strcat([Vname.x_rho,' ',Vname.y_rho,' ',          ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.TIC;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time did.s_rho did.eta_rho did.xi_rho];
Var.long_name       = 'total inorganic carbon';
Var.units           = 'millimole_carbon meter-3';
Var.time            = Vname.ocean_time;
if (spherical)
  Var.coordinates   = strcat([Vname.lon_rho,' ',Vname.lat_rho,' ',      ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
else
  Var.coordinates   = strcat([Vname.x_rho,' ',Vname.y_rho,' ',          ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.oxygen;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time did.s_rho did.eta_rho did.xi_rho];
Var.long_name       = 'dissolved oxygen concentration';
Var.units           = 'millimole_oxygen meter-3';
Var.time            = Vname.ocean_time;
if (spherical)
  Var.coordinates   = strcat([Vname.lon_rho,' ',Vname.lat_rho,' ',      ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
else
  Var.coordinates   = strcat([Vname.x_rho,' ',Vname.y_rho,' ',          ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

Var.name            = Vname.alkalinity;
Var.type            = nc_constant('nc_double');
Var.dimid           = [did.ocean_time did.s_rho did.eta_rho did.xi_rho];
Var.long_name       = 'total alkalinity';
Var.units           = 'milliequivalents meter-3';
Var.time            = Vname.ocean_time;
if (spherical)
  Var.coordinates   = strcat([Vname.lon_rho,' ',Vname.lat_rho,' ',      ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
else
  Var.coordinates   = strcat([Vname.x_rho,' ',Vname.y_rho,' ',          ...
                              Vname.s_rho,' ',Vname.ocean_time]); 
end
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end
clear Var

%--------------------------------------------------------------------------
%  Update global attributes.
%--------------------------------------------------------------------------

varid  = netcdf.getConstant('nc_global');

% type='ROMS INITIAL file';
% netcdf.putAtt(ncid, varid, 'type', type);

history = netcdf.getAtt(ncid, varid, 'history');

history=[history ...
    '; Initial bio variables using Matlab script: c_initial_bio_fennel, ',date_stamp];
netcdf.putAtt(ncid, varid, 'history', history);

%--------------------------------------------------------------------------
%  Leave definition mode and close NetCDF file.
%--------------------------------------------------------------------------

netcdf.endDef(ncid);
netcdf.close(ncid);

return

