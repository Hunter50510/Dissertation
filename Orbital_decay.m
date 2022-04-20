%% Orbital decay

fprintf('\n')
warning('off')

% Provides orbit of satellite in Circular orbit, i = 0

% Earth
mu = 398600; %km^3 s^-2
RE = 6371; %km
%g0 = 9.81e-3; %km/s
karmin_line = 100; %km

%Earth's rotation
T_Earth = 86164.0905; %One sidereal day
omega_Earth = 2*pi/T_Earth;
velo_Earth = omega_Earth*RE;

% input
inputs = {'Number of rotations per orbit: (Recommended 1-3)','Length of tether: (Recommended 100-1800)'};
title_input = 'Rotations and Length';
dims = [1 40];
hint = {'1','2500'};
control = inputdlg(inputs,title_input,dims,hint);
rot_orbit = str2double(control{1}); %Number of rotations per orbit
sat_length = str2double(control{2}); %Length of tether
fprintf('Number of rotations = %g\n',rot_orbit)
fprintf('Length of tether = %gkm\n\n',sat_length)
tic
% Sat information

sat_radius = sat_length/2; %km
a = RE + sat_radius + karmin_line; % Distance from centre of the Earth
theta_s = 0; %theta at the start of the run

% Orbit information
T = 2 * pi * sqrt(a^3/mu); % Time period
v_orb = (2*pi*a)/T; % Orbital velocity
delta_theta = 2*pi*rot_orbit; % Total angle change over orbit 

omega = delta_theta/T; % Angular velocity

% Orbital decay

r_sat_cs = 100e-3; %Radius of sat cross section (m)
rho_cbnt = 1600; %Density of carbon nano tubes (kg/m^3)

T_R = T;
rho_R = 5.77e-7;%*(100/sat_length); %http://www.braeunig.us/space/atmos.htm   https://www.ukssdc.ac.uk/cgi-bin/wdcc1/secure/msis_models.pl
v_R = omega_Earth*(RE+karmin_line)/1000;
cd_R = 0.1;
A = pi*r_sat_cs*sat_length*1000;
m = rho_cbnt*pi*(r_sat_cs^2*(.99*r_sat_cs)^2);

alpha_0_R = 0.5*rho_R*v_R^2*cd_R*(A/m);

dR_dt = alpha_0_R*T_R/pi;

fprintf('dR/dt = %gm/s\n',dR_dt)


