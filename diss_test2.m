function diss_test2
fprintf('\n')

% Provides orbit of satellite in Circular orbit, i = 0

% Earth
mu = 398600; %km^3 s^-2
RE = 6371; %km
%g0 = 9.81e-3; %km/s
karmin_line = 100; %km

% input
inputs = {'Number of rotations per orbit: (Recommended 1-3)','Length of tether: (Recommended 100-1800)'};
title_input = 'Rotations and Length';
dims = [1 40];
hint = {'1','1000'};
control = inputdlg(inputs,title_input,dims,hint);
rot_orbit = str2double(control{1}); %Number of rotations per orbit
sat_length = str2double(control{2}); %Length of tether
fprintf('Number of rotations = %g\n',rot_orbit)
fprintf('Length of tether = %gkm\n\n',sat_length)

% Sat information

alt = sat_length/2; %km
a = RE + alt + karmin_line; % Distance from centre of the Earth
theta_s = pi/4; %theta at the start of the run

% Orbit information
T = 2 * pi * sqrt(a^3/mu); % Time period
v_orb = (2*pi*a)/T; % Orbital velocity

% Cartesian Coordinates
% As the satellite is orbiting the equator, we will look in the [x y]
% coordinates
e_hat_r = [cos(theta_s); sin(theta_s)]; % Unit vector
r_ = a*e_hat_r; % Postions vector [X Y 0]
e_hat_theta = [-sin(theta_s); cos(theta_s)];
v_ = v_orb*e_hat_theta; % Velocity vector [x y 0]

% Initial conditions
initial_conditions = [r_(1) r_(2) v_(1) v_(2)];
U = initial_conditions;

[~,U] = ode45(@loving_life, [0 T], U);
x = U(:,1); y = U(:,2); u = U(:,3); v = U(:,4);

total__ = ones(length(x));
total = sum(total__);

total_store = [0 0];

for count = 2:total
    x_close = initial_conditions(1)-x(count);
	
    y_close = initial_conditions(2)-y(count);
  
    total_store(count-1) = sqrt(x_close^2  + y_close);
end
close_diff = min(total_store);

% The distance between initial and final point
% diff = sqrt(x(1)^2+y(1)^2+z(1)^2) - sqrt(x(end)^2+y(end)^2+z(end)^2);

fprintf('The difference between the inital and final position of the trajectory is %gkm\n\n',abs(close_diff))

tic
% In order to find the results to the desired accuracy, the ode45 is run
% again, with the accuracy of the run being refined each time
fprintf('Please wait while I run some calculations,\nyour buisness is important to us... *♫Hold Music♫*\n\n')
i = 0;
i2 = 0;
accuracy = 0.001;
while abs(close_diff) > accuracy
    i = i - 0.001;
    i2 = i;
    options = odeset('abstol',1*10^i,'reltol',1*10^i2);

    U2 = initial_conditions;

    [~,U2] = ode45(@loving_life, [0 T], U2, options);
    x = U2(:,1); y = U2(:,2); u = U2(:,3); v = U2(:,4);

    
total__ = ones(length(x));
total = sum(total__);
total_store = [0 0];
for count = 2:total
    x_close = initial_conditions(1)-x(count);
   
    y_close = initial_conditions(2)-y(count);

    total_store(count-1) = sqrt(x_close^2  + y_close);
end
close_diff = min(total_store);
end

fprintf('Thank you for waiting \n\n')
if i2 < 0
    fprintf('After adjusting the accuracy of ode45, the difference between the inital and final position of the trajectory is %0.2gkm (%0.4gkm)\n',close_diff,close_diff)
    fprintf('while the values of abstol and reltol are: 1e%g, and 1e%g\n\n',i,i2)
else
    fprintf('no\n')
end

figure(1)
plot(x,y,'color',[61/255 217/255 201/255])

xlabel('X (km)'),ylabel('Y (km)'),zlabel('Z (km)'),title('The trajectory of the spacecraft around Earth with recalculated error')
grid on
hold on
[X,Y,Z] = sphere;
surf(X*RE,Y*RE,Z*RE)
axis equal
toc

% Rotation of tip relative to surface
delta_theta = 2*pi*rot_orbit; % Total angle change (rad)
omega = delta_theta/T; %Angular velocity = delta theta/time period (rad/s)
tic
sat_COG = [x y];

P1 = zeros(size(x));
P2 = zeros(size(x));
i3 = 0;
sample_size = length(x);
P_angle = theta_s;
P_hold = zeros(size(x));
P_dist_hold = zeros(size(x));
P1_dist = zeros(size(x));
P2_dist = zeros(size(x));

sat_COG(sample_size+1,1:2) = initial_conditions(1,2); % take_out : forces Final P1P2 to be initial conditons. Would imply a perfect orbit, take out for real orbit

P_differ_angle = omega*(T/sample_size); 

e_hat_rP = [cos(theta_s); sin(theta_s)]; % So that P1(1) and P2(2) are intitial postions

while i3 < sample_size + 1 % take_out + 1 : forces Final P1P2 to be initial conditons. Would imply a perfect orbit, take out for real orbit
    i3 = i3 + 1;
    P1(i3, 1:2) = sat_COG(i3, 1:2) - (alt*e_hat_rP).';
    P2(i3, 1:2) = sat_COG(i3, 1:2) + (alt*e_hat_rP).';
    P_dist_hold(i3) = sqrt((P1(i3,1)-P2(i3,1))^2+(P1(i3,2)-P2(i3,2))^2);
    P1_dist(i3) = sqrt(P1(i3,1)^2+P1(i3,2)^2);
    P2_dist(i3) = sqrt(P2(i3,1)^2+P2(i3,2)^2);
    P_angle = P_angle + P_differ_angle;
    P_hold(i3) = P_angle;
    e_hat_rP = [cos(P_angle); sin(P_angle)];
end
P1X = P1(1:end,1); P1Y = P1(1:end,2);
P2X = P2(1:end,1); P2Y = P2(1:end,2);
toc
plot(P1X,P1Y)
plot(P2X,P2Y)
legend COG Earth P1 P2
hold off

figure(2)
plot(P_hold)
title('P angle')

figure(3)
plot(P_dist_hold)
title('P distance')
ylim([sat_length-1 sat_length+1])

figure(4)
plot(P1_dist)
title('P1 dist')

figure(5)
plot(P2_dist)
title('P2 dist')

function dU = loving_life(~,U)
    dU = zeros(size(U));
    x = U(1);
    y = U(2);
    u = U(3);
    v = U(4);
    % Equations for delta of each value
    r = sqrt(x^2 + y^2);
    
    dU(1) = u;
    dU(2) = v;
    dU(3) = -(mu*x)/r^3;
    dU(4) = -(mu*y)/r^3;
end

end