function diss_test2_2
fprintf('\n')
warning('off')

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

sat_radius = sat_length/2; %km
a = RE + sat_radius + karmin_line; % Distance from centre of the Earth
theta_s = 0; %theta at the start of the run

% Orbit information
T = 2 * pi * sqrt(a^3/mu); % Time period
v_orb = (2*pi*a)/T; % Orbital velocity
delta_theta = 2*pi*rot_orbit; % Total angle change over orbit 

omega = delta_theta/T; % Angular velocity

% Cartesian Coordinates
% As the satellite is orbiting the equator, we will look in the [x y]
% coordinates
e_hat_r = [cos(theta_s); sin(theta_s)]; % Unit vector
r_ = a*e_hat_r; % Postions vector [X Y 0]
e_hat_theta = [-sin(theta_s); cos(theta_s)];
v_ = v_orb*e_hat_theta; % Velocity vector [x y 0]

r_P1_orb = a-sat_radius; % v = d/t where d = pi*diameter
r_P2_orb = a+sat_radius;

r_P1 = r_P1_orb*e_hat_r;
r_P2 = r_P2_orb*e_hat_r;

theta = theta_s;

% Initial conditions
initial_conditions = [r_(1) r_(2) v_(1) v_(2) theta 0 r_P1(1) r_P1(2) r_P2(1) r_P2(2)];
tic
fprintf('Please wait while I run some calculations,\nyour buisness is important to us... *♫Hold Music♫*\n\n')

U = initial_conditions;
options = odeset('abstol',0.001,'reltol',0.000001);

[~,U] = ode45(@loving_life, [0 T], U, options);
x = U(:,1); y = U(:,2); u_CoG = U(:,3); v_CoG = U(:,4); theta = U(:,5); time = U(:,6); x_P1 = U(:,7); y_P1 = U(:,8); x_P2 = U(:,9); y_P2 = U(:,10);
 
fprintf('Thank you for waiting \n\n')

%fprintf('The difference between the inital and final position of the trajectory is %0.2gkm (%0.4gkm)\n',close_diff,close_diff)

figure(1)
plot(x,y,'color',[61/255 217/255 201/255])

xlabel('X (km)'),ylabel('Y (km)'),zlabel('Z (km)'),title('The trajectory of the spacecraft around Earth with recalculated error')
grid on
hold on
[X,Y,Z] = sphere;
surf(X*RE,Y*RE,Z*RE)
axis equal

toc

plot(x_P1, y_P1)

plot(x_P2, y_P2)
legend CoG Earth P1 P2
hold off

P_dist_hold = sqrt((x_P1-x_P2).^2+(y_P1-y_P2).^2);
P1_dist = sqrt(x_P1.^2+y_P1.^2);
P2_dist = sqrt(x_P2.^2+y_P2.^2);


figure(2)
plot(time,P_dist_hold)
title('P distance')

figure(3)
plot(time,P1_dist)
title('P1 distance')

figure(4)
plot(time,P2_dist)
title('P2 distance')


fprintf('%g\n',mean(P_dist_hold))
fprintf('%g\n',mean(P1_dist))
fprintf('%g\n',mean(P2_dist))

function dU = loving_life(t,U)
    dU = zeros(size(U));
    x = U(1);
    y = U(2);
    u_CoG = U(3);
    v_CoG = U(4);
    theta = U(5);
    time = U(6);
    x_P1 = U(7); 
    y_P1 = U(8); 
    x_P2 = U(9); 
    y_P2 = U(10); 
    
    
    % Equations for delta of each value
    r = sqrt(x^2 + y^2);
    
    dU(1) = u_CoG;
    dU(2) = v_CoG;
    dU(3) = -(mu*x)/r^3;
    dU(4) = -(mu*y)/r^3;

    dU(5) = omega*(t-time);
    dU(6) = t - time;
    theta = omega*t;
    e_hat_rP = [cos(theta); sin(theta)];

    dU(7) = x - (sat_radius*e_hat_rP(1)) - x_P1;
    dU(8) = y - (sat_radius*e_hat_rP(2)) - y_P1;
    dU(9) = x + (sat_radius*e_hat_rP(1)) - x_P2;
    dU(10) = y + (sat_radius*e_hat_rP(2)) - y_P2;
    
end
end