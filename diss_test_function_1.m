function [min_delta_v,max_delta_v,diff_velo_P1_Earth_removed,time_removed] = diss_test_function_1(rot,sat_l)
fprintf('\n')
warning('off')

% Provides orbit of satellite in Circular orbit, i = 0

% Earth
mu = 398600; %km^3 s^-2
RE = 6371; %km
karmin_line = 100; %km

%Earth's rotation
T_Earth = 86164.0905; %One sidereal day
omega_Earth = 2*pi/T_Earth;
velo_Earth = omega_Earth*RE;

% input
rot_orbit = rot; %Number of rotations per orbit
sat_length = sat_l; %Length of tether
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

v_P1_orb = (2*pi*(a-sat_radius))/T;
v_P2_orb = (2*pi*(a+sat_radius))/T;

v_P1 = v_P1_orb.*e_hat_theta;
v_P2 = v_P2_orb.*e_hat_theta;

theta = theta_s;

% Initial conditions
initial_conditions = [r_(1) r_(2) v_(1) v_(2) theta 0 r_P1(1) r_P1(2)...
    r_P2(1) r_P2(2) v_P1(1) v_P1(2) v_P2(1) v_P2(2)];

U = initial_conditions;
options = odeset('abstol',0.0001,'reltol',0.00001);

[~,U] = ode23(@spinny_boi, [0 3*T], U, options);
x = U(:,1); y = U(:,2); u_CoG = U(:,3); v_CoG = U(:,4); theta = U(:,5); ...
    time = U(:,6); x_P1 = U(:,7); y_P1 = U(:,8); x_P2 = U(:,9); ...
    y_P2 = U(:,10); u_P1 = U(:,11); v_P1 = U(:,12); u_P2 = U(:,13);...
    v_P2 = U(:,14);
 
figure(1)
plot(x,y,'color',[61/255 217/255 201/255])

xlabel('X (km)'),ylabel('Y (km)'),zlabel('Z (km)'),title('The trajectory of the spacecraft around Earth with recalculated error')
grid on
hold on
[X,Y,Z] = sphere;
surf(X*RE,Y*RE,Z*RE)
axis equal

plot(x_P1, y_P1)

plot(x_P2, y_P2)
legend CoG Earth P1 P2
hold off

velo_P1 = sqrt(u_P1.^2+v_P1.^2);
diff_velo_P1_Earth = velo_P1 - velo_Earth;

data_length = length(diff_velo_P1_Earth);
removal = data_length/3;
diff_velo_P1_Earth_removed = diff_velo_P1_Earth(ceil(removal):ceil(data_length-removal));
time_removed = time(ceil(removal):ceil(data_length-removal))-T;

min_delta_v = min(diff_velo_P1_Earth_removed);
max_delta_v = max(diff_velo_P1_Earth_removed);

    function dU = spinny_boi(t,U)
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
    u_P1 = U(11); 
    v_P1 = U(12); 
    u_P2 = U(13); 
    v_P2 = U(14); 
    
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
    
    dU(7) = (x - x_P1) - (sat_radius*e_hat_rP(1));
    dU(8) = (y - y_P1) - (sat_radius*e_hat_rP(2));
    dU(9) = (x- x_P2) + (sat_radius*e_hat_rP(1));
    dU(10) = (y - y_P2) + (sat_radius*e_hat_rP(2));

    r_P1 = sqrt(x_P1^2 + y_P1^2);
    r_P2 = sqrt(x_P2^2 + y_P2^2);

    dU(11) = dU(7) - u_P1;
    dU(12) = dU(8) - v_P1; 
    dU(13) = dU(9) - u_P2;
    dU(14) = dU(10) - v_P2;
end
end