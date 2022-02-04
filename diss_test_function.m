function [min_delta_v] = diss_test_function(rot,sat_l)
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
%inputs = {'Number of rotations per orbit: (Recommended 1-3)','Length of tether: (Recommended 100-1800)'};
%title_input = 'Rotations and Length';
%dims = [1 40];
%hint = {'2','1500'};
%control = inputdlg(inputs,title_input,dims,hint);
rot_orbit = rot; %str2double(control{1}); %Number of rotations per orbit
sat_length = sat_l; %str2double(control{2}); %Length of tether
fprintf('Number of rotations = %g\n',rot_orbit)
fprintf('Length of tether = %gkm\n\n',sat_length)
%tic
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

v_P1 = v_P1_orb.*e_hat_theta;%((2*pi*(a-sat_radius))/T)*e_hat_theta;
v_P2 = v_P2_orb.*e_hat_theta;%((2*pi*(a+sat_radius))/T)*e_hat_theta;

theta = theta_s;

% Initial conditions
initial_conditions = [r_(1) r_(2) v_(1) v_(2) theta 0 r_P1(1) r_P1(2)...
    r_P2(1) r_P2(2) v_P1(1) v_P1(2) v_P2(1) v_P2(2)];

%fprintf('Please wait while I run some calculations,\nyour buisness is important to us... *♫Hold Music♫*\n\n')

U = initial_conditions;
options = odeset('abstol',0.0001,'reltol',0.00001);

[~,U] = ode23(@spinny_boi, [0 3*T], U, options);
x = U(:,1); y = U(:,2); u_CoG = U(:,3); v_CoG = U(:,4); theta = U(:,5); ...
    time = U(:,6); x_P1 = U(:,7); y_P1 = U(:,8); x_P2 = U(:,9); ...
    y_P2 = U(:,10); u_P1 = U(:,11); v_P1 = U(:,12); u_P2 = U(:,13);...
    v_P2 = U(:,14);
 
%fprintf('Thank you for waiting \n\n')

%fprintf('The difference between the inital and final position of the trajectory is %0.2gkm (%0.4gkm)\n',close_diff,close_diff)

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

P_dist_hold = sqrt((x_P1-x_P2).^2+(y_P1-y_P2).^2);
P1_dist = sqrt(x_P1.^2+y_P1.^2);
P2_dist = sqrt(x_P2.^2+y_P2.^2);


%figure(2)
%plot(time,P_dist_hold)
%title('P distance')
%ylim([min(P_dist_hold)-1 max(P_dist_hold)+1])

%figure(3)
%plot(time,P1_dist)
%title('P1 distance')
%ylim([min(P1_dist)-1 max(P1_dist)+1])

%figure(4)
%plot(time,P2_dist)
%title('P2 distance')
%ylim([min(P2_dist)-1 max(P2_dist)+1])

velo_P1 = sqrt(u_P1.^2+v_P1.^2);
%velo_P2 = sqrt(u_P2.^2+v_P2.^2);
%velo_CoG = sqrt(u_CoG.^2+v_CoG.^2);

%figure(5)
%plot(time,velo_P1)
%title('P1 velocity')
%ylim([min(velo_P1)-1 max(velo_P1)+1])

%figure(6)
%plot(time,velo_P2)
%title('P2 velocity')
%ylim([min(velo_P2)-1 max(velo_P2)+1])

%figure(7)
%plot(time,velo_CoG)
%title('CoG velocity')
%ylim([min(velo_CoG)-1 max(velo_CoG)+1])

diff_velo_P1_Earth = velo_P1 - velo_Earth;

%figure(8)
%plot(time, diff_velo_P1_Earth)
%hold on
%title('Difference between P1 and Earth''s velocity')
%xlabel('Time_(_s_)'), ylabel('Velocity (km/s)')
%grid on
%plot([T,T], [min(diff_velo_P1_Earth)-0.1, max(diff_velo_P1_Earth)+0.1])
%plot([2*T,2*T], [min(diff_velo_P1_Earth)-0.1, max(diff_velo_P1_Earth)+0.1])
%legend({'Speed difference', 'End of first orbit','End of second orbit'}, 'Location','southeast')
%hold off

data_length = length(diff_velo_P1_Earth);
removal = data_length/10;
diff_velo_P1_Earth_removed = diff_velo_P1_Earth(ceil(removal):ceil(data_length-removal));
%time_removed = time(ceil(removal):ceil(data_length-removal));

%figure(9)
%plot(time_removed,diff_velo_P1_Earth_removed)
%hold on
%title('Difference between P1 and Earth''s velocity')
%xlabel('Time_(_s_)'), ylabel('Velocity (km/s)')
%grid on
%grid on
%plot([T,T], [min(diff_velo_P1_Earth_removed)-0.1, ...
%    max(diff_velo_P1_Earth_removed)+0.1])
%plot([2*T,2*T], [min(diff_velo_P1_Earth_removed)-0.1, ...
%    max(diff_velo_P1_Earth_removed)+0.1])
%legend({'Speed difference', 'End of first orbit','End of second orbit'}, ...
%    'Location','southeast')
%hold off

min_delta_v = min(diff_velo_P1_Earth_removed);
%fprintf('The minimum difference between the Earth and P1 = %gkm/s\n\n',min_delta_v)

%toc

%fprintf('%g\n',mean(P_dist_hold))
%fprintf('%g\n',mean(P1_dist))
%fprintf('%g\n',mean(P2_dist))

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

    %e_hat_thetaP = [-sin(theta); cos(theta)];

    dU(11) = dU(7) - u_P1;
    dU(12) = dU(8) - v_P1; 
    %((2*pi*(y_P1-sat_radius))/T)*e_hat_thetaP(2); %(y - y_P1)*e_hat_thetaP(2);
    dU(13) = dU(9) - u_P2;
    dU(14) = dU(10) - v_P2;
end
end