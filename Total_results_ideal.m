%Open diss_test_running_surf.mat
%A = load('diss_test_running_surf.mat')

figure(2)
surf(abs(min_store))
title('Absolute Minimum Speed Difference')
ylabel('Number of rotations per orbit')
zlabel('Minimum speed (km/s)')
xlabel('Length (km)')
set(gca,'YTick',1:10:91)
set(gca,'YTickLabel',{1:10}) %rot_start-rot_step:rot_step:rot_max)
set(gca,'XTick',0:5:50)
set(gca,'XTickLabel',{0:500:5000})%sat_start-sat_step:sat_step:sat_max)
min_v_c = colorbar;
min_v_c.Label.String = 'km/s';
colormap('turbo')
shading interp
view(90,-90)

clear
Re = 6371;
Karman_line = 100;
mu = 398600;

rot_min = 1;
rot_max = 10;
rot_step = 0.1;
rot = rot_min;

r_min = 100;
r_max = 5000;
r_step = 100;
r = r_min;

T = [];
i = 1;
j = 1;
k = 1;

for r = r_min:r_step:r_max
    T(i,j) = 2*pi*sqrt((Re+Karman_line+(r/2))^3/mu);
    for rot = rot_min:rot_step:rot_max
        RPM = rot/(T(i,j)/60);
        g_force(i,k) = (RPM/1000)^2*1.118*r*1000^2;
        k = k + 1;
    end
    i = i + 1;
    k = 1;
    if r == r_max
        i = 1;
        j = j + 1;
    end

end

figure(3)
surf(g_force)
title('Average g-force')
xlabel('Number of rotations per orbit'),ylabel('Length (km)'),zlabel('g force')
set(gca,'XTick',1:10:91)
set(gca,'XTickLabel',{1:10}) %rot_start-rot_step:rot_step:rot_max)
set(gca,'YTick',0:5:50)
set(gca,'YTickLabel',{0:500:5000})%sat_start-sat_step:sat_step:sat_max)
g_f_c = colorbar;
g_f_c.Label.String = 'g';
colormap('turbo')
view(0,90)
shading interp