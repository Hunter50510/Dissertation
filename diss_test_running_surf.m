%function diss_test_running_surf
tic
fprintf('Minimum speed difference achived =')
sat_start = 100;
sat = sat_start;
sat_max = 5000;
sat_1 = sat;
sat_step = 100;
rot_start = 1;
rot_step = 0.1;
rot_max = 10;
i = 1;
j = 0;
min_store = zeros(length(i));
max_store = zeros(length(i));
while sat <= sat_max
    j = j + 1;
    k = 1;
    rot = rot_start;
while rot <= rot_max
    [min,max] = diss_test_function(rot,sat);
    fprintf('%gkm/s\n',min)
    min_store(k,j) = min;
    max_store(k,j) = max;
    rot = rot + rot_step;
    k = k + 1;
end

sat = sat + sat_step;

end

while i <= j
    min_total = mink(min_store(1:end,i),2);
    fprintf('\nThe two smallest values for %gkm cable are:\n%gkm/s\n%gkm/s\n' ...
        ,sat_1,min_total)
    sat_1 = sat_1 + sat_step;
    i = i + 1;
end

figure(2)
surf(min_store)
title('Minimum speeds at different lengths and rotations')
ylabel('Number of rotations per orbit')
zlabel('Minimum speed (km/s)')
xlabel('Length (km)')
set(gca,'YTick',0:1:k)
set(gca,'YTickLabel',rot_start-rot_step:rot_step:rot_max)
set(gca,'XTick',0:1:j)
set(gca,'XTickLabel',sat_start-sat_step:sat_step:sat_max)

figure(3)
surf(max_store)
title('Maximum speeds at different lengths and rotations')
ylabel('Number of rotations per orbit')
zlabel('Minimum speed (km/s)')
xlabel('Length (km)')
set(gca,'YTick',0:1:k)
set(gca,'YTickLabel',rot_start-rot_step:rot_step:rot_max)
set(gca,'XTick',0:1:j)
set(gca,'XTickLabel',sat_start-sat_step:sat_step:sat_max)

toc
%end