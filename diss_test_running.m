function diss_test_running
tic
fprintf('Minimum speed difference achived =')
sat_start = 10000;
sat = sat_start;
sat_max = 10000;
sat_1 = sat;
sat_step = 100;
rot_start = 0;
rot_step = 1;
rot_max = 4;
i = 1;
j = 0;
min_store = zeros(length(i));
max_store = zeros(length(i));


while sat <= sat_max
    j = j + 1;
    k = 1;
    rot = rot_start;
while rot <= rot_max
    [min,max,P1_velo,time] = diss_test_function_1(rot,sat);
    fprintf('%gkm/s\n',min)
    min_store(k,j) = min;
    max_store(k,j) = max;
    rot = rot + rot_step;
    k = k + 1;
    figure(2)
    plot(time,P1_velo)
    hold on
    xlim([0 inf])
    grid on
end

sat = sat + sat_step;

end

while i <= j
%    min_total = mink(min_store(1:end,i),2);
 %   fprintf('\nThe two smallest values for %gkm cable are:\n%gkm/s\n%gkm/s\n' ...
 %       ,sat_1,min_total)
    sat_1 = sat_1 + sat_step;
    i = i + 1;
end

hold off

figure(3)
bar3(min_store)
title('Minimum speeds at different lengths and rotations')
ylabel('Number of rotations per orbit')
zlabel('Minimum speed (km/s)')
xlabel('Length (\times100km)')
set(gca,'YTick',0:1:k)
set(gca,'YTickLabel',rot_start-rot_step:rot_step:rot_max)
set(gca,'XTick',0:1:j)
set(gca,'XTickLabel',sat_start-sat_step:sat_step:sat_max)

figure(4)
bar3(max_store)
title('Maximum speeds at different lengths and rotations')
ylabel('Number of rotations per orbit')
zlabel('Minimum speed (km/s)')
xlabel('Length (\times100km)')
set(gca,'YTick',0:1:k)
set(gca,'YTickLabel',rot_start-rot_step:rot_step:rot_max)
set(gca,'XTick',0:1:j)
set(gca,'XTickLabel',sat_start-sat_step:sat_step:sat_max)
toc
end