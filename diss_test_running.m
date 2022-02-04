function diss_test_running
tic
fprintf('Minimum speed difference achived =')
sat = 1000;
i = 1;
j = 0;
min_store = zeros(length(i));
max_store = zeros(length(i));
while sat <= 2000
    j = j + 1;
    i = i + 1;
    k = 1;
    rot = 1;
while k <= 5
    [min,max] = diss_test_function(rot,sat);
    fprintf('%gkm/s\n',min)
    min_store(k,j) = min;
    max_store(k,j) = max;
    rot = rot + 1;
    k = k + 1;
end

sat = sat + 1000;

end

figure(2)
bar3(min_store)
title('Minimum speeds at different lengths and rotations')
ylabel('Number of rotations per orbit')
zlabel('Minimum speed (km/s)')
xlabel('Length (\times1000km)')

figure(3)
bar3(max_store)
title('Maximum speeds at different lengths and rotations')
ylabel('Number of rotations per orbit')
zlabel('Minimum speed (km/s)')
xlabel('Length (\times1000km)')

toc
end