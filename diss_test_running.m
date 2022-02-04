clear
clc

tic
fprintf('\nMinimum speed difference achived =\n')
sat = 1500;


rot = 1;
min_store = zeros(length(rot));
while rot <= 10
min = diss_test_function(rot,sat);
fprintf('%gkm/s\n',min)
min_store(rot) = min;
rot = rot+1;

end

bar(min_store)
title('Min speed difference at',sat)
xlabel('Number of rotations')
ylabel('Minimum speed difference (km/s)')
grid on
toc