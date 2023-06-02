[temp, heat2, time, type2]=importFSC_var("init_heat.txt");
[temp, heat4, time, type2]=importFSC_var("init_cool.txt");

figure(1)
hold on
plot(temp, heat2)
plot(temp, heat4)
hold off
xlabel('Temp [C]')
ylabel('Heatflow [mW]')
title('Unirradiated Epoxy (1st cycle)')
grid("on")