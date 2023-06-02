%%
[temp, heat2, time, type2]=importFSC_var("PCTFE_2.txt");
[temp, heat4, time, type4]=importFSC_var("PCTFE_4.txt");
[temp, heat6, time, type6]=importFSC_var("PCTFE_6.txt");
[temp, heat8, time, type8]=importFSC_var("PCTFE_8.txt");
[temp, heat10, time, type10]=importFSC_var("PCTFE_10.txt");
[temp, heat12, time, type12]=importFSC_var("PCTFE_12.txt");
[temp, heat14, time, type14]=importFSC_var("PCTFE_14.txt");
[temp, heat16, time, type16]=importFSC_var("PCTFE_16.txt");


%%
figure(1)
hold on
plot(temp, heat2)
plot(temp, heat4)
plot(temp, heat6)
plot(temp, heat8)
plot(temp, heat10)
plot(temp, heat12)
plot(temp, heat14)
plot(temp, heat16)
hold off
xlabel('Temp [C]')
ylabel('Heatflow [mW]')
title('PCTFE FSC Data')
grid("on")

%%
%enthvalheat = Spline_Baseline_Integral(temp,heat4,1)
