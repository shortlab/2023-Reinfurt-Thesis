%%
[temp, heat2, time, type2]=importFSC_var("54280_Slice3_cool3.txt");
[temp, heat4, time, type4]=importFSC_var("54280_Slice3_heat4.txt");
[temp, heat6, time, type6]=importFSC_var("54280_Slice4_cool3.txt");
[temp, heat8, time, type8]=importFSC_var("54280_Slice4_heat4.txt");
[temp, heat10, time, type10]=importFSC_var("54280_Slice5_cool3.txt");
[temp, heat12, time, type12]=importFSC_var("54280_Slice5_heat4.txt");
[temp, heat14, time, type14]=importFSC_var("54280_Slice6_cool3.txt");
[temp, heat16, time, type16]=importFSC_var("54280_Slice6_heat4.txt");


%%
figure(1)
hold on
plot(temp, heat2)
plot(temp, heat4)
%plot(temp, heat6)
%%plot(temp, heat8)
%plot(temp, heat10)
%plot(temp, heat12)
%plot(temp, heat14)
%plot(temp, heat16)
hold off
xlabel('Temp [C]')
ylabel('Heatflow [mW]')
title('PTFE Data (Connick)')
grid("on")


