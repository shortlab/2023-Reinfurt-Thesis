clear
clc

%% Set up filename variables

chip_number='68123';                                            % Choose '48901' or '48902' or '49242'
experiment_number=[1:13];
max_temperature='600';                                         % Choose '1000' or '800' or '600' in C
heating_rate=[1000,500,300,100,50,30,10,5,3,1,0.5,0.3,0.1];


%% Read in data from text file.
for sampnum=1:10
%if sampnum~=22||14
%sampnum=5;

for rate=13                                                     % Choose the heating rate to import here
    
    %data=readmatrix(strcat(chip_number,'_',num2str(experiment_number(rate)),'_Blank_',max_temperature,'_',num2str(heating_rate(rate)),'Kps.txt'));
     data=readmatrix(strcat('1E12_',num2str(sampnum),'.txt'));
     %data=readmatrix('1E11_2_option.txt');
end


%% Split concatenated data into segments

segment_number=1;                                                   
outputRow=0;

for inputRow=1:size(data,1)
    
    outputRow=outputRow+1;
    
    % Split into each segment
   
    if data(inputRow,1) == 'NaN'
        outputRow=outputRow-1;
        fprintf("test")
        continue
    elseif data(inputRow,1) == 0
        segment_number=segment_number+1;
        outputRow=1;
    end
    
    % Allocate data into variables 
    
    time_all(outputRow,segment_number)=data(inputRow,2);           % s
    temp_sample_all(outputRow,segment_number)=data(inputRow,3);    % C
    temp_ref_all(outputRow,segment_number)=data(inputRow,4);       % C
    heatflow_all(outputRow,segment_number)=data(inputRow,5);       % mW
    
end


%% Separate into heating segments

%[rowEnd]=find(isnan(time_all(:,2)),1)-1;

%time_heat=cat(2,time_all(1:rowEnd,2),time_all(1:rowEnd,4),time_all(1:rowEnd,6), time_all(1:rowEnd,8));
%temp_sample_heat=cat(2,temp_sample_all(1:rowEnd,3),temp_sample_all(1:rowEnd,5),temp_sample_all(1:rowEnd,7), temp_sample_all(1:rowEnd,9));
%temp_ref_heat=cat(2,temp_ref_all(1:rowEnd,2),temp_ref_all(1:rowEnd,4),temp_ref_all(1:rowEnd,6), temp_ref_all(1:rowEnd,8));
%heatflow_heat=cat(2,heatflow_all(1:rowEnd,3),heatflow_all(1:rowEnd,5),heatflow_all(1:rowEnd,7), heatflow_all(1:rowEnd,9));

%% Separate into heating segments

[rowEnd]=find(isnan(time_all(:,2)),1)-1;

time_heat=time_all(1:rowEnd,2);
temp_sample_heat=temp_sample_all(1:rowEnd,2);
temp_ref_heat=temp_ref_all(1:rowEnd,2);
heatflow_heat=heatflow_all(1:rowEnd,2);

for ii=3:segment_number
time_heat=cat(2,time_heat,time_all(1:rowEnd,ii));
temp_sample_heat=cat(2,temp_sample_heat,temp_sample_all(1:rowEnd,ii));
temp_ref_heat=cat(2,temp_ref_heat,temp_ref_all(1:rowEnd,ii));
heatflow_heat=cat(2,heatflow_heat,heatflow_all(1:rowEnd,ii));
end




%% Plot raw data for each heating segment

for heat=1:segment_number-1
    
    plot(temp_sample_heat(:,heat),heatflow_heat(:,heat),'-','LineWidth',0.5,'DisplayName',num2str(heat))
    hold on
    
end

set(gca,'FontSize',15)
xlabel('Temperature (\circC)','FontSize',20)
ylabel('Power (mW)','FontSize',20)
title(strcat('Epoxy Sample #',num2str(sampnum),' (500 Kps, 100 Cycles)'))
%title(legend,'Heat','FontWeight','normal')
%title(strcat('Heating rate =',{' '},num2str(heating_rate(rate)),' K/s'),'FontSize',20,'FontWeight','normal')
%legend('Location','best')


%% Eval Tg with Sigmoid Fit

Tg=zeros(1,49);
T_index=1;
testsse=zeros(1,49);

% fit sigmoid function
tempiso=temp_sample_heat(601:2001,2);

for ii=1:segment_number/2-1

    if mod(ii,2)==0
    heatiso=heatflow_heat(601:2001,ii);
    
    x=tempiso;
    y=heatiso;
    %f=fit(temp(900:2800), smoothheat,'exp2')
    %plot(f,temp(900:2800), smoothheat)
    
    ft = fittype('L - (U-L)/(1 + exp(-4*log(3)*(x-xmid)/xscale80))','indep','x');
         
    [mdl,gof] = fit(x,y,ft,'start',[-.5 ,.5, 100,20]);
    Tg(T_index)=mdl.xmid;
    CI=confint(mdl);
    model_err(T_index)=(CI(2,3)-CI(1,3))/2;
    testsse(T_index)=gof.sse;
    T_index=T_index+1;
    end
    
end

figure(2)
plot(x,y)
grid on
hold on
plot(mdl)
xlabel('Temp [C]')
ylabel('Heatflow [mW]')
title('20th Heating Cycle (Sigmoid Fit)')
xline(mdl.xmid);
%xline(mdl.xmid + mdl.xscale80/2);
%xline(mdl.xmid - mdl.xscale80/2);

figure(7)
plot(1:length(testsse),testsse)
grid on
xlabel('Cycle Number')
ylabel('Goodness of Fit (SSE)')
title(strcat('Sigmoid GOF (Sample # ',num2str(sampnum),')'), 'FontSize',18)

%% residual

resid=mdl(x)-y;

%figure(5)
%plot(x,resid);
%xlabel('Temp [C]')
%ylabel('Heatflow [mW]')
%title('Residual of Sigmoid Fit (Cycle 20)')
%grid("on")



%% Stats
%avg_Tg=mean(Tg);
%stan_dev=zeros(length(Tg),1);
%num_samp=1:length(Tg);
%for ii=2:length(Tg)
    %stan_dev(ii)=std(Tg(1:ii))/sqrt(ii);
%end

model_w=1./(model_err.^2);
avg_Tg=sum(model_w.*Tg)/sum(model_w);
stan_dev=zeros(length(Tg),1);
num_samp=1:length(Tg);
for ii=2:length(Tg)
    stan_dev(ii)=1/sqrt(sum(model_w(1:ii)));
end


stan_dev(1)=stan_dev(2);

figure(3)
plot(num_samp,stan_dev)
xlabel('Number of Cycles','FontSize',15)
ylabel('\sigma (K)','FontSize',15)
title(strcat('Standard Error of Tg (Sample # ',num2str(sampnum),')'), 'FontSize',18)
grid on

%% cycle Variance
cycvar=abs(Tg-avg_Tg);

figure(11)
plot(num_samp,cycvar)
xlabel('Cycle Number','FontSize',15)
ylabel('variance (K)','FontSize',15)
title(strcat('Variance of Tg (Sample # ',num2str(sampnum),')'), 'FontSize',18)
grid on

%% Outlier Removal
%count=0;
%for ii=1:length(Tg)
    
 %   if abs(Tg(ii-count)-avg_Tg)>=5*stan_dev(end)
 %       Tg(ii-count)=[];
 %       count=count+1;
 %   end
%end
%num_data=1:length(Tg)
%plot(num_data, Tg)

%% Multi_Samp_stats
%Tg_mult=zeros(1,11);
%Tg_mult_std=zeros(1,11);

Tg_mult(sampnum)=avg_Tg;
Tg_mult_std(sampnum)=stan_dev(end);

mult_gof(sampnum)=mean(testsse);
mult_var(sampnum)=mean(cycvar);
end
%end

%% more multi Stats
avg_Tg_mult=mean(Tg_mult);
stan_dev_mult=zeros(length(Tg_mult),1);
num_samp_mult=1:length(Tg_mult);

for ii=2:length(Tg_mult)
    stan_dev_mult(ii)=(sqrt(sum(Tg_mult_std(1:ii).^2)))./sqrt(ii);
end

stan_dev_mult(1)=Tg_mult_std(1);
stan_dev_mult(end);

%% Outlier Removal Again
count=0;
for ii=1:length(Tg_mult)
    
    if abs(mult_gof(ii)-median(mult_gof))>=2.4*std(mult_gof)
        Tg_mult(ii-count)=[];
        Tg_mult_std(ii-count)=[];
        count=count+1;
    end
end


%% Re-evaluate stats
avg_Tg_mult=mean(Tg_mult);
stan_dev_mult=zeros(length(Tg_mult),1);
num_samp_mult=1:length(Tg_mult);

for ii=2:length(Tg_mult)
    stan_dev_mult(ii)=(sqrt(sum(Tg_mult_std(1:ii).^2)));%./sqrt(ii);
end

stan_dev_mult(1)=Tg_mult_std(1);

%% Plot error
figure(4)
errorbar(num_samp_mult,Tg_mult,Tg_mult_std,'o-')
xlabel('Number of Samples','FontSize',15)
ylabel('Tg (\circC)','FontSize',15)
title('Average Tg for Multiple Samples','FontSize',18)
grid on



figure(5)
plot(num_samp_mult,stan_dev_mult)
xlabel('Number of Samples','FontSize',15)
ylabel('\sigma (K)','FontSize',15)
title('Standard Deviation of Tg for Multiple Samples', 'FontSize',18)
grid on

%% 
figure(6)
plot(num_samp,Tg)
grid on
xlabel('Number of Cycles','FontSize',15)
ylabel('Tg','FontSize',15)
title(strcat('Tg over total Heating Cycles (Sample # ',num2str(sampnum),')'), 'FontSize',18)

%% 
tot_samp=1:10;
figure(8)
plot(tot_samp,mult_gof)
grid on
xlabel('Sample Number','FontSize',15)
ylabel('Tg','FontSize',15)
title('Sigmoid GOF (All Samples)', 'FontSize',18)

figure(9)
plot(tot_samp,mult_var)
grid on
xlabel('Sample Number','FontSize',15)
ylabel('Variance (K)','FontSize',15)
title('Variance of Tg (All Samples)', 'FontSize',18)

%% Final Plot
irrad=[0,1*10^11, 3.5*10^10];
Tg_rad=[116.3183, 82.0996, 100.05];
y_err=[1.187, 2.0782, 1.3301];

figure(12)
semilogx(irrad,Tg_rad)
errorbar(irrad,Tg_rad,y_err,'o')
grid on
xlabel('Irradiation Fluence (ions/cm^2)','FontSize',15)
ylabel('Tg (\circC)','FontSize',15)
title('Tg change due to Irradiation', 'FontSize',18)
xlim([-.1*10^11, .12*10^12])
ylim([75,120])

%% New STD plot
stan_dev_mult_new=zeros(1,length(Tg_mult));
for ii=2:length(Tg_mult)
    stan_dev_mult_new(ii)=(std(Tg_mult(1:ii))/sqrt(ii));
end

%% std plot
figure(13)
plot(num_samp_mult,stan_dev_mult_new)
xlabel('Number of Samples','FontSize',15)
ylabel('\sigma (K)','FontSize',15)
title('Standard Error of Tg for Multiple Samples', 'FontSize',18)
grid on


%%
writematrix(Tg_mult, 'Tg_1E12.dat')
writematrix(Tg_mult_std, 'Tg_std_1E12.dat')