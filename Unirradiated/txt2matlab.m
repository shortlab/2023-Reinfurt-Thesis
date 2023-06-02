clear
clc

%% Set up filename variables

chip_number='68123';                                            % Choose '48901' or '48902' or '49242'
experiment_number=[1:13];
max_temperature='600';                                         % Choose '1000' or '800' or '600' in C
heating_rate=[1000,500,300,100,50,30,10,5,3,1,0.5,0.3,0.1];


%% Read in data from text file.
%for sampnum=1:30
%if sampnum~=22||14
sampnum=9;

for rate=13                                                     % Choose the heating rate to import here
    
    %data=readmatrix(strcat(chip_number,'_',num2str(experiment_number(rate)),'_Blank_',max_temperature,'_',num2str(heating_rate(rate)),'Kps.txt'));
     data=readmatrix(strcat('sig_red_100_',num2str(sampnum),'.txt'));
     %data=readmatrix('sigma_red_last_100.txt');
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

figure(1)
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
tempiso=temp_ref_heat(1201:2501,2);

for ii=1:segment_number/2-1

    if mod(ii,2)==0
    heatiso=heatflow_heat(1201:2501,ii);
    
    x=tempiso;
    y=heatiso;
    %f=fit(temp(900:2800), smoothheat,'exp2')
    %plot(f,temp(900:2800), smoothheat)
    
    ft = fittype('L - (U-L)/(1 + exp(-4*log(3)*(x-xmid)/xscale80))','indep','x');
         
    [mdl,gof] = fit(x,y,ft,'start',[-.5 ,.5, 130,20]);
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

figure(3)
plot(1:length(testsse),testsse)
grid on
xlabel('Cycle Number')
ylabel('Goodness of Fit (SSE)')
title(strcat('Sigmoid GOF (Sample # ',num2str(sampnum),')'), 'FontSize',18)

%% residual

resid=mdl(x)-y;

figure(4)
plot(x,resid);
xlabel('Temp [C]')
ylabel('Heatflow [mW]')
title('Residual of Sigmoid Fit (Cycle 20)')
grid("on")

%%
h=kstest(resid)
A=fft(resid);

%plot(1:length(A),A)
%%
time_mdl=time_heat(1201:2501,2);
time_tot=time_mdl(end)-time_mdl(1);
Fs=length(resid)/time_tot;
period=1/Fs;
%f = Fs*(0:1/time_tot:(time_tot/2))/time_tot;
f=Fs.*(0:2/length(time_mdl):1-1/length(time_mdl));
power=A.*conj(A);
C0=fftshift(A);
P1=(A/time_tot);
[F,f1]=easyFFT(resid,length(resid),1,Fs);
%%
figure(25)
plot(f1,F.*conj(F))
xlim([0,1000])
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('FFT of Residual of Sigmoid Fit')

%% Outlier Removal
count=0;
for ii=1:length(Tg)
    
    if abs(testsse(ii)-median(testsse))>=5*std(testsse)
        Tg(ii-count)=[];
        model_err(ii-count)=[];
        count=count+1;
    end
end
num_data=1:length(Tg);
%figure(20);
%plot(num_data, Tg)
%grid on
%xlabel('Number of Cycles','FontSize',15)
%ylabel('Tg','FontSize',15)
%title(strcat('Tg over total Heating Cycles (Sample # ',num2str(sampnum),')'), 'FontSize',18)


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

figure(5)
plot(num_samp,stan_dev)
xlabel('Number of Cycles','FontSize',15)
ylabel('\sigma (K)','FontSize',15)
title(strcat('Standard Deviation of Tg (Sample # ',num2str(sampnum),')'), 'FontSize',18)
grid on

%% cycle Variance
cycvar=abs(Tg-avg_Tg);

figure(6)
plot(num_samp,cycvar)
xlabel('Cycle Number','FontSize',15)
ylabel('variance (K)','FontSize',15)
title(strcat('Variance of Tg (Sample # ',num2str(sampnum),')'), 'FontSize',18)
grid on

%% Outlier Removal
%count=0;
%for ii=1:length(Tg)
    
    %if abs(testsse(ii)-median(testsse))>=5*std(testsse)
        %Tg(ii-count)=[];
        %count=count+1;
    %end
%end
%num_data=1:length(Tg)
%figure(20);
%plot(num_data, Tg)
%grid on
%xlabel('Number of Cycles','FontSize',15)
%ylabel('Tg','FontSize',15)
%title(strcat('Tg over total Heating Cycles (Sample # ',num2str(sampnum),')'), 'FontSize',18)


%% Multi_Samp_stats
%Tg_mult=zeros(1,11);
%Tg_mult_std=zeros(1,11);

Tg_mult(sampnum)=avg_Tg;
Tg_mult_std(sampnum)=stan_dev(end);

mult_gof(sampnum)=mean(testsse);
mult_var(sampnum)=mean(cycvar);
%end
%end

%% Outlier Removal Again
count=0;
for ii=1:length(Tg_mult)
    
    if abs(mult_gof(ii)-median(mult_gof))>=2.4*std(mult_gof)
        Tg_mult(ii-count)=[];
        Tg_mult_std(ii-count)=[];
        count=count+1;
    end
end

%% more multi Stats
avg_Tg_mult=mean(Tg_mult)
stan_dev_mult=zeros(length(Tg_mult),1);
num_samp_mult=1:length(Tg_mult);

%for ii=2:length(Tg_mult)
    %stan_dev_mult(ii)=(sqrt(sum(Tg_mult_std(1:ii).^2)))./sqrt(ii);
%end

for ii=2:length(Tg_mult)
    stan_dev_mult(ii)=(std(Tg_mult(1:ii)/sqrt(ii)));
end

stan_dev_mult(1)=Tg_mult_std(1);
stan_dev_mult(end)

%% Outlier Removal Again
%count=0;
%for ii=1:length(Tg_mult)
    
    %if abs(mult_gof(ii)-median(mult_gof))>=5*std(mult_gof)
        %Tg_mult(ii-count)=[];
        %Tg_mult_std(ii-count)=[];
        %count=count+1;
    %end
%end


%% Re-evaluate stats
%stan_dev_mult=zeros(length(Tg_mult),1);
%num_samp_mult=1:length(Tg_mult);

%for ii=2:length(Tg_mult)
    %stan_dev_mult(ii)=(sqrt(sum(Tg_mult_std(1:ii).^2)))./sqrt(ii);
%end

%stan_dev_mult(1)=Tg_mult_std(1);

%% Plot error
figure(7)
errorbar(num_samp_mult,Tg_mult,Tg_mult_std,'o-')
xlabel('Number of Samples','FontSize',15)
ylabel('Tg (K)','FontSize',15)
title('Average Tg for Multiple Samples','FontSize',18)
grid on



figure(8)
plot(num_samp_mult,stan_dev_mult)
xlabel('Number of Samples','FontSize',15)
ylabel('\sigma (K)','FontSize',15)
title('Standard Deviation of Tg for Multiple Samples', 'FontSize',18)
grid on

%% 
figure(9)
plot(num_samp,Tg)
grid on
xlabel('Number of Cycles','FontSize',15)
ylabel('Tg','FontSize',15)
title(strcat('Tg over total Heating Cycles (Sample # ',num2str(sampnum),')'), 'FontSize',18)

%% 
figure(10)
plot(1:length(mult_gof),mult_gof)
grid on
xlabel('Sample Number','FontSize',15)
ylabel('Goodness of Fit (SSE)','FontSize',15)
title('Sigmoid GOF (All Samples)', 'FontSize',18)

figure(11)
plot(1:length(mult_var),mult_var)
grid on
xlabel('Sample Number','FontSize',15)
ylabel('Variance (K)','FontSize',15)
title('Variance of Tg (All Samples)', 'FontSize',18)

%%
writematrix(Tg_mult, 'Tg_0.dat')
writematrix(Tg_mult_std, 'Tg_std_0.dat')
