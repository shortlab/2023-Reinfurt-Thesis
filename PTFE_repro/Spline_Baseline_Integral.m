function ints = Spline_Baseline_Integral(T,Q,x)
% For integrating Flash DSC PTFE cooling data (from run3). 
% Inputs: T,Q,x imported/loaded data in [C,mW,slice number]
% Outputs: the integrals for each sample based on spline baselines.
% Assumptions: Peaks are within [250,360] without hook contribution (don't 
% actually need to trim for spline however). The construction lines for the 
% splines are fit from [220,260] to [350,365]. 

%% Correctly orient data (assuming cooling)
Q = flipud(Q);
T = flipud(T);
n = length(x);

%% Spline baseline

ints = zeros(1,n);

% straight construction segments
T_ll = 75;
T_lr = 85;
T_rl = 135;
T_rr = 145;

for i = 1:n
    T_current = T(:,i);
    Q_current = Q(:,i);
    
    % Get correct indexes and data ranges
    left = find(T_current > T_ll & T_current < T_lr);
    right = find(T_current > T_rl & T_current < T_rr);
    corrected = find(T_current > T_ll & T_current < T_rr);
    T_left = T_current(left);
    T_right = T_current(right);
    Q_left = Q_current(left);
    Q_right = Q_current(right);
    T_corrected = T_current(corrected);
    Q_temp = Q_current(corrected);
    
    % Fit left and right straights to line
    T_fit = [T_right; T_left];
    fit_right = [ones(length(T_right),1) T_right]\Q_right;
    fit_left = [ones(length(T_left),1) T_left]\Q_left;
    Q_fit = [fit_right(1)+fit_right(2)*T_right; fit_left(1)+fit_left(2)*T_left];
    
    % Create spline
    Q_spline = spline(T_fit,Q_fit,T_corrected);
    Q_corr = Q_temp - Q_spline; 
    ints(i) = trapz(T_corrected,Q_corr);
       
end


%% Figure

figure;

hold on
plot(T_corrected,Q_temp,'b','LineWidth',0.5)
% plot(T_fit,Q_fit,'b.','LineWidth',1)
% plot(T_corrected,Q_spline,'b--','LineWidth',0.5)
plot(T_corrected,Q_temp,'c','LineWidth',0.5)
plot(T_fit,Q_fit,'c.','LineWidth',1)
plot(T_corrected,Q_spline,'c--','LineWidth',0.5)
% 
title('FDSC: Spline baseline defined for cooling data')
% % title('FDSC: Spline baseline defined for cooling data (set 4, slice 23)')
% % title('FDSC: Spline baseline defined for cooling data (set 4, slice 32)')
xlabel('Temperature [^{\circ}C]')
ylabel('Heat Flow [mW]')
% 
%  
%h = zeros(2, 1);
%h(1) = plot(NaN,NaN,'b-');
%h(2) = plot(NaN,NaN,'c-');
%legend(h, 'Set 4, slice 23','Set 4, slice 32')

end

