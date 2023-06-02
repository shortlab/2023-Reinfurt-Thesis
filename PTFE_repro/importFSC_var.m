function [T,Q,t,type] = importFSC_var(filename)
% Description: imports FSC data according to length of data.
% Inputs:  filename - string; file path to open.
% Outputs: T - column vector (double); Reference temperature [C].
%          Q - column vector (double); Heat flow [mW], exo up.
%          t - column vector (double); Time [s].
%          type - string; Checks if data is iso, heating, or cooling.
% Note: This function orients cooling data with T increasing.

%% Extract the data.
fid = fopen(filename,'r');
data = textscan(fid,'%f%f%f%f%f%[^\n\r]','HeaderLines',2);
fclose(fid);

%% Create output variable
T = data{5};
Q = data{3};
t = data{2};

%% Check iso, heating, or cooling.
% If cooling, will flip all arrays to order T as increasing

if ~any(T-T(1))         % For iso in FSC, all Tr should be identical.
    type = 'iso';       
elseif T(end)-T(1) > 0  % During heating, temperature increases.
    type = 'heat';
elseif T(end)-T(1) < 0  % During cooling, temperature decreases.
    type = 'cool';
    T = flipud(T);
    Q = flipud(Q);
    t = flipud(t);
end

end


