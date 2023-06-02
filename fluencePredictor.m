%% Quick calculator for flux from planar source
% Just need to change parameters in the first section and then run the
% script. (takes a few seconds to run)
% After the script finishes, the average flux over the sample surface and
% the time it takes to get to 10^11 alphas/cm^2 is stored in 
% flux_avg and time_days

%% Important dimenions

% Distance from source. If right up against the grate, it's ~0.45 cm.
d = 0.45;   % use cm.

% Estimated dimensions of the face of the sample. (width and height)
% Don't use anything less than 0.1 cm resolution.
w = 0.5;    % use cm. Dimension perpendicular to the grating. 
h = 1.5;    % use cm. Dimension parallel to the grating.

% Current activity of the source.
A = 417.98995;    % use uCi. 

%% Notes on assumptions
% Assumes center of sample is of center of source.
% Energy loss in air calculated using a stopping table from SRIM with an 
% approximation of air. 
% For comparison with Julie's GEANT simulation, her first sim came out with
% fraction of 0.0155923 hitting the surface of the sample, whereas I
% calculate a fraction of 0.015, which is a little low, but pretty close. I
% also have a slightly different energy distribution, but not totally off.

%% Stopping power: particle energy as a function of distance 

% Load stopping table

%load('StoppingTable.mat')

% Variable sTable has columns 
% Energy | Electronic stopping power | Nuclear Stopping Power
% keV    | MeV / (mg/cm2)            | MeV / (mg/cm2) 
% Multiply stopping power by 1.2300E-01 to convert keV / micron  

%dEdx_table = (sTable(:,2)+sTable(:,3))*1.2300E-01;   % total in keV/micron
%E_table = sTable(:,1);                               % energy in keV
%E_fine = linspace(min(E_table),max(E_table),600);
%dEdx_fine = interp1(E_table, dEdx_table,E_fine);
%ind = find(E_fine == 5300);
%Eback = E_fine(1:ind);
%dEdx = dEdx_fine(1:ind);

% Calculate the energy E of the 5.3 MeV alpha particles after traveling a
% distance D in air.

%D = zeros(size(Eback));
%for i = 2:length(Eback)
 %   D(i) = trapz(Eback(ind-i+1:end),1./dEdx(ind-i+1:end))*1e-4;
%end
%E = fliplr(Eback);

% Now, given a distance traveled in air, we can estimate the new energy of 
% that particle. D is in cm.

%% The source and sample setup

% Dimensions of source are 2x1 cm^2
% Splitting up the source into array of point sources, 0.01x0.01cm^2 each:
% Coordinates of each:
x_source = linspace(-1+0.005,1-0.005,200);
y_source = linspace(-0.5+0.005,0.5-0.005,100);

% Activity in uCi is converted to decays per second:
act = A*37000;
% Each of the 20,000 point sources has an activity of 
dA = act/20000;

% Doing the same for the sample:
x_sample = linspace(-w/2+0.005,w/2-0.005,w*100);
y_sample = linspace(-h/2+0.005,h/2-0.005,h*100);

%% Now to calculate the flux at each element of the sample surface
% Flux is in alphas/s/cm^2

% Each element is a distance r from each source at (x,y).
% The flux from each source is dA/(4*pi*r^2), which can then be summed up.
% And then we repeat for each element of the sample surface.

% Initialize a matrix to contain all the total fluxes at each element:
flux = zeros(length(x_sample),length(y_sample));
% Can also save all the energies at the sample surface, but it takes alot of time.
% Energy = zeros(1, length(x_sample)*length(y_sample)*length(x_source)*length(y_source));
% k = 1;
for i = 1:length(x_sample)
   for j = 1:length(y_sample)
       % Initialize a variable to hold the sum of fluxes.
       sum = 0;
       for a = 1:length(x_source)
           for b = 1:length(y_source)
               r = sqrt((x_source(a)-x_sample(i))^2+(y_source(b)-y_sample(j))^2+(0-d)^2)^0.5;
               f = dA/(4*pi*r^2);
               sum = sum+f;
%                Energy(k) = interp1(D,E,r);
%                k = k+1;
           end
       end
       flux(i,j) = sum;
   end
end

%% Average flux over the surface and the days to get to 10^11 alpha/cm^2

flux_avg = mean(mean(flux));
time_days = 3.5e11/flux_avg/3600/24;


