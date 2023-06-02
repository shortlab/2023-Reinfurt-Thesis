%% Script for loading hook recrystallization peak data
% I will load all five samples such that the (T,Q) data all has the same 
% length as decided by the shortest of the shortest measurements as measured
% from the hook of the cooling runs. This means I will need to check
% segment 28, or the 7th cooling segment's length for all the d2d's.

% For ease of enthalpy analysis, I will group each sample's peak area
% segments together so that the first 8 columns are the 8 d2d segments, and the
% ninth is the u2d cooling segment with the peak area (segment 4).

%% Find shortest segments:

% Note: the naming convention for "sample 5" (actually the first sample to
% be measured) is different. Also, it's stored in a different folder.

path_main = 'C:\Users\Rachel\Documents\Rachel\MIT\MIT Research\Short Lab\';
path_specific_1234 = 'Hook method\Exported data\Repeated data\';
path_specific_5 = 'Hook method\Exported data\Hook at Temp\';

% Example of file name for samples 1-4: 3_downto_down_28.txt
% Example of file name for samples 5: downto_down_28.txt
s3 = 'downto_down_28.txt';

% Just manually pull them in:
% Sample 1:
s1 = '1_';
filename = [path_main path_specific_1234 [s1 s3]];
[TA,~,~,~] = importFSC_var(filename);   % Automatically flips orders cooling data correctly (according to increasing temperature)

% Sample 2:
s1 = '2_';
filename = [path_main path_specific_1234 [s1 s3]];
[TB,~,~,~] = importFSC_var(filename);   

% Sample 3:
s1 = '3_';
filename = [path_main path_specific_1234 [s1 s3]];
[TC,~,~,~] = importFSC_var(filename);   

% Sample 4:
s1 = '4_';
filename = [path_main path_specific_1234 [s1 s3]];
[TD,~,~,~] = importFSC_var(filename);   

% Sample 5:
filename = [path_main path_specific_5 s3];
[TE,~,~,~] = importFSC_var(filename);   

% Check which of the vectors is the shortest:
L = min([length(TA) length(TB) length(TC) length(TD) length(TE)]);

%% Initialize arrays

T_PA_1 = zeros(L,9);
T_PA_2 = zeros(L,9);
T_PA_3 = zeros(L,9);
T_PA_4 = zeros(L,9);
T_PA_5 = zeros(L,9);

Q_PA_1 = zeros(L,9);
Q_PA_2 = zeros(L,9);
Q_PA_3 = zeros(L,9);
Q_PA_4 = zeros(L,9);
Q_PA_5 = zeros(L,9);

%% Set up for loops to load the rest of the data
% Going to be saving the indexes (end-L+1):end for each data set

% Not changing:
path_main = 'C:\Users\Rachel\Documents\Rachel\MIT\MIT Research\Short Lab\';
path_specific_1234 = 'Hook method\Exported data\Repeated data\';
path_specific_5 = 'Hook method\Exported data\Hook at Temp\';

% Segment numbers in file names for correct cooling segments
x = [4 8 12 16 20 24 28 32 4];

% Example of file name for samples 1-4: 3_downto_down_28.txt
% Example of file name for samples 5: downto_down_28.txt
s3 = '.txt';

% Sample 1
for i = 1:length(x)
    if i == 9
        s1 = '1_upto_down_'; 
    else
        s1 = '1_downto_down_'; 
    end
    s2 = num2str(x(i));
    filename = [path_main path_specific_1234 [s1 s2 s3]];
    [T,Q,~,~] = importFSC_var(filename);  
    T_PA_1(:,i) = T((end-L+1):end);
    Q_PA_1(:,i) = Q((end-L+1):end);  
end

% Sample 2
for i = 1:length(x)
    if i == 9
        s1 = '2_upto_down_'; 
    else
        s1 = '2_downto_down_'; 
    end
    s2 = num2str(x(i));
    filename = [path_main path_specific_1234 [s1 s2 s3]];
    [T,Q,~,~] = importFSC_var(filename);  
    T_PA_2(:,i) = T((end-L+1):end);
    Q_PA_2(:,i) = Q((end-L+1):end);  
end

% Sample 3
for i = 1:length(x)
    if i == 9
        s1 = '3_upto_down_'; 
    else
        s1 = '3_downto_down_'; 
    end
    s2 = num2str(x(i));
    filename = [path_main path_specific_1234 [s1 s2 s3]];
    [T,Q,~,~] = importFSC_var(filename);  
    T_PA_3(:,i) = T((end-L+1):end);
    Q_PA_3(:,i) = Q((end-L+1):end);  
end

% Sample 4
for i = 1:length(x)
    if i == 9
        s1 = '4_upto_down_'; 
    else
        s1 = '4_downto_down_'; 
    end
    s2 = num2str(x(i));
    filename = [path_main path_specific_1234 [s1 s2 s3]];
    [T,Q,~,~] = importFSC_var(filename);  
    T_PA_4(:,i) = T((end-L+1):end);
    Q_PA_4(:,i) = Q((end-L+1):end);  
end

% Sample 5
for i = 1:length(x)
    if i == 9
        s1 = 'upto_down_'; 
    else
        s1 = 'downto_down_'; 
    end
    s2 = num2str(x(i));
    filename = [path_main path_specific_5 [s1 s2 s3]];
    [T,Q,~,~] = importFSC_var(filename);  
    T_PA_5(:,i) = T((end-L+1):end);
    Q_PA_5(:,i) = Q((end-L+1):end);  
end


%% Close some variables before saving workspace

clear('filename','i','L','path_main','path_specific_1234','path_specific_5','Q','T','x','TA','TB','TC','TD','TE','s1','s2','s3')

















