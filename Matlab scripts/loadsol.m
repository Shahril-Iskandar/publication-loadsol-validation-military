clear
clc
close all
%% Change filename here
filename = 'P8 S3 T6 Flat 25kg 5min_22-08-23 18-15-41-203.txt';
%% Read participant info excel to get bodyweight

% Set directory to where the 'Participant info.xlsx' is located
cd 'C:\Users\14000\Downloads\Exoskeletal study\Data (Dynamic only)\Loadsol\new'

% Read Sheet1
info = readtable('Participant info.xlsx'); 

% Change directory to where the loadsol data are lcoated
% cd 'C:\Users\14000\Downloads\Exoskeletal study\Data (Dynamic only)\Loadsol\new'

% Get bodymass by matching the participant ID with info table
BodyMass = info{str2double(filename(2)),3};

aRate = 100;
%% Read force data
data = readmatrix(filename); % Read .txt file
time = data(:,1); % All rows, in the 1st column
leftforce = data(:,2); % All rows, in the 2nd column
rightforce = data(:,4); % All rows, in the 4th column

% Check for NaN values at the end

leftnan = sum(isnan(leftforce));
rightnan = sum(isnan(rightforce));

leftforce = leftforce(1:length(leftforce)-max(leftnan, rightnan));
rightforce = rightforce(1:length(rightforce)-max(leftnan, rightnan));
time = time(1:length(time)-max(leftnan, rightnan));
%% Offset if necessary

if strcmp(filename,'P1 S2 T1 Incline 35kg 9min_22-06-02 14-56-40-724.txt')
    offsetvalue = 150;
    leftforce = leftforce - offsetvalue; % Add the offset value to existing value  
elseif strcmp(filename,'P1 S3 T1 Incline 35kg 9min_22-06-06 10-20-28-511.txt')
    offsetvalue = 120;
    leftforce = leftforce - offsetvalue;
end
%% Apply filter to force data
% [b,a] = butter(4,25/(aRate/2),'low'); % 25Hz low-pass Butterworth-Filter (4th order)
% % filter left
% leftforce = filtfilt(b,a,leftforce); 
% % filter right
% rightforce = filtfilt(b,a,rightforce);
%% Assign treadmill speed at m/s for calculation later
if filename(10) == 'F' % Check if filename is 'Flat'
    speed = 1.11; % 4km/h
elseif filename(10) == 'I' % Check if filaname is 'Incline'
    speed = 0.56; % 2km/h
else % Check if filaname is 'Decline'
    speed = 0.86;% 3km/h
end
%% Loop through the whole raw force data

% Find index of touchdown and toeoff 
thre = 40; % Entre VGRF threshold here, e.g. 10N, 20N

% Left leg
% Getting index of touchdown and toeoff
tdownL = [];
toffL = [];
for i = 3:length(leftforce)
    if i == length(leftforce)-2
        break
    end
    if (leftforce(i) <= thre) && (leftforce(i+1) >= thre) && (leftforce(i+2) >= thre) && (leftforce(i-1) <= thre)
        tdownL(end+1) = i + 1;
    end
    if (leftforce(i) >= thre) && (leftforce(i+1) <= thre) && (leftforce(i+2) <= thre) && (leftforce(i-1) >= thre) && length(tdownL)>=length(toffL)
        toffL(end+1) = i;
        if strcmp(filename,'P3 S3 T6 Incline 35kg 9min_22-06-14 11-50-48-354.txt')
            toffL = toffL(toffL~=1610); % Remove the extra toff detected near touchdown
        elseif strcmp(filename,'P6 S3 T3 Incline 25kg 9min_22-07-14 11-49-10-261.txt')
            toffL = toffL(toffL~=391);
            toffL = toffL(toffL~=3052);
        end
    end
end

% Right leg
% Getting index of touchdown and toeoff
tdownR = [];
toffR = [];

for i = 3:length(rightforce)
    if i == length(rightforce)-2
        break
    end
    if (rightforce(i) <= thre) && (rightforce(i+1) >= thre) && (rightforce(i+2) >= thre) && (rightforce(i-1) <= thre) && (rightforce(i-2) <= thre)
        tdownR(end+1) = i + 1;
    end
    if (rightforce(i) >= thre) && (rightforce(i+1) <= thre) && (rightforce(i+2) <= thre) && (rightforce(i-1) >= thre) && length(tdownR)>=length(toffR)
        toffR(end+1) = i;
        if strcmp(filename,'P3 S2 T3 Incline 25kg 9min_22-06-09 11-58-09-067.txt')
            toffR = toffR(toffR~=1164); % Remove the extra toff detected near touchdown
        elseif strcmp(filename, 'P7 S2 T1 Incline 35kg 9min_22-06-27 10-09-08-958.txt')
            toffR = toffR(toffR~=2491);
        elseif strcmp(filename, 'P2 S2 T6 Incline 35kg 9min_22-06-09 14-31-00-024.txt')
            toffR = toffR(toffR~=1585);
        elseif strcmp(filename, 'P8 S2 T1 Incline 35kg 9min_22-08-17 15-26-55-032.txt')
            toffR = toffR(toffR~=281);
            toffR = toffR(toffR~=476);
        end
    end
end

% Remove data if touchdown/takeoff is detected again within 5 frames
for i=1:length(tdownL)
    if i==length(tdownL)
        break
    end
    if ((tdownL(i+1) - tdownL(i)) < 10)
        tdownL(i+1) = [];
    end
end

for i=1:length(tdownR)
    if i==length(tdownR)
        break
    end
    if ((tdownR(i+1) - tdownR(i)) < 10)
        tdownR(i+1) = [];
    end
end

for i=1:length(toffL)
    if i==length(toffL)
        break
    end
    if ((toffL(i+1) - toffL(i)) < 10)
        toffL(i+1) = [];
    end
end

for i=1:length(toffR)
    if i==length(toffR)
        break
    end
    if ((toffR(i+1) - toffR(i)) < 10)
        toffR(i+1) = [];
    end
end

% Fixing missing touchdown/takeoff
if strcmp(filename,'P3 S2 T3 Incline 25kg 9min_22-06-09 11-58-09-067.txt')
    tdownR=sort([tdownR 426]);
elseif strcmp(filename, 'P6 S3 T3 Incline 25kg 9min_22-07-14 11-49-10-261.txt')
    tdownL=sort([tdownL 1312]);
elseif strcmp(filename, 'P7 S3 T5 Decline 25kg 9min_22-06-30 11-37-04-250.txt')
    tdownL=sort([tdownL 533]);
elseif strcmp(filename, 'P7 S3 T2 Decline 35kg 9min_22-06-30 10-26-04-740.txt')
    toffL=sort([toffL 1277]);
% elseif strcmp(filename, 'P8 S2 T1 Incline 35kg 9min_22-08-17 15-26-55-032')
%     tdownR=sort([tdownR 411]);
end

% Taking data from touchdown to toeoff as some data may start from toeoff
if toffL(1) < tdownL(1) % Force > threshold from start
    toffL(1) = []; % Remove first takeoff
end

if length(tdownL) > length(toffL) % Touchdown but no takeoff
    tdownL(end) = []; % Remove last touchdown
end
 
% Taking data from touchdown to toeoff as some data may start from toeoff
if toffR(1) < tdownR(1) % Force > threshold from start
    toffR(1) = []; % Remove first takeoff
end

if length(tdownR) > length(toffR) % Touchdown but no takeoff
    tdownR(end) = []; % Remove last touchdown
end

% Plot graph
figure('Name', 'Both legs raw data')
tiledlayout(2,1)
nexttile
plot(leftforce); hold on; plot(tdownL, leftforce(tdownL), 'go')
plot(toffL, leftforce(toffL), 'ro'); title('Left leg'); hold off;
nexttile
plot(rightforce); hold on; plot(tdownR, rightforce(tdownR), 'go')
plot(toffR, rightforce(toffR), 'ro'); title('Right leg'); hold off;
% i=0;
% %% Check if peaks are identified correctly for all steps
% close all
% for i=1:length(tdownL)
%     findpeaks(leftforce(tdownL(i):toffL(i)),'MinPeakDistance', 1, 'MinPeakHeight', 600, 'SortStr', 'none');
%     title(sprintf('Left %d',i));
%     pause(0.5);
% end
% 
% for i=1:length(tdownR)
%     findpeaks(rightforce(tdownR(i):toffR(i)),'MinPeakDistance', 1, 'MinPeakHeight', 700, 'SortStr', 'none');
%     title(sprintf('Right %d',i));
%     pause(0.5);
% end
%% columns to remove

columns_to_remove_L = [];
columns_to_remove_R = [];

tdownL(:, columns_to_remove_L) = [];
toffL(:, columns_to_remove_L) = [];
tdownR(:, columns_to_remove_R) = [];
toffR(:, columns_to_remove_R) = [];

%% Plotting graph of each stance phase
figure('Name','Stance Phase of Each Leg')
tiledlayout(1,2)
nexttile
% Left leg
stancetimeL = [];
for i=1:length(tdownL)
    stancetimeL(end+1) = (toffL(i) - tdownL(i))/aRate;
    plot(leftforce(tdownL(i):toffL(i)))
    hold on
    if i == length(tdownL)
        ylabel('Force (N)'); xlabel('No. of frames'); title('Left leg');
        hold off
    end
end
nexttile
% Right leg
stancetimeR = [];
for i=1:length(tdownR)
    stancetimeR(end+1) = (toffR(i) - tdownR(i))/aRate;
    plot(rightforce(tdownR(i):toffR(i)))
    hold on
    if i == length(tdownR)
        ylabel('Force (N)'); xlabel('No. of frames'); title('Right leg');
        hold off
    end
end

% Statistic
mean_stanceL = mean(stancetimeL);
mean_stanceR = mean(stancetimeR);
%% Find swing time = toeoff to touchdown
% Left
swingtime_L = [];
for i=1:length(tdownL)
   if (i == length(tdownL))
       break
   end
   swingtime_L(end+1) = time(tdownL(i+1)) - time(toffL(i));
end

% Right
swingtime_R = [];
for i=1:length(tdownR)
   if (i == length(tdownR))
       break
   end
   swingtime_R(end+1) = time(tdownR(i+1)) - time(toffR(i));
end
%% Compute time of each step (stance + swing) = touchdown to touchdown
% Left
steptime_L = [];
for i=1:length(tdownL)
   if (i == length(tdownL))
       break
   end
   steptime_L(end+1) = time(tdownL(i+1)) - time(tdownL(i));
end

% Right
steptime_R = [];
for i=1:length(tdownR)
   if (i == length(tdownR))
       break
   end
   steptime_R(end+1) = time(tdownR(i+1)) - time(tdownR(i));
end
%% Identify stance time phase outliers and remove them
disp('In Progress - Checking for outlier')
% Left leg
% Get index of outlier from stance time
% Detecting outlier based on mean more than 4 SD
outindexstanceL = find(isoutlier(stancetimeL, 'median', 'ThresholdFactor', 3.5));

% Remove touchdown and toeoff outliers
tdownL(outindexstanceL) = [];
toffL(outindexstanceL) = [];

% Plot left step
figure('Name','Every stance (after deleting)')
tiledlayout(1,2)
nexttile
for i=1:length(tdownL)
    plot(leftforce(tdownL(i):toffL(i)))
    hold on
    if i == length(tdownL)
        ylabel('Force (N)'); xlabel('No. of frames'); title('Left leg');
        hold off
    end
end

% Right leg

% Get index of outlier from stance time
outindexstanceR = find(isoutlier(stancetimeR, 'median', 'ThresholdFactor', 3.5));

% Remove touchdown and toeoff outliers
tdownR(outindexstanceR) = [];
toffR(outindexstanceR) = [];

% Plot right step
nexttile
for i=1:length(tdownR)
    plot(rightforce(tdownR(i):toffR(i)))
    hold on
    if i == length(tdownR)
        ylabel('Force (N)'); xlabel('No. of frames'); title('Right leg')
        hold off
    end
end
%% Identify step time outlier and remove

% Left leg

% Method 1: Using the same index identified in the stance time
% If stance outlier is remove, remove step time accordingly
if outindexstanceL <= length(steptime_L) % Last index has been removed
    steptime_L(outindexstanceL) = [];
end

% Right leg

% Method 1: Using the same index identified in the stance time
if outindexstanceR <= length(steptime_R) % Last index has been removed
    steptime_R(outindexstanceR) = [];
end

% Statistics
mean_steptimeL = mean(steptime_L);
mean_steptimeR = mean(steptime_R);

%% Update stance time after removing outliers

% Left leg
if (isempty(outindexstanceL)) ~= 1 && (length(outindexstanceL(end)) == length(stancetimeL))
    outindexstanceL(end) = []; % Remove the last index because less 1 index
    stancetimeL(outindexstanceL) = [];
else
    stancetimeL(outindexstanceL) = [];
end

% Right leg
if (isempty(outindexstanceR)) ~= 1 && (length(outindexstanceR(end)) == length(stancetimeR))
    outindexstanceR(end) = []; % Remove the last index
    stancetimeR(outindexstanceR) = [];
else
    stancetimeR(outindexstanceR) = [];
end

% Statistic
mean_stanceL = mean(stancetimeL);
mean_stanceR = mean(stancetimeR);
%% Update swing time 

% Left leg
% Using the same index identified in the stance time
if outindexstanceL <= length(swingtime_L) % Last index removed
    swingtime_L(outindexstanceL) = [];
end

% Right leg
if outindexstanceR <= length(swingtime_R)
    swingtime_R(outindexstanceR) = []; % Last index removed
end
disp('Completed - Outlier(s) removed')
% Statistic
mean_swingL = mean(swingtime_L);
mean_swingR = mean(swingtime_R);

stepcountL = length(steptime_L);
stepcountR = length(steptime_R);
%% Find peaks in each step and loading rate (new method)
disp('In Progress - Calculating loading rate')
% Left
fz_peak1_L = []; % Impact peak
fz_peak2_L = []; % Active peak
time_to_impact_peak_L = [];
loadingrate_L = [];

for i=1:length(tdownL)
%     disp(i);
    [pksL,locsL] = findpeaks(leftforce(tdownL(i):toffL(i)),'MinPeakDistance', 1, 'MinPeakHeight', 600, 'SortStr', 'none');
    fz_peak1_L(end+1) = pksL(1); % value of 1st peak force, impact peak
    if length(pksL) >= 2
        fz_peak2_L(end+1) = pksL(end); % value of last peak force, active peak
    end
    time_to_impact_peak_L(end+1) = locsL(1)/aRate;
    % Loading rate based on 20-80% gradient,  because already index tdown to toff, cannot find value directly based on the index
    loadingrate_L(end+1) = (leftforce(tdownL(i) + round(0.8*locsL(1)) - 1) - leftforce(tdownL(i) + round(0.2*locsL(1)) - 1)) / ((0.8*locsL(1)/aRate)-(0.2*locsL(1)/aRate));
end

% Right
fz_peak1_R = [];
fz_peak2_R = [];
time_to_impact_peak_R = [];
loadingrate_R = [];

for i=1:length(tdownR)
%     disp(i);
    [pksR,locsR] = findpeaks(rightforce(tdownR(i):toffR(i)),'MinPeakDistance', 1, 'MinPeakHeight', 700, 'SortStr', 'none');
    fz_peak1_R(end+1) = pksR(1); % value of 1st peak force
    if length(pksR) >= 2
        fz_peak2_R(end+1) = pksR(end); % value of 2nd peak force
    end
    time_to_impact_peak_R(end+1) = locsR(1)/aRate;
    % Loading rate based on 20-80% gradient
    loadingrate_R(end+1) = (rightforce(tdownR(i) + round(0.8*locsR(1)) - 1) - rightforce(tdownR(i) + round(0.2*locsR(1)) - 1)) / ((0.8*locsR(1)/aRate)-(0.2*locsR(1)/aRate));
end
disp('Completed - Loading rate calculated')
% Statistics
mean_impact_peak_L = mean(fz_peak1_L);
mean_impact_peak_R = mean(fz_peak1_R);
mean_active_peak_L = mean(fz_peak2_L);
mean_active_peak_R = mean(fz_peak2_R);
mean_loadingrate_L = mean(loadingrate_L);
mean_loadingrate_R = mean(loadingrate_R);
%% Normalised Fz

% Left leg
for i=1:length(tdownL)
    x = time(tdownL(i):toffL(i));
    y = [leftforce(tdownL(i):toffL(i))]';
    yy(i,:) = spline(x,[0 y 0]);
    xx(i,:) = linspace(time(tdownL(i)),time(toffL(i)),101);  % create 101 points (0-100%) 
    Fz_L_norm(i,:) = ppval(yy(i,:),xx(i,:)); % rename the normliased Fz_1 to Fz_1_norm
end

% Normalise Fz to Body Weight for graph plotting
Fz_L_BW = Fz_L_norm/BodyMass/9.81;

% Calcualte avearege force-time history of all rearfoot steps
Fz_L_mean = mean(Fz_L_BW,'omitnan');

% To show plot on its own:
% figure('Name','All steps - Left')
% for i = 2:length(tdownL)
%        plot(Fz_L_BW(i,:),'b'); ylabel('Force (Body Weight)'); xlabel('Stance Phase (%)'); hold on
% end

% Right leg
for i=1:length(tdownR)
    x = time(tdownR(i):toffR(i));
    y = [rightforce(tdownR(i):toffR(i))]';
    yy(i,:) = spline(x,[0 y 0]);
    xx(i,:) = linspace(time(tdownR(i)),time(toffR(i)),101);  % create 101 points (0-100%) 
    Fz_R_norm(i,:) = ppval(yy(i,:),xx(i,:)); % rename the normliased Fz_1 to Fz_1_norm
end

% Normalise Fz to Body Weight for graph plotting
Fz_R_BW = Fz_R_norm/BodyMass/9.81;

% Calcualte avearege force-time history of all rearfoot steps
Fz_R_mean = mean(Fz_R_BW,'omitnan');

% To show plot on its own:
% figure('Name','All steps - Right')
% for i = 2:length(tdownR)
%        plot(Fz_R_BW(i,:),'b'); ylabel('Force (Body Weight)'); xlabel('Stance Phase (%)'); hold on
% end

% Fz
figure('Name','Normalised Mean GRF')
plot(Fz_L_mean, 'b','LineWidth',2); title('Normalised Mean Vertical GRF', 'FontSize', 20)
ylabel('Force (Body Weight)','FontSize',16); xlabel('Stance Phase (%)','FontSize',16);
xlim([0 101]);
hold On
plot(Fz_R_mean, ':r','LineWidth',2);
ylabel('Force (Body Weight)','FontSize',16); xlabel('Stance Phase (%)','FontSize',16);
legend('Left leg', 'Right leg')
hold Off

%% Store stance data into array
% Left leg
for i=1:length(tdownL)
    stanceFzL{i} = leftforce(tdownL(i):toffL(i)); % Store each fz stance data into a cell array
end
% Right leg
for i=1:length(tdownR)
    stanceFzR{i} = rightforce(tdownR(i):toffR(i)); % Store each fz stance data into a cell array
end
%% Store raw data of all stance phase into Excel by column, after removing crossover
% Name the file according to participant ID
stancefilename = sprintf('%s_Stance_Raw_Data.xlsx', filename(1:26));
disp('In Progress - Saving raw stance data (after removing crossover) into Excel')
% Left stance
j = 65; % char(65) = A
if length(stanceFzL) >= 53 % char(90) = Z
    for i = 1:26
        str1 = char(j);
        writematrix(stanceFzL{1,i}, stancefilename, 'Sheet', 'Fz_Left', 'Range', append(str1,num2str(1)));
        j = j + 1;
    end
    j = 65; % Reset back to A, Excel columns after Z is AA, AB, AC...
    for i = 27:52
        str1 = char('A');
        str2 = char(j);
        writematrix(stanceFzL{1,i}, stancefilename, 'Sheet', 'Fz_Left', 'Range', append(str1,str2,num2str(1)));
        j = j + 1;    
    end
    j = 65; % Reset back to A, Excel columns after Z is AA, AB, AC...
    for i = 53:length(stanceFzL)
        str1 = char('B');
        str2 = char(j);
        writematrix(stanceFzL{1,i}, stancefilename, 'Sheet', 'Fz_Left', 'Range', append(str1,str2,num2str(1)));
        j = j + 1;          
    end
elseif (length(stanceFzL) >= 26) && (length(stanceFzL) <= 52) 
    for i = 1:26
        str1 = char(j);
        writematrix(stanceFzL{1,i}, stancefilename, 'Sheet', 'Fz_Left', 'Range', append(str1,num2str(1)));
        j = j + 1;
    end
    j = 65; % Reset back to A, Excel columns after Z is AA, AB, AC...
    for i = 27:length(stanceFzL)
        str1 = char('A');
        str2 = char(j);
        writematrix(stanceFzL{1,i}, stancefilename, 'Sheet', 'Fz_Left', 'Range', append(str1,str2,num2str(1)));
        j = j + 1;    
    end    
else
    for i = 1:length(stanceFzL)
        str1 = char(j);
        writematrix(stanceFzL{1,i}, stancefilename, 'Sheet', 'Fz_Left', 'Range', append(str1,num2str(1)));
        j = j + 1;
    end
end
disp('Completed - Saved left raw stance data into Excel')

% Right stance
j = 65; % char(j) = A
if length(stanceFzR) >= 53 % char(90) = Z
    for i = 1:26
        str1 = char(j);
        writematrix(stanceFzR{1,i}, stancefilename, 'Sheet', 'Fz_Right', 'Range', append(str1,num2str(1)));
        j = j + 1;
    end
    j = 65; % Reset back to A, Excel columns after Z is AA, AB, AC...
    for i = 27:52
        str1 = char('A');
        str2 = char(j);
        writematrix(stanceFzR{1,i}, stancefilename, 'Sheet', 'Fz_Right', 'Range', append(str1,str2,num2str(1)));
        j = j + 1;
    end
    j = 65; % Reset back to A, Excel columns after Z is AA, AB, AC...
    for i = 53:length(stanceFzR)
        str1 = char('B');
        str2 = char(j);
        writematrix(stanceFzR{1,i}, stancefilename, 'Sheet', 'Fz_Right', 'Range', append(str1,str2,num2str(1)));
        j = j + 1;          
    end
elseif (length(stanceFzR) >= 26) && (length(stanceFzR) <= 52)
    for i = 1:26
        str1 = char(j);
        writematrix(stanceFzR{1,i}, stancefilename, 'Sheet', 'Fz_Right', 'Range', append(str1,num2str(1)));
        j = j + 1;
    end
    j = 65; % Reset back to A, Excel columns after Z is AA, AB, AC...
    for i = 27:length(stanceFzR)
        str1 = char('A');
        str2 = char(j);
        writematrix(stanceFzR{1,i}, stancefilename, 'Sheet', 'Fz_Right', 'Range', append(str1,str2,num2str(1)));
        j = j + 1;
    end    
else
    for i = 1:length(stanceFzR)
        str1 = char(j);
        writematrix(stanceFzR{1,i}, stancefilename, 'Sheet', 'Fz_Right', 'Range', append(str1,num2str(1)));
        j = j + 1;
    end       
end
disp('Completed - Saved raw stance data into Excel')
%% Store mean data into Excel
disp('In Progress - Saving mean stance data into Excel')
% Create Mean excel data is don't have yet
% Append data to next row

a = filename(1:2); % Get participant ID
% b = filename(4:5); % Get session number

% Get from Participant Info if condition is with Exos or without
if filename(5) == '2'
    Exo = info{str2double(filename(2)),5}; % Session 2
else
    Exo = info{str2double(filename(2)),6}; % Session 3
end

c = filename(7:8); % Get trial number

if filename(10) == 'F' % Check if filename is 'Flat'
    d = filename(10:13); % Get full name
    Load = filename(15:18); % Get the load in kg
elseif filename(10) == 'I' % Check if filaname is 'Incline'
    d = filename(10:17); % Get full name
    Load = filename(18:21); % Get the load in kg
else % Check if filaname is 'Decline'
    d = filename(10:16); % Get full name
    Load = filename(18:21); % Get the load in kg
end

norm_fz_L = string(Fz_L_mean);
norm_fz_R = string(Fz_R_mean);
meanfzdata_L = [a Exo d Load norm_fz_L];
meanfzdata_R = [a Exo d Load norm_fz_R];
writematrix(meanfzdata_L, 'Mean_Data_All_Participant_Updated.xlsx', 'Sheet', 'Fz_Left', 'WriteMode', 'append')
writematrix(meanfzdata_R, 'Mean_Data_All_Participant_Updated.xlsx', 'Sheet', 'Fz_Right', 'WriteMode', 'append')
disp('Completed - Saved mean stance data into Excel')

%% Impulse in each stance 
% Left
fz_impulse_L = [];
for i=1:length(tdownL)
    fz_impulse_L(end+1) = trapz(leftforce(tdownL(i):toffL(i)))/aRate;
end

% Right
fz_impulse_R = [];
for i=1:length(tdownR)
    fz_impulse_R(end+1) = trapz(rightforce(tdownR(i):toffR(i)))/aRate;
end

% Statistic
mean_impulseL = mean(fz_impulse_L);
mean_impulseR = mean(fz_impulse_R);

%% Calculate cadence and stride lengths
cadence = (stepcountL+stepcountR)/time(end)*60;
step_frequency = cadence/60;
step_length = speed/step_frequency;
stride_length = step_length*2;

%%  Write summary table for outputs
Participant_ID = filename(1:2);

% Dec = Decline, Inc = Incline, Fla = Flat
if filename(10:12) == "Fla"
    Condition = "Flat";
elseif filename(10:12) == "Inc"
    Condition = "Incline";
else
    Condition = "Decline";
end

% % Get from Participant Info if condition is with Exos or without
% if filename(5) == '2'
%     Exo = info{str2double(filename(2)),5}; % Session 2
% else
%     Exo = info{str2double(filename(2)),6}; % Session 3
% end

participant = {Participant_ID, Exo, Condition, Load};
T = cell2table(participant, 'VariableNames', ["ID" "Exo" "Condition" "Load"]);

Stance_Time_L = [mean_stanceL];
Stance_Time_R = [mean_stanceR];
Step_Time_L = [mean_steptimeL];
Step_Time_R = [mean_steptimeR];
Impact_PeakForce_L = [mean_impact_peak_L];
Impact_PeakForce_R = [mean_impact_peak_R];
Active_PeakForce_L = [mean_active_peak_L];
Active_PeakForce_R = [mean_active_peak_R];
Impulse_L = [mean_impulseL];
Impulse_R = [mean_impulseR];
Loading_Rate_L = [mean_loadingrate_L];
Loading_Rate_R = [mean_loadingrate_R];
Cadence = [cadence];
Step_Frequency = [step_frequency];
Step_Length = [step_length];
Stride_Length = [stride_length];

T1 = table(Stance_Time_L, Stance_Time_R, Step_Time_L, Step_Time_R, ...
Impact_PeakForce_L, Impact_PeakForce_R, Active_PeakForce_L, Active_PeakForce_R, ...
Impulse_L, Impulse_R, Loading_Rate_L, Loading_Rate_R, ...
Cadence, Step_Frequency, Step_Length, Stride_Length);
T1 = [T,T1];
T1

% Save T1 to Excel without overwriting previous saved data in the same sheet
writetable(T1, 'Loadsol GRF variables-updated.xlsx', 'WriteMode', 'Append')
disp('Completed - Saved Loadsol GRF variables into Excel')