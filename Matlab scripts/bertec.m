clear
clc
close all
%%
cd 'C:\Users\14000\Downloads\Exoskeletal study\Data (Dynamic only)'

% [To change filename accordingly]
filename = 'P1_S3_T4_Incline_25kg_9min.c3d'; % Rename file with c3d at the end
disp(filename);
% emgnormfile = 'P1_S2_EMG_Norm_4min.c3d'; % If no emg data, leave blank
%========================================================================
% fps_fz = 1000; % frequency for grf data
% fps_emg = 2000; % frequency for emg data
info = readtable('Participant info.xlsx'); % Read Sheet1
BodyMass = info{str2double(filename(2)),3}; % Check ID in the filename entered and get the corresponding bodymass accordingly
disp('In Progress - Reading force data in c3d file')

% Import data using c3d file
c3d = c3dserver;
pRet = c3d.Open(filename,3); % Load File in C3D Server
% nFirst = c3d.GetVideoFrame(0); % First Video Frame
% nLast = c3d.GetVideoFrame(1); % Last Video Frame
% nFrame = c3d.nframes(); % Number of frames
vRate = c3d.GetVideoFrameRate(); % Video sampling frequenzy
nRatio = c3d.GetAnalogVideoRatio(); % Sampling ratio Video/Force
aRate = vRate * nRatio; % Sampling frequenzy analoge data
fprintf('c3d frequency: %s Hz \n', num2str(aRate))
% aRate = 2000;
nAnalog = c3d.GetAnalogChannels; % Number of used analog Channels

% Getting force platform data,
% Check if file data is for 'Decline' from filename
% Decline = participant walking in opposite direction, meaning belt 1 will
% be right leg instead of left leg, vice versa
if filename(10) == 'D'
    fz_right = double(getanalogchannel(c3d,'Force.Fz1'))*-1;
    fx_right = double(getanalogchannel(c3d,'Force.Fx1'))*-1;
    fy_right = double(getanalogchannel(c3d,'Force.Fy1'))*-1;
    fz_left = double(getanalogchannel(c3d,'Force.Fz2'))*-1;
    fx_left = double(getanalogchannel(c3d,'Force.Fx2'))*-1;
    fy_left = double(getanalogchannel(c3d,'Force.Fy2'))*-1;
else 
    fz_left = double(getanalogchannel(c3d,'Force.Fz1'))*-1;
    fx_left = double(getanalogchannel(c3d,'Force.Fx1'))*-1;
    fy_left = double(getanalogchannel(c3d,'Force.Fy1'))*-1;
    fz_right = double(getanalogchannel(c3d,'Force.Fz2'))*-1;
    fx_right = double(getanalogchannel(c3d,'Force.Fx2'))*-1;
    fy_right = double(getanalogchannel(c3d,'Force.Fy2'))*-1;    
end

% If c3d is in 2000 Hz, change back to 1000 Hz
% To standardised interpolation using MATLAB
if (str2double(filename(2)) == 3 || str2double(filename(2)) == 4 || str2double(filename(2)) == 7) && (aRate == 2000) % If id is 3,4 or 7, and is 2000 Hz, run code
    fz_right = fz_right(1:2:end);
    fx_right = fx_right(1:2:end);
    fy_right = fy_right(1:2:end);
    fz_left = fz_left(1:2:end);
    fx_left = fx_left(1:2:end);
    fy_left = fy_left(1:2:end);
    disp("Data was changed to 1000 Hz");
    aRate = 1000;
end

% Assign treadmill speed at m/s for calculation later
if filename(10) == 'F' % Check if filename is 'Flat'
    speed = 1.11; % 4km/h
elseif filename(10) == 'I' % Check if filaname is 'Incline'
    speed = 0.56; % 2km/h
else % Check if filaname is 'Decline'
    speed = 0.86;% 3km/h
end

% Interpolate force data for participants without EMG data
if (str2double(filename(2)) == 3 || str2double(filename(2)) == 4 || str2double(filename(2)) == 7) && (aRate == 1000) % If id is 3,4 or 7, and is 1000 Hz, run code
    filteredfz_left = interp(fz_left,2); % Interpolate by twice its frequency
    filteredfx_left = interp(fx_left,2)*-1;% Changing the values to match direction
    filteredfy_left = interp(fy_left,2);
    filteredfz_right = interp(fz_right,2);
    filteredfx_right = interp(fx_right,2);
    filteredfy_right = interp(fy_right,2);
    disp("Data was interpolated");
else
    filteredfz_left = fz_left; % Just replace the name but haven't filter
    filteredfx_left = fx_left*-1; % Changing the values to match direction
    filteredfy_left = fy_left;
    filteredfz_right = fz_right;
    filteredfx_right = fx_right;
    filteredfy_right = fy_right ;    
    disp("Data was NOT interpolated");
end

aRate = 2000; % Change the frequency to 2000

% Apply filter to force data
[b,a] = butter(4,100/(aRate/2),'low'); % 100Hz low-pass Butterworth-Filter (4th order)

% filter left
filteredfz_left = filtfilt(b,a,filteredfz_left); 
filteredfx_left = filtfilt(b,a,filteredfx_left); 
filteredfy_left = filtfilt(b,a,filteredfy_left); 
% filter right
filteredfz_right = filtfilt(b,a,filteredfz_right);
filteredfx_right = filtfilt(b,a,filteredfx_right); 
filteredfy_right = filtfilt(b,a,filteredfy_right);

% Compute time
increment = 1/aRate;
starttime = 0; timeend = length(filteredfz_left);
totaltime = 0:increment:(timeend-1)*increment;
fprintf('Total time: %s seconds\n', num2str(totaltime(end))) % Show total time in seconds to confirm
% totaltime = starttime + increment*(0:timeend-1);
disp('Completed - Extracted force data from c3d file')

%% Offset files affected
disp('In Progress - Checking if Fz offset is required')
% Checking if filename match
if strcmp(filename,'P1_S3_T1_Incline_35kg_9min.c3d')
    offsetvalue = 120;
    filteredfz_right = filteredfz_right + offsetvalue; % Add the offset value to existing value
    disp('Offset done');
elseif strcmp(filename, 'P2_S3_T6_Incline_35kg_9min.c3d')
    offsetvalue = 255;
    filteredfz_right = filteredfz_right + offsetvalue;
    disp('Offset done');    
elseif strcmp(filename, 'P5_S3_T2_Incline_25kg_9min.c3d')
    offsetvalue = 60;
    filteredfz_right = filteredfz_right + offsetvalue;
    disp('Offset done');
elseif strcmp(filename, 'P6_S2_T2_Decline_25kg_9min.c3d')
    offsetvalue = 60;
    filteredfz_left = filteredfz_left + offsetvalue;
    disp('Offset done');
elseif strcmp(filename, 'P8_S2_T4_Incline_25kg_9min.c3d')
    offsetvalue = 50;
    filteredfz_right = filteredfz_right + offsetvalue;
    disp('Offset done');    
end
disp('Completed - Fz offset checked')

%% Offset Fx (mediolateral)
disp('In Progress - Checking if Fx offset is required')

% Checking if filename match
if strcmp(filename,'P2_S2_T2_Decline_25kg_9min.c3d')
    offsetvalue = 60;
    filteredfx_left = filteredfx_left - offsetvalue; % Add the offset value to existing value
    disp('Offset done');
elseif strcmp(filename,'P2_S2_T3_Incline_25kg_9min.c3d')
    offsetvalue = 50;
    filteredfx_right = filteredfx_right + offsetvalue; % Add the offset value to existing value
    disp('Offset done');
elseif strcmp(filename,'P2_S2_T6_Incline_35kg_9min.c3d')
    offsetvalue = 65;
    filteredfx_right = filteredfx_right + offsetvalue; % Add the offset value to existing value
    disp('Offset done');
elseif strcmp(filename,'P3_S2_T2_Decline_25kg_9min.c3d')
    offsetvalue = 80;
    filteredfx_left = filteredfx_left - offsetvalue; % Add the offset value to existing value
    disp('Offset done');
elseif strcmp(filename,'P3_S2_T3_Incline_25kg_9min.c3d')
    offsetvalue = 45;
    filteredfx_right = filteredfx_right + offsetvalue; % Add the offset value to existing value
    disp('Offset done');
elseif strcmp(filename,'P3_S2_T5_Decline_35kg_9min.c3d')
    offsetvalue = 60;
    filteredfx_left = filteredfx_left - offsetvalue; % Add the offset value to existing value
    disp('Offset done');
elseif strcmp(filename,'P3_S2_T6_Incline_35kg_9min.c3d')
    offsetvalue = 80;
    filteredfx_right = filteredfx_right + offsetvalue; % Add the offset value to existing value
    disp('Offset done');    
elseif strcmp(filename,'P3_S3_T2_Decline_25kg_9min.c3d')
    offsetvalue = 50;
    filteredfx_left = filteredfx_left - offsetvalue; % Add the offset value to existing value
    disp('Offset done');
elseif strcmp(filename,'P3_S3_T4_Flat_35kg_9min.c3d')
    offsetvalue = 70;
    filteredfx_right = filteredfx_right + offsetvalue; % Add the offset value to existing value
    disp('Offset done');
elseif strcmp(filename, 'P2_S3_T6_Incline_35kg_9min.c3d')
    offsetvalue = 28;
    filteredfx_right = filteredfx_right - offsetvalue; % Add the offset value to existing value
    disp('Offset done');
elseif strcmp(filename, 'P5_S3_T2_Incline_25kg_9min.c3d')
    offsetvalue = 40;
    filteredfx_right = filteredfx_right + offsetvalue; % Add the offset value to existing value
    disp('Offset done');
elseif strcmp(filename, 'P8_S2_T4_Incline_25kg_9min.c3d')
    offsetvalue = 40;
    filteredfx_right = filteredfx_right + offsetvalue; % Add the offset value to existing value
    disp('Offset done'); 
end
disp('Completed - Fx offset checked')


%%
% plot(filteredfx_right); hold on; plot(filteredfx_left); hold off; legend('Right', 'Left');
%% Loop through the whole raw force data

% Find index of touchdown and toeoff 
thre = 40; % Entre VGRF threshold here, e.g. 10N, 20N

% Left leg
% Getting index of touchdown and toeoff
tdownL = [];
toffL = [];
for i = 1:length(filteredfz_left)
    if i == length(filteredfz_left)
        break
    end
    if (filteredfz_left(i) < thre) && (filteredfz_left(i+1) > thre) && (filteredfz_left(i+2) > thre) && (filteredfz_left(i-1) < thre)
        tdownL(end+1) = i + 1;
    end
    if (filteredfz_left(i) > thre) && (filteredfz_left(i+1) < thre) && (filteredfz_left(i+2) < thre) && (filteredfz_left(i-1) > thre) && length(tdownL)>=length(toffL)
        toffL(end+1) = i;
    end
end

% Taking data from touchdown to toeoff as some data may start from toeoff
if toffL(1) < tdownL(1) % Force > threshold from start
    toffL(1) = []; % Remove first takeoff
end

if length(tdownL) > length(toffL) % Touchdown but no takeoff
    tdownL(end) = []; % Remove last touchdown
end

% Right leg
% Getting index of touchdown and toeoff
tdownR = [];
toffR = [];

for i = 1:length(filteredfz_right)
    if i == length(filteredfz_right)
        break
    end
    if (filteredfz_right(i) < thre) && (filteredfz_right(i+1) > thre) && (filteredfz_right(i+2) > thre) && (filteredfz_right(i-1) < thre)
        tdownR(end+1) = i + 1;
    end
    if (filteredfz_right(i) > thre) && (filteredfz_right(i+1) < thre) && (filteredfz_right(i+2) < thre) && (filteredfz_right(i-1) > thre) && length(tdownR)>=length(toffR)
        toffR(end+1) = i;
    end
end

% Taking data from touchdown to toeoff as some data may start from toeoff
if toffR(1) < tdownR(1) % Force > threshold from start
    toffR(1) = []; % Remove first takeoff
end

if length(tdownR) > length(toffR) % Touchdown but no takeoff
    tdownR(end) = []; % Remove last touchdown
end

% Visually check for any crossover (force did not return to 0)
figure('Name', 'Both legs raw data')
tiledlayout(2,1)
nexttile
plot(filteredfz_left); hold on; plot(tdownL, filteredfz_left(tdownL), 'go')
plot(toffL, filteredfz_left(toffL), 'ro'); title('Left leg'); hold off;
nexttile
plot(filteredfz_right); hold on; plot(tdownR, filteredfz_right(tdownR), 'go')
plot(toffR, filteredfz_right(toffR), 'ro'); title('Right leg'); hold off;
%% Plotting graph of each stance before removing crossover steps
% Can visually inspect for crossover steps
% Remember to check at touchdown area
outindexstanceL = [];
outindexstanceR = [];

figure('Name','Stance Phase of Each Leg (Before removing crossover)')
tiledlayout(1,2)
nexttile
% Left leg
stancetimeL = [];
for i=1:length(tdownL)
    if length(filteredfz_left(tdownL(i):toffL(i))) <= 300
        outindexstanceL(end+1) = i;
        continue
    end
%     if length(filteredfz_left(tdownL(i):toffL(i))) >= 2000
%         outindexstanceL(end+1) = i;
%         continue
%     end 
    stancetimeL(end+1) = (toffL(i) - tdownL(i))/aRate; % 2000 Hz
    plot(filteredfz_left(tdownL(i):toffL(i)))
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
    if length(filteredfz_right(tdownR(i):toffR(i))) <= 300
        outindexstanceR(end+1) = i;
        continue
    end
%     if length(filteredfz_right(tdownR(i):toffR(i))) >= 2000
%         outindexstanceR(end+1) = i;
%         continue
%     end    
    stancetimeR(end+1) = (toffR(i) - tdownR(i))/aRate;
    plot(filteredfz_right(tdownR(i):toffR(i)))
    hold on
    if i == length(tdownR)
        ylabel('Force (N)'); xlabel('No. of frames'); title('Right leg');
        hold off;
    end
end

% Remove touchdown and toeoff outliers
tdownL(outindexstanceL) = [];
toffL(outindexstanceL) = [];

% Remove touchdown and toeoff outliers
tdownR(outindexstanceR) = [];
toffR(outindexstanceR) = [];
%% Identify stance time phase outliers and remove them
disp('In Progress - Checking for outlier')
% Left leg
% Get index of outlier from stance time
% Detecting outlier based on mean more than 4 SD
outindexstanceL = find(isoutlier(stancetimeL, 'median', 'ThresholdFactor', 4));

% Remove touchdown and toeoff outliers
tdownL(outindexstanceL) = [];
toffL(outindexstanceL) = [];

% Plot left step
figure('Name','Every stance (after deleting)')
tiledlayout(1,2)
nexttile
for i=1:length(tdownL)
    plot(filteredfz_left(tdownL(i):toffL(i)))
    hold on
    if i == length(tdownL)
        ylabel('Force (N)'); xlabel('No. of frames'); title('Left leg');
        hold off
    end
end

% Right leg

% Get index of outlier from stance time
outindexstanceR = find(isoutlier(stancetimeR, 'median', 'ThresholdFactor', 4));

% Remove touchdown and toeoff outliers
tdownR(outindexstanceR) = [];
toffR(outindexstanceR) = [];

% Plot right step
nexttile
for i=1:length(tdownR)
    plot(filteredfz_right(tdownR(i):toffR(i)))
    hold on
    if i == length(tdownR)
        ylabel('Force (N)'); xlabel('No. of frames'); title('Right leg')
        hold off
    end
end
%% Update stance time after removing outliers

% Statistic
mean_stanceL = mean(stancetimeL);
mean_stanceR = mean(stancetimeR);
disp('Completed - Outlier(s) removed')
%% Find swing time = toeoff to touchdown
% Left
swingtime_L = [];
for i=1:length(tdownL)
   if (i == length(tdownL))
       break
   end
   swingtime_L(end+1) = totaltime(tdownL(i+1)) - totaltime(toffL(i));
end

% Right
swingtime_R = [];
for i=1:length(tdownR)
   if (i == length(tdownR))
       break
   end
   swingtime_R(end+1) = totaltime(tdownR(i+1)) - totaltime(toffR(i));
end

% Get index of outlier from stance time
% Detecting outlier based on mean more than 4 SD
outindexswingL = find(isoutlier(swingtime_L, 'median', 'ThresholdFactor', 4));
outindexswingR = find(isoutlier(swingtime_R, 'median', 'ThresholdFactor', 4));

% Remove outlier
swingtime_L(outindexswingL) = [];
swingtime_R(outindexswingR) = [];

% Statistics
mean_swingL = mean(swingtime_L);
mean_swingR = mean(swingtime_R);
%% Compute time of each step (stance + swing) = touchdown to touchdown
% Left
steptime_L = [];
for i=1:length(tdownL)
   if (i == length(tdownL))
       break
   end
   steptime_L(end+1) = totaltime(tdownL(i+1)) - totaltime(tdownL(i));
end

% Right
steptime_R = [];
for i=1:length(tdownR)
   if (i == length(tdownR))
       break
   end
   steptime_R(end+1) = totaltime(tdownR(i+1)) - totaltime(tdownR(i));
end

% Statistics
stepcountL = length(steptime_L);
stepcountR = length(steptime_R);
mean_steptimeL = mean(steptime_L);
mean_steptimeR = mean(steptime_R);

%% columns to remove

columns_to_remove_L = [];
columns_to_remove_R = [];

tdownL(:, columns_to_remove_L) = [];
toffL(:, columns_to_remove_L) = [];
tdownR(:, columns_to_remove_R) = [];
toffR(:, columns_to_remove_R) = [];
%% Find peaks in each step and loading rate (new method)
disp('In Progress - Calculating loading rate')
% Left
fz_peak1_L = []; % Impact peak
fz_peak2_L = []; % Active peak
time_to_impact_peak_L = [];
loadingrate_L = [];

for i=1:length(tdownL)
    [pksL,locsL] = findpeaks(filteredfz_left(tdownL(i):toffL(i)),'MinPeakDistance',400, 'MinPeakHeight', BodyMass*9.81, 'SortStr', 'none');
    % if first peak not measured, take note of the index
    if isempty(pksL)
        sprintf("Left %d", i)
    else
        fz_peak1_L(end+1) = pksL(1); % value of 1st peak force, impact peak
        time_to_impact_peak_L(end+1) = locsL(1)/aRate;
        % Loading rate based on 20-80% gradient,  because already index tdown to toff, cannot find value directly based on the index
        loadingrate_L(end+1) = (filteredfz_left(tdownL(i) + round(0.8*locsL(1)) - 1) - filteredfz_left(tdownL(i) + round(0.2*locsL(1)) - 1)) / ((0.8*locsL(1)/aRate)-(0.2*locsL(1)/aRate));
    end
    if length(pksL) >= 2
        fz_peak2_L(end+1) = pksL(end); % value of last peak force, active peak
    end
end

% Right
fz_peak1_R = [];
fz_peak2_R = [];
time_to_impact_peak_R = [];
loadingrate_R = [];

for i=1:length(tdownR)
    [pksR,locsR] = findpeaks(filteredfz_right(tdownR(i):toffR(i)),'MinPeakDistance',400, 'MinPeakHeight', BodyMass*9.81, 'SortStr', 'none');
    if isempty(pksR)
        sprintf("Right %d", i)
    else
        fz_peak1_R(end+1) = pksR(1); % value of 1st peak force
        time_to_impact_peak_R(end+1) = locsR(1)/aRate;
        % Loading rate based on 20-80% gradient
        loadingrate_R(end+1) = (filteredfz_right(tdownR(i) + round(0.8*locsR(1)) - 1) - filteredfz_right(tdownR(i) + round(0.2*locsR(1)) - 1)) / ((0.8*locsR(1)/aRate)-(0.2*locsR(1)/aRate));
    end
    if length(pksR) >= 2
        fz_peak2_R(end+1) = pksR(end); % value of 2nd peak force
    end
end
disp('Completed - Loading rate calculated')

%% Normalised Fz
disp('In Progress - Normalising Fz')
% Left leg
for i=1:length(tdownL)
    x = totaltime(tdownL(i):toffL(i));
    y = [filteredfz_left(tdownL(i):toffL(i))]';
    yy(i,:) = spline(x,[0 y 0]);
    xx(i,:) = linspace(totaltime(tdownL(i)),totaltime(toffL(i)),101);  % create 101 points (0-100%) 
    Fz_L_norm(i,:) = ppval(yy(i,:),xx(i,:)); % rename the normliased Fz_1 to Fz_1_norm
end

% Normalise Fz to Body Weight for graph plotting
Fz_L_BW = Fz_L_norm/BodyMass/9.81;

% Calcualte avearege force-time history of all rearfoot steps
Fz_L_mean = mean(Fz_L_BW,'omitnan');

% Unnormalised force
Fz_L_unnormalised = mean(Fz_L_norm);

% To show plot on its own:
% figure('Name','All steps - Left')
% for i = 2:length(tdownL)
%        plot(Fz_L_BW(i,:),'b'); ylabel('Force (Body Weight)'); xlabel('Stance Phase (%)'); hold on
% end

% Right leg
for i=1:length(tdownR)
    x = totaltime(tdownR(i):toffR(i));
    y = [filteredfz_right(tdownR(i):toffR(i))]';
    yy(i,:) = spline(x,[0 y 0]);
    xx(i,:) = linspace(totaltime(tdownR(i)),totaltime(toffR(i)),101);  % create 101 points (0-100%) 
    Fz_R_norm(i,:) = ppval(yy(i,:),xx(i,:)); % rename the normliased Fz_1 to Fz_1_norm
end

% Normalise Fz to Body Weight for graph plotting
Fz_R_BW = Fz_R_norm/BodyMass/9.81;

% Calcualte avearege force-time history of all rearfoot steps
Fz_R_mean = mean(Fz_R_BW,'omitnan');

% Unnormalised force
Fz_R_unnormalised = mean(Fz_R_norm);

% To show plot on its own:
% figure('Name','All steps - Right')
% for i = 2:length(tdownR)
%        plot(Fz_R_BW(i,:),'b'); ylabel('Force (Body Weight)'); xlabel('Stance Phase (%)'); hold on
% end
disp('Completed - Normalising Fz')
%% Calculate mean impact and active peak forces from mean unnormalised data
% Left
[pksL,locsL] = findpeaks(Fz_L_unnormalised, 'MinPeakHeight', BodyMass*9.81, 'SortStr', 'none');


% Right
[pksR,locsR] = findpeaks(Fz_R_unnormalised, 'MinPeakHeight', BodyMass*9.81, 'SortStr', 'none');

% Statistics
mean_impact_peak_L = pksL(1);
mean_impact_peak_R = pksR(1);
mean_active_peak_L = pksL(2);
mean_active_peak_R = pksR(2);
mean_loadingrate_L = mean(loadingrate_L);
mean_loadingrate_R = mean(loadingrate_R);
%% Normalised Fx (Mediolateral)
disp('In Progress - Normalising Fx')
% Left leg
for i=1:length(tdownL)
    x = totaltime(tdownL(i):toffL(i));
    y = [filteredfx_left(tdownL(i):toffL(i))]';
    yy(i,:) = spline(x,[0 y 0]);
    xx(i,:) = linspace(totaltime(tdownL(i)),totaltime(toffL(i)),101);  % create 101 points (0-100%) 
    Fx_L_norm(i,:) = ppval(yy(i,:),xx(i,:)); % rename the normliased Fz_1 to Fz_1_norm
end

% Normalise Fz to Body Weight for graph plotting
Fx_L_BW = Fx_L_norm/BodyMass/9.81;

% Calcualte avearege force-time history of all rearfoot steps
Fx_L_mean = mean(Fx_L_BW,'omitnan');

% Unnormalised force
Fx_L_unnormalised = mean(Fx_L_norm);

% To show plot on its own:
% figure('Name','All steps - Left')
% for i = 2:length(tdownL)
%        plot(Fz_L_BW(i,:),'b'); ylabel('Force (Body Weight)'); xlabel('Ground Contact Time (%)'); hold on
% end

% Right leg
for i=1:length(tdownR)
    x = totaltime(tdownR(i):toffR(i));
    y = [filteredfx_right(tdownR(i):toffR(i))]';
    yy(i,:) = spline(x,[0 y 0]);
    xx(i,:) = linspace(totaltime(tdownR(i)),totaltime(toffR(i)),101);  % create 101 points (0-100%) 
    Fx_R_norm(i,:) = ppval(yy(i,:),xx(i,:)); % rename the normliased Fz_1 to Fz_1_norm
end

% Normalise Fz to Body Weight for graph plotting
Fx_R_BW = Fx_R_norm/BodyMass/9.81;

% Calcualte avearege force-time history of all rearfoot steps
Fx_R_mean = mean(Fx_R_BW,'omitnan');

% Unnormalised force
Fx_R_unnormalised = mean(Fx_R_norm);

% To show plot on its own:
% figure('Name','All steps - Right')
% for i = 2:length(tdownR)
%        plot(Fz_R_BW(i,:),'b'); ylabel('Force (Body Weight)'); xlabel('Ground Contact Time (%)'); hold on
% end
disp('Completed - Normalising Fx')
%% Normalised Fy (Anterior-posterior)
disp('In Progress - Normalising Fy')
% Left leg
for i=1:length(tdownL)
    x = totaltime(tdownL(i):toffL(i));
    y = [filteredfy_left(tdownL(i):toffL(i))]';
    yy(i,:) = spline(x,[0 y 0]);
    xx(i,:) = linspace(totaltime(tdownL(i)),totaltime(toffL(i)),101);  % create 101 points (0-100%) 
    Fy_L_norm(i,:) = ppval(yy(i,:),xx(i,:)); % rename the normliased Fz_1 to Fz_1_norm
end

% Normalise Fz to Body Weight for graph plotting
Fy_L_BW = Fy_L_norm/BodyMass/9.81;

% Calcualte avearege force-time history of all rearfoot steps
Fy_L_mean = mean(Fy_L_BW,'omitnan');

% Unnormalised force
Fy_L_unnormalised = mean(Fy_L_norm);

% To show plot on its own:
% figure('Name','All steps - Left')
% for i = 2:length(tdownL)
%        plot(Fz_L_BW(i,:),'b'); ylabel('Force (Body Weight)'); xlabel('Stance Phase (%)'); hold on
% end

% Right leg
for i=1:length(tdownR)
    x = totaltime(tdownR(i):toffR(i));
    y = [filteredfy_right(tdownR(i):toffR(i))]';
    yy(i,:) = spline(x,[0 y 0]);
    xx(i,:) = linspace(totaltime(tdownR(i)),totaltime(toffR(i)),101);  % create 101 points (0-100%) 
    Fy_R_norm(i,:) = ppval(yy(i,:),xx(i,:)); % rename the normliased Fz_1 to Fz_1_norm
end

% Normalise Fz to Body Weight for graph plotting
Fy_R_BW = Fy_R_norm/BodyMass/9.81;

% Calcualte avearege force-time history of all rearfoot steps
Fy_R_mean = mean(Fy_R_BW,'omitnan');

% Unnormalised force
Fy_R_unnormalised = mean(Fy_R_norm);

% To show plot on its own:
% figure('Name','All steps - Right')
% for i = 2:length(tdownR)
%        plot(Fz_R_BW(i,:),'b'); ylabel('Force (Body Weight)'); xlabel('Stance Phase (%)'); hold on
% end
disp('Completed - Normalising Fy')
%% To compare left and right mean directly in the same plot
disp('In Progress - Plotting normalised graph')
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

% Fy
figure('Name','Normalised Mean GRF')
plot(Fy_L_mean, 'b','LineWidth',2); title('Normalised Mean Anterior-Posterior GRF', 'FontSize', 20);
ylabel('Force (Body Weight)','FontSize',16); xlabel('Stance Phase (%)','FontSize',16);
xlim([0 101]);
hold On
plot(Fy_R_mean, ':r','LineWidth',2);
ylabel('Force (Body Weight)','FontSize',16); xlabel('Stance Phase (%)','FontSize',16);
legend('Left leg', 'Right leg')
hold Off

% Fx
figure('Name','Normalised Mean GRF')
plot(Fx_L_mean, 'b','LineWidth',2); title('Normalised Mean Mediolateral GRF', 'FontSize', 20);
ylabel('Force (Body Weight)','FontSize',16); xlabel('Stance Phase (%)','FontSize',16);
xlim([0 101]);
hold On
plot(Fx_R_mean, ':r','LineWidth',2);
ylabel('Force (Body Weight)','FontSize',16); xlabel('Stance Phase (%)','FontSize',16);
legend('Left leg', 'Right leg')
hold Off
disp('Completed - Normalised graph plotted')
%% Store stance data into array, after removing crossover
disp('In Progress - Storing stance data into array')
% Left leg
for i=1:length(tdownL)
    stanceFzL{i} = filteredfz_left(tdownL(i):toffL(i)); % Store each fz stance data into a cell array
    stanceFxL{i} = filteredfx_left(tdownL(i):toffL(i));
    stanceFyL{i} = filteredfy_left(tdownL(i):toffL(i));
end
% Right leg
for i=1:length(tdownR)
    stanceFzR{i} = filteredfz_right(tdownR(i):toffR(i)); % Store each fz stance data into a cell array
    stanceFxR{i} = filteredfx_right(tdownR(i):toffR(i));
    stanceFyR{i} = filteredfy_right(tdownR(i):toffR(i));
end
disp('Completed - Stored stance data into array')
%% Store mean data into Excel
disp('In Progress - Saving mean stance data into Excel')
% Create Mean excel data is don't have yet
% Append data to next row

a = filename(1:2); % Get participant ID
% b = filename(4:5); % Get session number

% Get from Participant Info if condition is with Exos or without
if filename(5) == '2'
    Exo = info{str2double(filename(2)),5}; % Session 2
elseif filename(5) == '3'
    Exo = info{str2double(filename(2)),6}; % Session 3
else
    Exo = info{str2double(filename(2)),7}; % Session 4
end

c = filename(7:8); % Get trial number

if filename(10) == 'F' % Check if filename is 'Flat'
    d = filename(10:13); % Get full name
elseif filename(10) == 'I' % Check if filaname is 'Incline'
    d = filename(10:16); % Get full name
else % Check if filaname is 'Decline'
    d = filename(10:16); % Get full name
end

e = filename(end-12:end-9); % Get the load in kg

norm_fz_L = string(Fz_L_mean);
norm_fz_R = string(Fz_R_mean);
norm_fx_L = string(Fx_L_mean);
norm_fx_R = string(Fx_R_mean);
norm_fy_L = string(Fy_L_mean);
norm_fy_R = string(Fy_R_mean);
meanfzdata_L = [a Exo d e norm_fz_L];
meanfzdata_R = [a Exo d e norm_fz_R];
meanfxdata_L = [a Exo d e norm_fx_L];
meanfxdata_R = [a Exo d e norm_fx_R];
meanfydata_L = [a Exo d e norm_fy_L];
meanfydata_R = [a Exo d e norm_fy_R];
writematrix(meanfzdata_L, 'Mean_Data_All_Participant.xlsx', 'Sheet', 'Fz_Left', 'WriteMode', 'append')
writematrix(meanfzdata_R, 'Mean_Data_All_Participant.xlsx', 'Sheet', 'Fz_Right', 'WriteMode', 'append')
writematrix(meanfxdata_L, 'Mean_Data_All_Participant.xlsx', 'Sheet', 'Fx_Left', 'WriteMode', 'append')
writematrix(meanfxdata_R, 'Mean_Data_All_Participant.xlsx', 'Sheet', 'Fx_Right', 'WriteMode', 'append')
writematrix(meanfydata_L, 'Mean_Data_All_Participant.xlsx', 'Sheet', 'Fy_Left', 'WriteMode', 'append')
writematrix(meanfydata_R, 'Mean_Data_All_Participant.xlsx', 'Sheet', 'Fy_Right', 'WriteMode', 'append')

% Unnormalised data
unnorm_fz_L = string(Fz_L_unnormalised);
unnorm_fz_R = string(Fz_R_unnormalised);
unnorm_fx_L = string(Fx_L_unnormalised);
unnorm_fx_R = string(Fx_R_unnormalised);
unnorm_fy_L = string(Fy_L_unnormalised);
unnorm_fy_R = string(Fy_R_unnormalised);
unnormfzdata_L = [a Exo d e unnorm_fz_L];
unnormfzdata_R = [a Exo d e unnorm_fz_R];
unnormfxdata_L = [a Exo d e unnorm_fx_L];
unnormfxdata_R = [a Exo d e unnorm_fx_R];
unnormfydata_L = [a Exo d e unnorm_fy_L];
unnormfydata_R = [a Exo d e unnorm_fy_R];
writematrix(unnormfzdata_L, 'Mean_Data_All_Participant_Unnormalised.xlsx', 'Sheet', 'Fz_Left', 'WriteMode', 'append')
writematrix(unnormfzdata_R, 'Mean_Data_All_Participant_Unnormalised.xlsx', 'Sheet', 'Fz_Right', 'WriteMode', 'append')
writematrix(unnormfxdata_L, 'Mean_Data_All_Participant_Unnormalised.xlsx', 'Sheet', 'Fx_Left', 'WriteMode', 'append')
writematrix(unnormfxdata_R, 'Mean_Data_All_Participant_Unnormalised.xlsx', 'Sheet', 'Fx_Right', 'WriteMode', 'append')
writematrix(unnormfydata_L, 'Mean_Data_All_Participant_Unnormalised.xlsx', 'Sheet', 'Fy_Left', 'WriteMode', 'append')
writematrix(unnormfydata_R, 'Mean_Data_All_Participant_Unnormalised.xlsx', 'Sheet', 'Fy_Right', 'WriteMode', 'append')


disp('Completed - Saved mean stance data into Excel')
%% Impulse in each stance (not normalised)
disp('In Progress - Calculating impulse')
% Left
fz_impulse_L = [];
for i=1:length(tdownL)
    fz_impulse_L(end+1) = trapz(filteredfz_left(tdownL(i):toffL(i)))/aRate;
end

% Right
fz_impulse_R = [];
for i=1:length(tdownR)
    fz_impulse_R(end+1) = trapz(filteredfz_right(tdownR(i):toffR(i)))/aRate;
end

% Statistic
mean_impulseL = mean(fz_impulse_L);
mean_impulseR = mean(fz_impulse_R);
disp('Completed - Impulse calculated')

%% Find Fx and Fy max and min in each stance phase
disp('In Progress - Calculating braking and propulsive forces')
Fx_max_L = [];
Fx_min_L = [];
Fy_max_L = [];
Fy_min_L = [];
for i=1:length(tdownL)
    Fx_max_L(end+1) = max(filteredfx_left(tdownL(i):toffL(i)));
    Fx_min_L(end+1) = min(filteredfx_left(tdownL(i):toffL(i)));
    Fy_max_L(end+1) = max(filteredfy_left(tdownL(i):toffL(i)));
    Fy_min_L(end+1) = min(filteredfy_left(tdownL(i):toffL(i)));
end

Fx_max_R = [];
Fx_min_R = [];
Fy_max_R = [];
Fy_min_R = [];
for i=1:length(tdownR)
    Fx_max_R(end+1) = max(filteredfx_right(tdownR(i):toffR(i)));
    Fx_min_R(end+1) = min(filteredfx_right(tdownR(i):toffR(i)));
    Fy_max_R(end+1) = max(filteredfy_right(tdownR(i):toffR(i)));
    Fy_min_R(end+1) = min(filteredfy_right(tdownR(i):toffR(i)));
end

mean_braking_force_L = mean(Fy_min_L);
mean_braking_force_R = mean(Fy_min_R);
mean_propulsive_force_L = mean(Fy_max_L);
mean_propulsive_force_R = mean(Fy_max_R);
disp('Completed - Braking and propulsive forces calculated')
%% Calculate cadence and stride lengths
disp('In Progress - Calculating spatialtemporal variables')
cadence = (stepcountL+stepcountR)/totaltime(end)*60;
step_frequency = cadence/60;
step_length = speed/step_frequency;
stride_length = step_length*2;
disp('Completed - Spatialtemporal variables calculated')
%%  Write summary table for outputs
disp('In Progress - Writing summary table')
Participant_ID = filename(1:2);
Load = filename(end-12:end-9); % To ensure filename has no 'Retry'

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
if strcmp(filename(1:5),'P6_S4')
    T = cell2table(participant, 'VariableNames', ["ID" "Customised_Exo" "Condition" "Load"]);
else
    T = cell2table(participant, 'VariableNames', ["ID" "Exo" "Condition" "Load"]);
end

Stance_Time_L = [mean_stanceL];
Stance_Time_R = [mean_stanceR];
Step_Time_L = [mean_steptimeL];
Step_Time_R = [mean_steptimeR];
Impact_PeakForce_L = [mean_impact_peak_L];
Impact_PeakForce_R = [mean_impact_peak_R];
Active_PeakForce_L = [mean_active_peak_L];
Active_PeakForce_R = [mean_active_peak_R];
Braking_Force_Fy_L = [mean_braking_force_L];
Braking_Force_Fy_R = [mean_braking_force_R];
Propulsive_Force_Fy_L = [mean_propulsive_force_L];
Propulsive_Force_Fy_R = [mean_propulsive_force_R];
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
Braking_Force_Fy_L, Braking_Force_Fy_R, Propulsive_Force_Fy_L, Propulsive_Force_Fy_R, ...
Impulse_L, Impulse_R, Loading_Rate_L, Loading_Rate_R, ...
Cadence, Step_Frequency, Step_Length, Stride_Length);
T1 = [T,T1];
T1
disp('Completed - Summary table')

% Save T1 to Excel without overwriting previous saved data in the same sheet
disp('In Progress - Saving GRF variables into Excel')
writetable(T1, 'Bertec GRF variables_6may_manual.xlsx', 'WriteMode', 'Append')
disp('Completed - Saved GRF variables into Excel')
%% Run EMG data

% If participant ID is 1,2,5,6,8 (with EMG data), run emgdata file
% if filename(2) == '1' || filename(2) == '2' || filename(2) == '5' || filename(2) == '6' || filename(2) == '8'
%     emgdata; % Run emgdata matlab code
% end
disp('Completed - All code run successfully')

% Check if everything is okay, then run savingmeantoexcel.ean