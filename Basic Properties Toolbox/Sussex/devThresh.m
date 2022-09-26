% devThresh
% development version of threshold detector taken from Istep (26.09.22)

% load file of interest
[file,path] = uigetfile('*.ma'); % select file of interest
clampfile = fullfile(path,file); % get full filepath of file
S = ephysIO(clampfile);
Time = S.array(:,1);
Waves = S.array(:,2:end);

% cd to path
splitPath = split(clampfile,filesep);
newPath = char(string(join(splitPath(1:end-2),filesep)));
cd(newPath)

%% run analysis or not?
% gives the user the option to abort if the data looks horrendous

% select region to search for action potentials
detStart = 0.5/S.xdiff; % start of detection in data points
detEnd = 1.5/S.xdiff; % end of detection in data points

% preallocate pks, locs, w, p, numSpikes and mute error here
warning('off','signal:findpeaks:largeMinPeakHeight');
pks = cell(size(Waves,2),1);
locs = cell(size(Waves,2),1);
difflocs = cell(size(Waves,2),1); % difference between loc in dp
normlocs = cell(size(Waves,2),1); % difference between loc in dp
w = cell(size(Waves,2),1);
p = cell(size(Waves,2),1);
numSpikes = zeros(size(Waves,2),1);
nlocs = cell(size(Waves,2),1); % interval index


% findpeaks to determine number and location of AP's per wave
for i = 1:size(Waves,2)
    % pks = value of peak
    % locs = index of peak
    % w = width
    % p = prominence
    % hard coded to only look for AP's between detStart and detEnd
    [pks{i},locs{i},w{i},p{i}] = findpeaks(Waves(round(detStart):round(detEnd),i),...
        'MinPeakHeight',0,...
        'MinPeakProminence',0.01,...
        'MaxPeakWidth',0.02*10^4,...
        'MinPeakDistance',0.02*10^4,...
        'WidthReference','halfheight');
    numSpikes(i) = size(pks{i},1);
    difflocs{i} = diff(locs{i}); % difference between loc in dp
    normlocs{i} = normalize(difflocs{i},'scale','first'); % difflocs norm to first value
    nlocs{i} = 1:size(difflocs{i}); % interval index
end

% First AP of interest values
C = locs;
idx = ~cellfun('isempty',C);
outs = zeros(size(C));
outs(idx) = cellfun(@(v)v(1),C(idx));
for lp = 1:size(outs,1)
    if outs(lp) > 0
        outs(lp) = outs(lp) + detStart; % account for the detection zone
    end
end
logicalIndexes =  outs < 27500 & outs > 1;
wavenum_first = (size(locs,1)-(sum(logicalIndexes)-1)); % wave that the first action potential following stimulus injection occurs on

%% AP analysis - Surpasses 4 mV/ms

% plot action potential
warning('off','MATLAB:colon:nonIntegerIndex')
AP_Window = Waves(outs(wavenum_first)-500:outs(wavenum_first)+500,wavenum_first)*1000;
warning('on','MATLAB:colon:nonIntegerIndex')

% Overshoot in mV
[Overshoot,ind_o] = max(AP_Window);
% Afterhyperpolarisation in mV
[Afterhyperpolarisation,ind_a] = min(AP_Window(ind_o:end));
% Baseline
Base = mean(AP_Window(1:350));

% Action potential halfwidth
% Halfwidth in ms
figure;
% subplot(7,4,[4,8,12]);
[~,~,Halfwidth,~] = findpeaks(AP_Window-Base,'MinPeakHeight',0,...
        'MaxPeakWidth',0.02*10^4,...
        'MinPeakDistance',600,...
        'WidthReference','halfheight',...
        'Annotate','extent');
findpeaks(AP_Window-Base,'MinPeakHeight',0,...
        'MaxPeakWidth',0.02*10^4,...
        'MinPeakDistance',600,...
        'WidthReference','halfheight',...
        'Annotate','extent');
box off; grid off; xlabel('Data Points'); ylabel('Adjusted membrane potential');
set(gca,'linewidth',2); set(gcf,'color','white'); title('4 mV/ms');
ylim([-40 120])
Halfwidth = (Halfwidth*Time(2))*1000; % Halfwidth in ms

% Action Potential Threshold

% dv/dt surpasses 4mV/ms (dv/dt > 1)
N = diff(AP_Window(1:ind_o-20)); % period start to just before the peak
[closestValue, closestIndex] = min(abs(N - 1.'));
hold on; plot(closestIndex,AP_Window(closestIndex)-Base,'or')
ind_t = closestIndex; % threshold index=
Threshold = AP_Window(ind_t); % in mV

hold on; plot((ind_a+ind_o),Afterhyperpolarisation-Base,'* y') % plot after

% recalculate (and overwrite) the amplitude now you have a threshold
Amplitude = abs((Overshoot-Base)-(Threshold-Base)); % in mV

% Depolarisation Rate
% identify the closest value (and index) to the threshold value after peak
N = AP_Window(ind_t:ind_o)-Base; % period thresh - peak
[~, closestIndex_20] = min(abs(N - 0.2*Amplitude.'));
[~, closestIndex_80] = min(abs(N - 0.8*Amplitude.'));

% between 20 % and 80 % of the rise phase (based on amplitude, not index)
rise_20_ind = closestIndex_20 + ind_t;
rise_80_ind = closestIndex_80 + ind_t;
Rise = mean(gradient(AP_Window(rise_20_ind:rise_80_ind)));
hold on; plot((rise_20_ind:rise_80_ind),AP_Window(rise_20_ind:rise_80_ind)-Base,'linewidth',2)

% Repolarisation Rate

% identify the closest value (and index) to the threshold value after peak
N = AP_Window(ind_o:ind_o+ind_a)-Base; % period peak - hyper
[~, closestIndex_20] = min(abs(N - 0.2*Amplitude.'));
[~, closestIndex_80] = min(abs(N - 0.8*Amplitude.'));

% between 20 % and 80 % of the rise phase (based on amplitude, not index)
fall_20_ind = closestIndex_80 + ind_o;
fall_80_ind = closestIndex_20 + ind_o;
Fall = mean(gradient(AP_Window(fall_20_ind:fall_80_ind)));
hold on; plot((fall_20_ind:fall_80_ind),AP_Window(fall_20_ind:fall_80_ind)-Base,'linewidth',2)

legend('trace','peak', ...
    'amplitude','halfwidth', ...
    'border','threshold', ...
    'dv/dt',...
    'hyperpolar.',...
    'rise','fall',...
    'Location','northeast')
%% AP analysis - Surpasses 4 mV/ms


% Overshoot in mV
[Overshoot,ind_o] = max(AP_Window);
% Afterhyperpolarisation in mV
[Afterhyperpolarisation,ind_a] = min(AP_Window(ind_o:end));
% Baseline
Base = mean(AP_Window(1:350));

% Action potential halfwidth
% Halfwidth in ms
figure;
% subplot(7,4,[4,8,12]);
[~,~,Halfwidth,~] = findpeaks(AP_Window-Base,'MinPeakHeight',0,...
        'MaxPeakWidth',0.02*10^4,...
        'MinPeakDistance',600,...
        'WidthReference','halfheight',...
        'Annotate','extent');
findpeaks(AP_Window-Base,'MinPeakHeight',0,...
        'MaxPeakWidth',0.02*10^4,...
        'MinPeakDistance',600,...
        'WidthReference','halfheight',...
        'Annotate','extent');
box off; grid off; xlabel('Data Points'); ylabel('Adjusted membrane potential');
set(gca,'linewidth',2); set(gcf,'color','white'); title('triple ndiff peak');
ylim([-40 120])
Halfwidth = (Halfwidth*Time(2))*1000; % Halfwidth in ms

% Action Potential Threshold

% triple ndiff
diffWin = AP_Window;
for diffnum = 1:3
    diffWin = ndiff(diffWin,Time(1:size(diffWin,1)));
end

% plot the differntial below the trace
hold on; plot((diffWin)/100-20);
[val,ind] = max(diffWin); xline(ind)

ind_t = closestIndex; % threshold index
% plot the differntial below the trace
Threshold = AP_Window(ind_t); % in mV

hold on; plot((ind_a+ind_o),Afterhyperpolarisation-Base,'* y') % plot after

% recalculate (and overwrite) the amplitude now you have a threshold
Amplitude = abs((Overshoot-Base)-(Threshold-Base)); % in mV

% Depolarisation Rate
% identify the closest value (and index) to the threshold value after peak
N = AP_Window(ind_t:ind_o)-Base; % period thresh - peak
[~, closestIndex_20] = min(abs(N - 0.2*Amplitude.'));
[~, closestIndex_80] = min(abs(N - 0.8*Amplitude.'));

% between 20 % and 80 % of the rise phase (based on amplitude, not index)
rise_20_ind = closestIndex_20 + ind_t;
rise_80_ind = closestIndex_80 + ind_t;
Rise = mean(gradient(AP_Window(rise_20_ind:rise_80_ind)));
hold on; plot((rise_20_ind:rise_80_ind),AP_Window(rise_20_ind:rise_80_ind)-Base,'linewidth',2)

% Repolarisation Rate

% identify the closest value (and index) to the threshold value after peak
N = AP_Window(ind_o:ind_o+ind_a)-Base; % period peak - hyper
[~, closestIndex_20] = min(abs(N - 0.2*Amplitude.'));
[~, closestIndex_80] = min(abs(N - 0.8*Amplitude.'));

% between 20 % and 80 % of the rise phase (based on amplitude, not index)
fall_20_ind = closestIndex_80 + ind_o;
fall_80_ind = closestIndex_20 + ind_o;
Fall = mean(gradient(AP_Window(fall_20_ind:fall_80_ind)));
hold on; plot((fall_20_ind:fall_80_ind),AP_Window(fall_20_ind:fall_80_ind)-Base,'linewidth',2)

legend('trace','peak', ...
    'amplitude','halfwidth', ...
    'border','threshold', ...
    'dv/dt',...
    'hyperpolar.',...
    'rise','fall',...
    'Location','northeast')
