function [output] = IStep(LPF_Hz, winSize_ms,adaptI)
%% Decription
% Current step protocol analysis, featuring user defined input at various
% points. Due to built in ephysIO functionality, the user only needs select
% the recording file for the first wave in the set of recordings and
% ephysIO will automatically load the remaining files taking the order from
% the folder name (ie, '000', '001, '002', etc).
%
% IStep v1.3.2 (last updated: 03/11/23)
% Author: OGSteele
%
% example use;
%       output = IStep(LPF_Hz, winSize_ms);
%       eg. cellOne = IStep(330,25,300)
%
% output 
%   is a matlab structure containing the following;
%       filepath - filepath of the first wave
%       ID - slice_ID on plate
%       steps - current amplitude of each step in pA
%       waveform - current step protocol waveform
%       time - time in s
%       waves - waves in V
%       ephysIO - ephysIO output format, with associated metadata if needed
%       numSpikes - number of spikes per wave
%       Rh - Rheobase in pA
%       sag_mV - Ih Sag Amplitude in mV
%       sag_ratio - Ih Sag : steady state ratio
%       peak - Overshoot in mV (note, not same as figure)
%       afterhyp - Afterhyperpolarisation value in mV
%       amp - Action potential amplitude in mV
%       thresh - Threshold potential in mV
%       half - Halfwidth in ms
%       rise - Depolarisation rate in mV/s
%       fall - Repolarisation rate in mV/s
%       IR - Input Resistance in MOhm
%       Vm - resting membrane potential
%       Vm_stability - does Vm fluctuate more than 3stdev (caution if so)
%       subAP_Vm - sub action potential Vm in mV following stimulation
%       Rs_Init - Access resistance approximated from initial step  in MOhm
%       peakTimes - time of peaks in s
%       peakAmps - overshoot of peaks in mV
%       ISI - interspike interval in s
%       peakHz - instantaneous spike frequency in Hz
%       adaptIndex - adaptation index
%       Offline_BB - Vm adjustments for offline bridge balance in V
%       Online_BB_performed = Yes/No 
%       Offline_BB_performed = Yes/No 
%       notes - notes input at the end
%       LPF_Hz - Low pass filter cut off in Hz applied to AP wave
%       winSize_ms - Window size in ms for detection of action potentials
%
% inputs
%       LPF_Hz is the low pass filter cut off in Hz (only ndiff filtered)
%       winSize_ms is the AP window size in ms (peak index +/- winsize/2)
%       adaptI is the current step value (pA) to analyse adaptation at
%
% -----
% Note on Inputs
% Generally, the faster the AP, the lower the LPF_Hz cutoff required
%   tested on (approximate appropriate filter cutoff); 
%       - organotypic DIV21 CA1 pyramidal neurons @ RT (330 Hz)
%       - cultured iPSC derived neurons @ RT (1000 Hz)
% Inverse is true for window size to capture the afterhyperpolarisation;
%   tested on (approximate appropriate window size); 
%       - organotypic DIV21 CA1 pyramidal neurons @ RT (25 ms)
%       - cultured iPSC derived neurons @ RT (50 ms)
% Advised to use devThresh.m script to determine optimal values
%
% -----
% Dependancies
%   - Signal Processing Toolbox (mathworks)
%
% -----
% Notes on paths
% Users should avoid having data stored on the path as this can confuse
% ephysIO. Users should avoid storing data on the path if that data is 
% stored in sequential folders ('001', '002' etc) that has a higher folder 
% number than the data you're telling it to analyse (ie, the data on the 
% path has 31 folders and the data you're interested in analysing has 30 
% folders) as this will cause ephysIO, and therefore IStep, to 
% throw an error ('cannot change to folder xxx as it is non-existant'
% and terminate the function prematurely.

%% Update Log

% 17.09.22 [OGS] v1.1
%   - improve cross platform functionality (filesep throughout)
%   - move saved figure and output to the root folder of the recording
%   - correct file naming bug
%   - include subAP_Vm figure too (renamed old fig to master fig)
%   - fixed lazy rheobase annotation alignment
%   - read in command potential waveform data from .ma recordings
%   - allowed UI selection of where to detect action potentials from
%   - increased robustness of the Ih Sag plotting
%   - included the use of ephysIO inside the script for ease of use
%   - tidied up description

% 04.10.22 [OGS] v1.2.1
%   - implemented robust threshold detection (initial peak of 2div)
%   - tidied up devThresh and IStep as a package
%   - introduced LPF_Hz and winSize_ms input arguments
%   - updated IStep description
%   - tidied up labelling on subgraphs 3 (Rh) and 7 (IR)
%   - updated baseline selection to be a user input
%   - plotted detection regions on the overall traces (subgraph 1)
%   - corrected AP Analysis Labelling
%   - adjusted minPeakHeight in AP Analysis to prevent double spikes
%   - included 7 pole medianf to decrease impact of fast noise in Ih Sag

% 20.10.22 [OGS] v1.2.2
%   - updated dependancies to include readMeta.m and the signal processign
%   toolbox
%   - all figures closed on running

% 21.10.22 [OGS] v1.2.3
%   - clearer description of ephysIO functionality and multiple succesive
%   file loading (also put in display line about this)
%   - included dependancies list in devThresh for that to work
%   independantly also
%   - introduced catch me if user doesn't select input arguments

% 09.11.22 [OGS] v1.2.4
%   - transparent labelling of rise and fall times, however benefit of this
%   is limited due it being an overlay of the actual trace that is thicker.
%   Consider amending this in future. 
%   - corrected IR bug where offline bridge balance correction was
%   performed before the output structure was created, but after the figure
%   was plotted. Shouldn't change the analysis, but confusing none the
%   less. 
%   - optional save loop at the end

% 15.11.22 [OGS] v1.2.5
%   - inclusion of waitbar to tell users not to close figures during saving

% 07.12.22 [OGS] v1.2.6
%   - clarification of point relating to peak in output and peak in figure.
%   The figure value is not real, as it's adjusted to place the baseline at
%   0 mV. The overall amplitude value is therefore not affected. The peak
%   found in the output corresponds to the non-baseline adjusted value, ie
%   the true most positive value. This value is adjusted during offline
%   bridge balancing. 
%   - Note: Offline bridge balancing not always appropriate when the offset
%   is positive rather than negative. 
%   - Also noted strange occurance of pA corruption in the waveform on
%   certain recordings. Appears inconsistent, may be affected by OneDrive?
%   Changes the order of magnitude of the pA waveform (-2e-10 to
%   -4e-9). Will continue to investigate. 

% 12.07.23 [ACP] v1.2.7
%   - fixed bug in Ih sag calculation that estimated the Ih sag from
%   assumed Vm of -65 mV rather than the measured value
%   - calculation of Ih Sag ratio as 1 minus this value

% 20.09.23 [OGS] v1.2.8
%   - amended description to highlight known issue with path setting that
%   relates back to ephysIO, an expected behaviour. 
%   - included reporting of Vm as well as flagging potential issue with
%   seal quality

% 27.09.23 [OGS] v1.2.9
%   - Corrected Figure 3 to be the correct number of SDs on the legend
%   - Ensure Figure 3 is now saving, and corrected the relative order of
%   saving. 
%   - Removed todolist, see Github issues
%   - added hard toggle for saving of figures to save on time taken and
%   storage space required (on by default)
%   - stopped code from crashing if cancelled file selection

% 05.10.23 [OGS] v1.3
%   - Inclusion of spike frequency adaptation at a given stimulus
%   intesnsity

% 1.11.23 [OGS] v1.3.1
%   - Increase window to 1s following current injection to catch delayed
%   charging of membrane. Not sure this is long term fix, and will look to
%   include plotting of this region
%   - Made the calculation of IR much more robust with averaging of three
%   different IR's before the onset of Ih current contamination. Will need
%   to consider display of this information.
%   - Made searching for the zerowave more robust, as it now finds the
%   absolute minimum rather than hard coding zero.

% 03.11.23 [OGS] v1.3.2
%   - minor improvement to robustness, removing the hardcoding on 'wave 11'
%   as 'zerowave'
%   - removed the zerowave Ih sag as this made no sense and made it hard to
%   visualise the the average

%% Code

% figure save toggle
figsave = 1; % 1 = true, 0 = false

% close figures
close all

% determine the number of input arguments
numInputs = nargin;

% load file of interest
disp('---------')
disp('Select the first file (from the first wave) in the recording, ephysIO will load the rest automatically')
disp('---------')
[file,path] = uigetfile('*.ma'); % select file of interest
clampfile = fullfile(path,file); % get full filepath of file
if size(clampfile,2) < 5 % if the filepath is too short to make any sense
    disp('File selection cancelled, code aborted')
    return
else % if filepath is longer than 5 characters, file selected
    S = ephysIO(clampfile);
    Time = S.array(:,1);
    Waves = S.array(:,2:end);
    
    % cd to path
    splitPath = split(clampfile,filesep);
    newPath = char(string(join(splitPath(1:end-2),filesep)));
    cd(newPath)
    
    %% run analysis or not?
    % gives the user the option to abort if the data looks horrendous
    
    % plot overall figure
    fh = figure();
    plot(Time,Waves(:,1)*1000,'color','black')
    hold on; plot(Time,Waves*1000,'color','black','HandleVisibility','off')
    box off; set(gcf,'color','white'); set(gca,'linewidth',2)
    ylabel('Membrane Potential (mV)'); xlabel('Time (s)')
    title('Current Step Waveform')
    
    % ask the user 
    dlgTitle    = 'Run Analysis';
    dlgQuestion = 'Would you like to run the Current Step Analysis?';
    run = questdlg(dlgQuestion,dlgTitle,'Yes','No','No');
    
    if run == "Yes"
        disp('---------')
        disp('Performing Analysis, please wait ...')
        disp('---------')
    
        % select region to search for action potentials
        title('Select detection region')
        [detection,~] = ginput(2);
        detStart = round(detection(1)/S.xdiff); % start of detection in data points
        detEnd = round(detection(2)/S.xdiff); % end of detection in data points
        xline(detection(1),'linestyle','--','color','blue','linewidth',2,'HandleVisibility','off')
        xline(detection(2),'linestyle','--','color','blue','linewidth',2)
    
        % select region to detect the baseline from
        title('Select baseline region')
        [baseline,~] = ginput(2);
        baseStart = round(baseline(1)/S.xdiff); % start of detection in data points
        baseEnd = round(baseline(2)/S.xdiff); % end of detection in data points
        xline(baseline(1),'linestyle','--','color','green','linewidth',2,'HandleVisibility','off')
        xline(baseline(2),'linestyle','--','color','green','linewidth',2)
        legend('Recording Data','AP Detection Region','Baseline Detection Region','linewidth',1)
        title('Detection Regions')
        Vbase = mean(Waves(baseStart:baseEnd,:)); % determine the baseline (intra-step)
    
        % pause for 2 seconds to allow user to visualise regions, then close
        pause(2)
        % close figure
        close(fh)
    
    
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
        
        % plot whole of current step protocol
        fh = figure();
        fh.WindowState = 'maximized'; subplot(7,4,[1,5]); plot(Time,Waves*1000,'color','black','HandleVisibility','off')
        box off; set(gcf,'color','white'); set(gca,'linewidth',2)
        ylabel('Membrane Potential (mV)');
        title('Current Step Waveform')
        ax = gca; xax = ax.XAxis; set(xax,'visible','off')
        hold on
        xline(detection(1),'linestyle','--','color','blue','linewidth',2,'HandleVisibility','off') % detection start
        xline(detection(2),'linestyle','--','color','blue','linewidth',2) % detection end
        xline(baseline(1),'linestyle','--','color','green','linewidth',2,'HandleVisibility','off') % baseline start
        xline(baseline(2),'linestyle','--','color','green','linewidth',2) % baseline end
        hold off
        legend('Detection Region','Baseline Region','linewidth',1)
        
        % plot waveform of current step read from second channel of clampfit
        SI = ephysIO({clampfile,3}); % load in the current data
        x = SI.array(:,1); % x is time here
        pA_waveform = SI.array(:,2:end); % y is the array of current data
        Ibase = mean (mean(pA_waveform(baseStart:baseEnd,:))); % determine the baseline (intra-step)
        Vbase = mean(Waves(baseStart:baseEnd,:)); % resting membrane potential in Vm of each step
        N = size(SI.array,2) - 1; % get the number of waves in the array (minus time)
        lo = dsearchn(x,detection(1)) + 1; % start of the test pulse
        hi = dsearchn(x,detection(2)) - 1; % end of the test pulse
        pA = zeros(1,N); % preallocate a blank series of steps
        for i = 1:N
           pA(i) = mean (pA_waveform(round(lo+(detEnd-detStart)*0.1):...
               round(hi-(detEnd-detStart)*0.1),i)); % fill the steps with the mean of each step
        end
        pA = fix((pA - Ibase) * 1e+12); % round to zero, baseline subtract and put into pA
        subplot(7,4,9); plot(x,(pA_waveform-Ibase)*1e12,'linewidth',1,'color','black','HandleVisibility','off')
        box off; set(gca,'linewidth',2); set(gcf,'color','white');
        xlabel('Time (s)'); ylabel('Command(pA)'); ylim([min(min(pA_waveform*1e12))-50,max(max(pA_waveform*1e12))+50])
        ax = gca; yax = ax.YAxis; set(yax,'TickDirection','out')
        hold on
        xline(detection(1),'linestyle','--','color','blue','linewidth',2,'HandleVisibility','off') % detection start
        xline(detection(2),'linestyle','--','color','blue','linewidth',2) % detection end
        xline(baseline(1),'linestyle','--','color','green','linewidth',2,'HandleVisibility','off') % baseline start
        xline(baseline(2),'linestyle','--','color','green','linewidth',2) % baseline end
        hold off
        
        % First AP of interest values
        [v,zerowave] = min(abs(pA)); % find wave closest to zero
        C = locs;
        idx = ~cellfun('isempty',C);
        outs = zeros(size(C));
        outs(idx) = cellfun(@(v)v(1),C(idx));
        for lp = 1:size(outs,1)
            if outs(lp) > 0
                outs(lp) = outs(lp) + detStart; % account for the detection zone
            end
        end
        logicalIndexes =  outs < 40000 & outs > 1;
        wavenum_first = (size(locs,1)-(sum(logicalIndexes)-1)); % wave that the first action potential following stimulus injection occurs on
        subplot(7,4,[2,6,10]);
        plot(Time,Waves(:,wavenum_first)*1000,'color','red'); hold on; plot(Time,Waves(:,zerowave)*1000,'color','black')
        box off; set(gcf,'color','white'); set(gca,'linewidth',2)
        xlabel('Time (s)'); ylabel('Membrane Potential (mV)');
        lgd = legend(char(string(pA(wavenum_first))),char(string(pA(zerowave))),'linewidth',1);
        title(lgd,'Current (pA)')
        title('Exemplary Waves')
        
        %% Membrane Excitability and Rheobase
        % determine rheobase as the first amount of current to induce APs in the
        % first 25% of the current step rather than the first current step value to
        % elcit any AP at all. 
        
        Rh = pA(wavenum_first); % Rheobase in pA
        
        % plot membrane excitability and Rheobase
        subplot(7,4,[3,7,11])
        plot(pA,numSpikes,'color','red','linewidth',3); box off; set(gcf,'color','white'); set(gca,'linewidth',2);
        title('Membrane Excitability'); xlabel('Current Step (pA)'); ylabel('Number of Action Potentials')
        hold on; xline(Rh,'--','linewidth',1.5); 
        txt_Rh_1 = {['\bf Rheobase:  '] [num2str(Rh) ' pA  \rightarrow']};
        %text(Rh-135,8,txt_Rh_1);
        %h = annotation('textbox', [0 1 0 0], 'String', 'YourString', 'FitBoxToText', true);
        legend('Number of APs','Rheobase','linewidth',1,'Location','northwest');
        
        %% Ih Sag Values
        % measured after a depolarising current injection of -100 pA from -65 mV
        % Calculates both Amplitude (relative to steady state) and ratio (to steady state)
        
        % median filter Ih Sag values to remove any noise
        for i = 1:size(Waves,2)
            mf_Waves(:,i) = medianf(Waves(:,i),Time,9);
        end
    
        % Plot Ih Sag Waves
        subplot(7,4,[17,21,25]); % creat subplot
        steadyStateWaveNum = zerowave-1; % find the wave number with the closest to zerohold on; 
        plot(Time,mf_Waves(:,2:steadyStateWaveNum-1)*1000, 'color',[0.8,0.8,0.8],'HandleVisibility','off') % plot the other of the waves in gray in the background
        hold on; plot(Time,mf_Waves(:,steadyStateWaveNum)*1000,'color','black'); % plot steady state waveth wave (ie, zero input) 
        hold on; plot(Time,mf_Waves(:,1)*1000,'color','red') % plot first wave in red state waveth wave (ie, zero input) 
        box off; set(gcf,'color','white'); set(gca,'linewidth',2)
        xlabel('Time (s)'); ylabel('Membrane Potential (mV)');
        lgd = legend(char(string(pA(zerowave))),char(string(pA(1))),...
            'linewidth',1,...
            'location','northwest',...
            'AutoUpdate','off');
        title(lgd,'Current (pA)')
        title('Ih Sag Calculation')
        
        % Calculate Ih Sag 
        % steady state during the final quarter and sag during the first third
        detDur = detEnd-detStart; % calculate detection duration
        SS_start = round(detEnd - (detDur/4)); % start of the steady state zone
        SS_end = detEnd; % end of steady state zone
        sagStart = detStart; % start of sag region
        sagEnd = round(detStart + (detDur/3)); % end of sag region
        % preallocate for loop variables
        % plotting should be done to where the current input is zero
        y = zeros(1,steadyStateWaveNum);
        x = zeros(1,steadyStateWaveNum);
        SS_Value = zeros(1,steadyStateWaveNum);
        Ih_Sag_Amp = zeros(1,steadyStateWaveNum);
        Ih_Sag_Percentage = zeros(1,steadyStateWaveNum);
        for i = 1:steadyStateWaveNum
            [y(i),x(i)] = min(mf_Waves(sagStart:sagEnd,i));
            SS_Value(i) = mean(mf_Waves(SS_start:SS_end,i));
            Ih_Sag_Amp(i) = (SS_Value(i)-y(i))*1000; % Sag amplitude in mV
            Ih_Sag_Percentage(i) = (1-(SS_Value(i)-Vbase(i))/(y(i)-Vbase(i)))*100; % Sag percentage
        end 
        hold on; yline((SS_Value(1))*1000,'--r'); yline((y(1)*1000),'--r');
        ylim([(min(y)*1000) - 20 , max(max(Waves(:,1:steadyStateWaveNum))*1000) + 20]); % make sure graph fits neatly
        
        subplot(7,4,[18,22,26]);
        plot(pA(1:steadyStateWaveNum),Ih_Sag_Percentage,'color','black','linewidth',3); box off; title('Ih Sag')
        set(gcf,'color','white'); xlabel('Current Step (pA)'); ylabel('Ih Sag - Steady State Ratio (%)');
        ylim([0,round(max(Ih_Sag_Percentage)*1.2)]); yline(trimmean(Ih_Sag_Percentage,10),'--b')
        legend('Ih sag percentage', 'Combined average','linewidth',1)
        set(gca,'linewidth',2)
        
        %% AP analysis
        
        % catch me if user doesn't select input arguments at the start
        if numInputs < 2
            % prompt the user to select the two inputs here
            disp('enter input arguments into dialog box here')
            prompt = {...
                'Enter Low Pass Filter Cutoff (Hz):',...
                'Enter Window Size (ms):'};
            dlg_title = 'Input parameters';
            num_lines = 1;
            def = {'330','25'};
            answer  = inputdlg(prompt,dlg_title,num_lines,def,'on');
            answer = str2double(answer);
            LPF_Hz = answer(1);
            winSize_ms = answer(2);
        end
    
        % determine window size
        % 25 is good for organotypic
        % 50 is good for slower APs
        %winSize_ms = inputdlg('Choose AP window size (ms)','winSize_ms',1,{'longer window for slow APs'});
        %winSize_ms = str2num(cell2mat(winSize_ms)); % convert to number
        winSize_ms = (winSize_ms*1e-3)/S.xdiff; % (ms to s)to data points
        
        
        % plot action potential
        warning('off','MATLAB:colon:nonIntegerIndex')
        AP_Window = Waves(outs(wavenum_first)-(winSize_ms/2):outs(wavenum_first)+(winSize_ms/2),wavenum_first)*1000;
        warning('on','MATLAB:colon:nonIntegerIndex')
        
        
        % Overshoot in mV
        [Overshoot,ind_o] = max(AP_Window);
        % Afterhyperpolarisation in mV
        [Afterhyperpolarisation,ind_a] = min(AP_Window(ind_o:end));
        % Baseline
        Base = mean(AP_Window(1:350));
        
        % Action potential halfwidth
        % Halfwidth in ms
        subplot(7,4,[4,8,12]);
        [~,~,Halfwidth,~] = findpeaks(AP_Window-Base,'MinPeakHeight',15,...
                'MaxPeakWidth',0.02*10^4,...
                'MinPeakDistance',600,...
                'WidthReference','halfheight',...
                'Annotate','extent');
        findpeaks(AP_Window-Base,'MinPeakHeight',15,...
                'MaxPeakWidth',0.02*10^4,...
                'MinPeakDistance',600,...
                'WidthReference','halfheight',...
                'Annotate','extent');
        box off; grid off; xlabel('Data Points'); ylabel('Adjusted membrane potential');
        set(gca,'linewidth',2); set(gcf,'color','white'); title('Action Potential Analysis');
        ylim([-40 120])
        xlim([300 1000]) % Fix for Kate
        Halfwidth = (Halfwidth*Time(2))*1000; % Halfwidth in ms
        
        % Action Potential Threshold
        
        % triple ndiff
        diffWin = AP_Window;
        for diffnum = 1:2
            diffWin = ndiff(diffWin,Time(1:size(diffWin,1)));
        end
        
        % LPF_Hz Note <-- faster the action potential, lower the threshold
        % 330 works well for organotypic, 1000 works well for iPSCs
        
        % plot the differntial below the trace
        y = filter1(diffWin,Time(1:size(diffWin,1)),0,LPF_Hz)*5e-8-20; 
        hold on; plot(y); % plotting of the actual triple diff trace
        
        % find peaks in the double diff trace
        [~,dd_locs,~,~] = findpeaks(y,"MinPeakProminence",0.5); 
        % make sure to understand peak prominence properly
        
        % Discover the intial peak of the ndiff^2 trace
        ind_t = dd_locs(1); % threshold index is the first peak detected
        Threshold = AP_Window(ind_t); % in mV <-- to go to the output
        hold on; plot(ind_t,y(ind_t),'*b','LineWidth',3) % plot the initial peak
        hold on; plot(ind_t,AP_Window(ind_t)-Base,'or','LineWidth',3) % plot the threshold
        
        % plot the afterhyperpolarisation
        hold on; plot((ind_a+ind_o),Afterhyperpolarisation-Base,'* y','LineWidth',3) % plot after
        
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
        hold on; plot((rise_20_ind:rise_80_ind),AP_Window(rise_20_ind:rise_80_ind)-Base,'linewidth',3,'color',[0,0,0,0.5])
        
        % Repolarisation Rate
        
        % identify the closest value (and index) to the threshold value after peak
        
        N = AP_Window(ind_o:ind_o+ind_a)-Base; % period peak - hyper
        [~, closestIndex_20] = min(abs(N - 0.2*Amplitude.'));
        [~, closestIndex_80] = min(abs(N - 0.8*Amplitude.'));
        
        % between 20 % and 80 % of the rise phase (based on amplitude, not index)
        fall_20_ind = closestIndex_80 + ind_o;
        fall_80_ind = closestIndex_20 + ind_o;
        Fall = mean(gradient(AP_Window(fall_20_ind:fall_80_ind)));
        hold on; plot((fall_20_ind:fall_80_ind),AP_Window(fall_20_ind:fall_80_ind)-Base,'linewidth',3,'color',[1,0,0,0.5])
        
        legend('trace','peak', ...
            'amplitude','halfwidth', ...
            'double ndiff','initial max','threshold',...
            'hyperpolar.',...
            'rise','fall',...
            'Location','northeast')
        
        % plot action potential waveform
        subplot(7,4,[20,24,28]); plot(AP_Window, gradient(AP_Window),'linewidth',2,'color','black')
        box off; title('Action Potential Waveform')
        set(gca,'linewidth',2); set(gcf,'color','white'); xlabel('Membrane Potential (mV)'); ylabel('dV/dt');
        
        %% Input Resistance
        % Calculates input resistance in MOhm from the difference in voltage, divided by
        % the current applied, to the the steady state potentials on the last two
        % waves
        
        % plot waves used for input resistance calculation
        subplot(7,4,[19,23,27]);

        plot(Time,Waves(:,zerowave)*1000,'color','red'); 
        hold on; plot(Time,Waves(:,zerowave-3:zerowave-1)*1000,'color','black')
        box off; set(gcf,'color','white'); set(gca,'linewidth',2)
        xlabel('Time (s)'); ylabel('Membrane Potential (mV)');
        lgd = legend(char(string(pA(zerowave))), ...
            char(string(pA(zerowave-1))),...
            char(string(pA(zerowave-2))),...
            char(string(pA(zerowave-3))),...
            'linewidth',1,...
            'location','southeast',...
            'AutoUpdate','off');
        title(lgd,'Current (pA)')
        title('Input Resistance Calculation')
        %show lines for region of IR determination
        hold on; xline(1.5,'--'); xline(1.3,'--')
        
        % Calculate Input Resistance in MegaOhms
        IR_start = round(detEnd - (detDur/4)); % start of the steady state zone
        IR_end = detEnd; % end of steady state zone

        for i = 1:3
            deltaV(i) = abs(mean(Waves(IR_start:IR_end,zerowave)) - mean(Waves(IR_start:IR_end,zerowave-i))); % Delta_Voltage (Volts)
            I(i) = (pA(zerowave-i)-pA(zerowave))*1e-12; % I (Amps)
            R(i) = deltaV(i) / I(i); % R (Ohms)
            IR(i) = R(i) / 1e6; % R (MegaOhms)
        end
        aveIR = mean(IR);

        deltaV = abs(mean(Waves(IR_start:IR_end,1)) - mean(Waves(IR_start:IR_end,3))); % Delta_Voltage (Volts)
        I = (pA(3)-pA(1))*1e-12; % I (Amps)
        R = deltaV / I; % R (Ohms)
        IR = R / 1e6; % R (MegaOhms)
        % temporarily hiding this until sorting
        %txt = {['\bf Input Resistance: '],[num2str(aveIR) ' M\Omega']};
        %hold on; annotation('textbox',[.2 .5 .3 .3],'String',txt,'FitBoxToText','on');
        
        %% Sub-AP Vm values
        Vm_start = round(detEnd - (detDur/2)); % start of the steady state zone
        Vm_end = detEnd; % end of steady state zone
        
        C = locs;
        idx = ~cellfun('isempty',C);
        outs = zeros(size(C));
        outs(idx) = cellfun(@(v)v(1),C(idx));
        logicalIndexes =  outs > 1;
        wavenum_first = (size(locs,1)-(sum(logicalIndexes)-1)); % wave that the first action potential following stimulus injection occurs on
        subAP_Vm = mean(Waves(Vm_start:Vm_end,1:wavenum_first-1));
        
        
        % plot the subAP_Vm values
        fh2 = figure; subplot(1,2,1)
        plot(Time,Waves(:,1),'color','black'); hold on
        plot(Time,Waves(:,2:wavenum_first-1),'color','black','HandleVisibility','off')
        xline(1,'--r'); xline(1.5,'--r','HandleVisibility','off')
        plot(1.25,subAP_Vm(1),'ob')
        for i = 2:size(subAP_Vm,2)
        hold on; plot(1.25,subAP_Vm(i),'ob','HandleVisibility','off'); hold off
        end
        box off; set(gcf,'color','white'); set(gca,'linewidth',2)
        legend('Waves','Averaged Period','Average Vm','linewidth',1,'autoupdate','off','location','southeast')
        xlabel('Time (s)'); ylabel('Membrane Potential (mV)');
        
        subplot(1,2,2); plot(pA(1:wavenum_first-1),subAP_Vm,'-o','color','blue','linewidth',3)
        box off; set(gcf,'color','white'); set(gca,'linewidth',2)
        xlabel('Current Step (pA)'); ylabel('Membrane Potential (mV)');
    
        %% Vm reporting
        Vbase = Vbase*1000;  % convert to mV
        fh3 = figure; plot(Vbase,'-ob','linewidth', 2)
        hold on; yline(mean(Vbase)+3*std(Vbase),'--r')
        hold on; yline(mean(Vbase)-3*std(Vbase),'--r', ...
            'HandleVisibility','off')
        box off; set(gcf,'color','white'); set(gca,'linewidth',2)
        legend('Vm','+/- 3 SD','linewidth',1,'autoupdate','off','location','northeast')
        xlabel('I Step (#)'); ylabel('Membrane Potential (mV)');
    
        % save average Vm for output
        Vm = trimmean(Vbase,10);
    
        % flagging if out of range
        stable = any(~inrange(Vbase,[(mean(Vbase) - 3*std(Vbase)) (mean(Vbase) + 3*std(Vbase)) ]));
        if stable == 0
        title('Suggested stability check passed')
        Vm_stability = "stable";
        else
        title('Suggested stability check failed')
        f = msgbox("Caution: Possible unstable patch","Stability Issue","error");
        Vm_stability = "caution";
        end

        %% Spike Frequency Adaptation
        % plot whole of current step protocol
        fh4 = figure();
        fh4.WindowState = 'maximized'; subplot(7,3,[1,4]); plot(Time,Waves*1000,'color','black','HandleVisibility','off')
        box off; set(gcf,'color','white'); set(gca,'linewidth',2)
        ylabel('Membrane Potential (mV)');
        title('Current Step Waveform')
        ax = gca; xax = ax.XAxis; set(xax,'visible','off')
        hold on
        xline(detection(1),'linestyle','--','color','blue','linewidth',2,'HandleVisibility','off') % detection start
        xline(detection(2),'linestyle','--','color','blue','linewidth',2) % detection end
        xline(baseline(1),'linestyle','--','color','green','linewidth',2,'HandleVisibility','off') % baseline start
        xline(baseline(2),'linestyle','--','color','green','linewidth',2) % baseline end
        hold off
        legend('Detection Region','Baseline Region','linewidth',1)
    
        % plot waveform of current step read from second channel of clampfit
        SI = ephysIO({clampfile,2}); % load in the current data
        x = SI.array(:,1); % x is time here
        pA_waveform = SI.array(:,2:end); % y is the array of current data
        Ibase = mean (mean(pA_waveform(baseStart:baseEnd,:))); % determine the baseline (intra-step)
        N = size(SI.array,2) - 1; % get the number of waves in the array (minus time)
        lo = dsearchn(x,detection(1)) + 1; % start of the test pulse
        hi = dsearchn(x,detection(2)) - 1; % end of the test pulse
        pA = zeros(1,N); % preallocate a blank series of steps
        for i = 1:N
           pA(i) = mean (pA_waveform(round(lo+(detEnd-detStart)*0.1):...
               round(hi-(detEnd-detStart)*0.1),i)); % fill the steps with the mean of each step
        end
        pA = fix((pA - Ibase) * 1e+12); % round to zero, baseline subtract and put into pA
        subplot(7,3,7); plot(x,(pA_waveform-Ibase)*1e12,'linewidth',1,'color','black','HandleVisibility','off')
        box off; set(gca,'linewidth',2); set(gcf,'color','white');
        xlabel('Time (s)'); ylabel('Command(pA)'); ylim([min(min(pA_waveform*1e12))-50,max(max(pA_waveform*1e12))+50])
        ax = gca; yax = ax.YAxis; set(yax,'TickDirection','out')
        hold on
        xline(detection(1),'linestyle','--','color','blue','linewidth',2,'HandleVisibility','off') % detection start
        xline(detection(2),'linestyle','--','color','blue','linewidth',2) % detection end
        xline(baseline(1),'linestyle','--','color','green','linewidth',2,'HandleVisibility','off') % baseline start
        xline(baseline(2),'linestyle','--','color','green','linewidth',2) % baseline end
        hold off
    
        % catch for scaling error
        if sum(pA) > 10000
            pA = [-200:20:400];
        else
        end

         % catch me if user doesn't select input arguments at the start
        if numInputs < 3
            % prompt the user to select the input here
            disp('enter current value (pA) to analyse spike frequency adaptation at')
            prompt = {...
                'Enter adaptation stimulation value (pA):'};
            dlg_title = 'Adaptation stimulation value input';
            num_lines = 1;
            def = {'300'};
            answer  = inputdlg(prompt,dlg_title,num_lines,def,'on');
            answer = str2double(answer);
            adaptI = answer(1);
        end
        % Plot the wave of interest, in our case closest to 300 pA
        % stimulation
        [val,ind] = min(abs(pA-adaptI));
        subplot(7,3,[2,5,8]);
        plot(Time,Waves(:,ind)*1000,'color','red'); hold on; plot(Time,Waves(:,zerowave)*1000,'color','black')
        box off; set(gcf,'color','white'); set(gca,'linewidth',2)
        xlabel('Time (s)'); ylabel('Membrane Potential (mV)');
        lgd = legend(char(string(pA(ind))),char(string(pA(zerowave))),'linewidth',1);
        title(lgd,'Current (pA)')
        title('Exemplary Waves')
    
        % extract and plot peak coordinates
        peakInd = cell2mat(locs(ind));
        peakVal = cell2mat(pks(ind));
        peakTimes = zeros(size(peakInd));
        peakAmps = zeros(size(peakVal));
        for pI = 1:size(peakTimes,1)
            peakTimes(pI) = Time(peakInd(pI)+detStart);
            peakAmps(pI) = peakVal(pI)*1000;
        end
        plot(peakTimes,peakAmps,'--ob','DisplayName','Peak Intervals')
    
        % convert to instantaneous frequency and plot
        ISI = diff(peakTimes);
        peakHz = 1./ISI;
        subplot(7,3,[3,6,9])
        plot(peakHz,'-o',...
            'LineWidth',1,...
            'MarkerSize',10,...
            'MarkerEdgeColor','r')
        box off; set(gcf,'color','white'); set(gca,'linewidth',2)
        xlabel('Interspike Interval #'); ylabel('Instantaneous Frequency (Hz)');
        legend('Spike Frequency','linewidth',1)
    
        % calculate adaptation index and plot
        adaptIndex = ISI(1)/ISI(end);
        subplot(7,3,[13,16,19]);
        names = char('','Initial','Final','');
        x = [2,3];
        y = [ISI(1),ISI(end)];
        plot(x,y,'--ob',...
            'LineWidth',1.5,...
            'MarkerSize',6,...
            'MarkerEdgeColor','r')
        box off; set(gcf,'color','white'); set(gca,'linewidth',2)
        set(gca,'XTick',1:4,'XTickLabel',names)
        xlim([1,4]);
        xlabel('Interspike Interval #'); ylabel('Insterspike Interval (s)');
        txt = {['Adaptation Index: '],[num2str(adaptIndex)]};
        annotation('textbox','String',txt,'Position',subplot(7,3,[13,16,19]).Position,'VerticalAlignment','top','HorizontalAlignment','right','FitBoxToText','on','linewidth',1.5);
    
        % plot as inst freq over time with raster and current injection
        subplot(7,3,[14:15,17:18,20:21]);
        fTime(1:2) = [0,peakTimes(1)-0.0001];
        for fT = 1:size(peakTimes,1)
            fTime(2+fT) = peakTimes(fT);
        end
        fVals(1:3) = [0,0,0];
        for fV = 1:size(peakHz,1)
            fVals(3+fV) = peakHz(fV);
        end
        plot(fTime,fVals,'linewidth',1.5,'DisplayName','Instantaneous Frequency (Hz)')
        xlabel('Time (s)'); ylabel('Instantaneous Frequency (Hz)');
        hold on; plot(peakTimes(1),65,"|",'color','r','linewidth',3,'MarkerSize',10,'DisplayName','Spikes')
        hold on; plot(peakTimes,65,"|",'color','r','linewidth',3,'MarkerSize',10,'HandleVisibility','off')
        box off; set(gcf,'color','white'); set(gca,'linewidth',2)
        x = [0.5,1.5]; y = [-10,-10]; hold on; line(x,y,'linewidth',10,'color','black','DisplayName','Current Injection')
        legend('linewidth',1)
        xlim([0.4,1.4]); ylim([-20,70])
   
        %% Bridge Balance adjustments
        % apply necessary bridge balance adjustments if the user requested offline
        
        % plot the initial current step so user can see if recording was balanced
        t = figure; plot(Time(1:151),Waves(300:450,1)); box off; xlabel('Time (s)');
        title('Initial Current Step [Zoomed]'); set(gca,'linewidth',2); 
        set(gcf,'color','white'); ylabel('Membrane Potential (mV)'); 
        
        % was this recording bridge balanced appropriately?
        dlgTitle    = 'Bridge Balance';
        dlgQuestion = 'Was this recording appropriately bridge balanced?';
        balanced = questdlg(dlgQuestion,dlgTitle,'Yes','No','Yes');
        
        close(t) % closes the test fig as it's not interesting anymore
        
        % logical fork following balancing
        if balanced == "Yes"
            disp('Recording appropriately balanced, no further action required')
            Vm_adjust = NaN(size(Waves,2),1);
            Rs_Init = NaN;
            offline_BB_performed = "No";
        elseif balanced == "No"
            disp('Recording not appropriately balanced, performing offline bridge balance')
            offline_BB_performed = "Yes";
            % calculate Rs values from the initial current step
            Istep = -60; % pA
            % find the triple diff peak
            [~, ind] = findpeaks((gradient(gradient(Waves(:,1))))*-1,...
                'minPeakProminence',2.5e-4,'NPeaks',1);
            % preallocate
            Vm_1 = zeros(size(Waves,2),1);
            Vm_2 = zeros(size(Waves,2),1);
            delta_Vm = zeros(size(Waves,2),1);
            Vm_adjust = zeros(size(Waves,2),1);
            for i = 1:size(Waves,2)
                Vm_1(i) = mean(Waves(1:ind,i)); % in V
                Vm_2(i) = mean(Waves(405:410,i)); % in V
                delta_Vm(i) = abs(Vm_1(i)-Vm_2(i))*1000; % in mV
            end
            Rs_Init = abs(median(delta_Vm/Istep)*1000); % in MOhm
            for i = 1:size(Waves,2)
                Vm_adjust(i) = ((Rs_Init*1e6)*(pA(i)*1e-12)); % in V (initial)
            end
            
            % plot the required adjustment
            figure; plot(pA, Vm_adjust*1000,'-o'); box off; set(gca,'linewidth',2); 
            set(gcf,'color','white'); xlabel('pA'); ylabel('Adjustment req. (mV)')
            title('Offline Bridge Balance Adjustment Required')
            
            % balance the data here
            % IR adjustment
            deltaV = abs((mean(Waves(IR_start:IR_end,1))-Vm_adjust(1)) ...
                - (mean(Waves(IR_start:IR_end,2))-Vm_adjust(2))); % Delta_Voltage (Volts)
            I = abs(pA(1)*1e-12 - pA(2)*1e-12); % I (Amps)
            R = deltaV / I; % R (Ohms)
            IR = R / 1e6; % R (MegaOhms)
            % replot the whole of that IR figure
            figure(fh); delete(IRtext)
                subplot(7,4,[19,23,27]); plot(Time,Waves(:,1)*1000,'color','black'); hold on; plot(Time,Waves(:,3)*1000,'color','red')
                box off; set(gcf,'color','white'); set(gca,'linewidth',2)
                xlabel('Time (s)'); ylabel('Membrane Potential (mV)');
                lgd = legend(char(string(pA(1))),char(string(pA(3))),...
                    'linewidth',1,...
                    'location','southeast',...
                    'AutoUpdate','off');
                title(lgd,'Current (pA)')
                title('Input Resistance Calculation')
                %show lines for region of IR determination
                hold on; xline(1.5,'--'); xline(1.3,'--')
                
                % Calculate Input Resistance in MegaOhms
                txt = {['\bf Input Resistance: '],[num2str(IR) ' M\Omega \rightarrow']};
                hold on; IRtext = text(0.7,(mean(Waves(IR_start:IR_end,2))*1000) + 15,txt);
            % thresh
            Threshold = Threshold-(Vm_adjust(wavenum_first)*1000);
            % afterhyp
            Afterhyperpolarisation = Afterhyperpolarisation-(Vm_adjust(wavenum_first)*1000);
            % peak
            Overshoot = Overshoot-(Vm_adjust(wavenum_first)*1000);
            % subAP_Vm 
            subAP_Vm = subAP_Vm-(Vm_adjust(1:size(subAP_Vm,2)));
        
        end
        %% Output
        
        % Would you like to add any notes? 
        prompt = {'Would you like to add any notes'};
        dlgtitle = 'Notes';
        definput = {'eg. dodgy input'};
        dims = [1 40];
        notes = inputdlg(prompt,dlgtitle,dims,definput);
        
        % create output structure
        output.filepath = path;
        output.steps = pA;
        output.waveform = pA_waveform;
        output.time = Time;
        output.waves = Waves;
        output.ephysIO = S;
        output.numSpikes = numSpikes;
        output.Rh = Rh;
        output.sag_mV = Ih_Sag_Amp;
        output.sag_ratio = Ih_Sag_Percentage;
        output.peak = Overshoot;
        output.afterhyp = Afterhyperpolarisation;
        output.amp = Amplitude;
        output.thresh = Threshold;
        output.half = Halfwidth;
        output.rise = Rise;
        output.fall = Fall;
        output.IR = IR;
        output.Vm = Vm;
        output.Vm_stability = Vm_stability;
        output.subAP_Vm = subAP_Vm;
        output.peakTimes = peakTimes;
        output.peakAmps = peakAmps;
        output.ISI = ISI;
        output.peakHz = peakHz;
        output.adaptIndex = adaptIndex;
        output.Offline_BB = Vm_adjust;
        output.Online_BB_performed = balanced;
        output.Offline_BB_performed = offline_BB_performed;
        output.Notes = notes;
        output.LPF_Hz = LPF_Hz;
        output.winSize_ms = winSize_ms;
        
        % navigate to root dir
        cd(newPath) 
        %chdir(fullfile('..','..'))
        
        dlgTitle    = 'Save output';
        dlgQuestion = 'Would you like to save the outputs?';
        saveout = questdlg(dlgQuestion,dlgTitle,'Yes','No','Yes');
        
        if saveout == "Yes"
            f = waitbar(0,'Saving  ...');
            set(f,'Name','Saving output, do not close figures');
            % save output
            outname = split(strtrim(clampfile),filesep);
            outname = char(string(outname(end-2))); % name the output the same as the folder the recording came from
            if figsave == 1
                saveas(fh,[outname,'_master.fig']); % save the master fig
                    waitbar(.2,f,'Saving master figure ...');
                saveas(fh2,[outname,'_subAP.fig']); % save the subAP_Vm fig
                    waitbar(.4,f,'Saving sub action potential Vm figure ...');
                saveas(fh3,[outname,'_Vm_stability.fig']); % save the subAP_Vm fig
                    waitbar(.6,f,'Saving Vm stabilitiy figure ...');
                saveas(fh4,[outname,'_adapt.fig']);
                    waitbar(.8,f,'Saving adaptation figure ...')
            else
            end
            save([outname,'.mat'],'output')
                waitbar(1,f,'Saving final data structure ... ');
                close(f)
        else
        end
        % return to if loop from the top 
    else
        close(fh) % close figure
        disp('---------')
        disp('Analysis Aborted, please find a new recording')
        disp('---------')
end
end
end


%% Dependencies
%  Function File: eventer
%
%  peak = eventer(file,TC,s,SF)
%  peak = eventer(file,TC,s,SF,...,'method',method)
%  peak = eventer(file,TC,s,SF,...,'exclude',exclude)
%  peak = eventer(file,TC,s,SF,...,'criterion',criterion)
%  peak = eventer(file,TC,s,SF,...,'rmin',rmin)
%  peak = eventer(file,TC,s,SF,...,'hpf',hpf)
%  peak = eventer(file,TC,s,SF,...,'lpf',lpf)
%  peak = eventer(file,TC,s,SF,...,'taus',taus)
%  peak = eventer(file,TC,s,SF,...,'baseline',baseline)
%  peak = eventer(file,TC,s,SF,...,'win',win)
%  peak = eventer(file,TC,s,SF,...,'lambda',lambda)
%  peak = eventer(file,TC,s,SF,...,'config',config)
%  peak = eventer(file,TC,s,SF,...,'average',average)
%  peak = eventer(file,TC,s,SF,...,'export',format)
%  peak = eventer(file,TC,s,SF,...,'channel',channel)
%  peak = eventer(file,TC,s,SF,...,'wave',wave)
%  peak = eventer(file,TC,s,SF,...,'exmode',wave)
%  peak = eventer(file,TC,s,SF,...,'figure',format)
%  peak = eventer(file,TC,s,SF,...,'globvar',globvar)
%  [peak,IEI] = eventer(...)
%
%  peak = eventer(file,TC,s,SF) returns the amplitudes of spontaneous EPSC-
%    or EPSP-like waveforms, detected either from the first derivative or by
%    FFT-based deconvolution [1]. In both cases, a template modeled by the
%    sum of two exponentials is required, whose time constants must be
%    provided in units seconds as a vector (TC). The sign (s) of the event
%    peak deflections in the wave file must be specified as '-' or '+'.
%    Event times are then defined where there are maxima of delta-like waves
%    exceeding a threshold, which is set to the standard deviation of the
%    noise of the first derivative (or deconvoluted wave) multiplied by a
%    scale factor (SF). Events are then extracted as episodic data and the
%    factor parameter of least-squares fits with the template define the
%    peak amplitudes. See the associated input-output function ephysIO for
%    details of supported input file formats. The file extension must be
%    included in the filename of the file input argument. The template
%    kinetics is modeled by the difference of 2 exponentials:
%      f(t) = exp ( - t / tau_decay ) - exp ( - t / tau_rise )
%
%  peak = eventer(file,TC,s,SF,...,'method',method) sets the detection
%    method of eventer. The options available are 'deconvolution' (default)
%    or 'derivative'. At very low sampling frequencies, the derivative
%    method is preferable, in which case the template is only used for
%    refitting events and applying the criterion for event screening.
%
%  peak = eventer(file,TC,s,SF,...,'exclude',exclude) sets exclusion zones,
%    which must be specified as a 2-column matrix, where the first and
%    second columns define the start and end times of the excluded zones
%    (in seconds). Events detected in these regions are discarded. By
%    default, there are no exclusion zones.
%
%  peak = eventer(file,TC,s,SF,...,'criterion',criterion) sets the type of
%    correlation coefficient used when comparing the fitted model template
%    with each event. Options are 'Pearson', 'Spearman' or 'Kendall', or
%    at an object of class TreeBagger for event classification (machine
%    learning).
%
%  peak = eventer(file,TC,s,SF,...,'rmin',rmin) sets the correlation
%    coefficient for fitted model template and the event. Events with a
%    correlation coefficient (r) less than the value defined in rmin are
%    discarded [2]. The default rmin is 0.4, where r < 0.4 is considered
%    a very weak correlation.
%
%  peak = eventer(file,TC,s,SF,...,'hpf',hpf) sets the -3 dB cut-off (in
%    Hz) of the low-pass median filter applied to the deconvoluted wave,
%    where the filtered wave is then subtracted from the deconvoluted
%    wave. Thus, this is essentially a high-pass filter. The algorithm
%    implements a bounce correction to avoid end effects. The default
%    cut-off is 1 Hz.
%
%  peak = eventer(file,TC,s,SF,...,'lpf',lpf) sets the -3 dB cut-off (in
%    Hz) of the low-pass binomial filter applied to the deconvoluted wave.
%    The default cut-off is 200 Hz.
%
%  peak = eventer(file,TC,s,SF,...,'taus',taus) sets the number of time
%    constants after the peak of the template to use when fitting the
%    template to the detected events. The default is 0.4, which corresponds
%    to the time for the template to decay by 25 % of the peak. This
%    limited fit ensures that peak measurements are not compromised at high
%    event frequencies where there is event overlap.
%
%  peak = eventer(file,TC,s,SF,...,'baseline',baseline) sets the length
%    of time to use as the pre-event baseline for the template fit in
%    seconds.
%
%  peak = eventer(file,TC,s,SF,...,'win',win) sets the event window limits
%    in seconds for the conversion of the continuous wave to episodic
%    data. By default the limits are set at [-0.01 0.04].
%
%  peak = eventer(file,TC,s,SF,...,'lambda',lambda) sets the damping
%    factor used in the Levenberg-Marquardt ordinary non-linear least-
%    squares fitting procedures. The default is 1. The higher the value
%    of lambda, the more robust the fitting procedure is to the initial
%    values but the greater the number of iterations it will use.
%
%  peak = eventer(file,TC,s,SF,...,'config',config) sets the configuration
%    of the recording wave to either 'VC' (for voltage clamp) or 'CC'
%    (for current-clamp).The default is 'VC'.
%
%  peak = eventer(file,TC,s,SF,...,'average',average) sets eventer to
%    compute either the ensemble 'mean' or 'median' of the merged event
%    data. The default is 'mean'.
%
%  peak = eventer(file,TC,s,SF,...,'export',format) sets eventer to
%    export the episodic wave data of all detected events in the
%    specified format.
%
%  peak = eventer(file,TC,s,SF,...,'channel',channel) sets eventer to
%    select the recording channel from the file. If none is specified,
%    eventer analyses channel number 1. This argument is ignored for
%    loaded filetypes that do not support multiple recording channels.
%
%  peak = eventer(file,TC,s,SF,...,'wave',wave) sets eventer to select
%    the wave number from the file. If none is specified, eventer analyses
%    the first wave.
%
%  peak = eventer(file,TC,s,SF,...,'exmode',exmode) tells eventer what to
%    do with the first event proceeding each exclusion zone. For mode = 1
%    (default) IEIs are calculated for these events from the last event
%    preceeding the exclusion zone under the assumption that no events
%    occurred during the exclusion zone. For mode = 2, these events are
%    assigned an IEI of NaN. Note that events with NaN values are excluded
%    during the merge (i.e. that are in the ALL_events output directory).
%
%  peak = eventer(file,TC,s,SF,...,'figure',format) tells eventer what
%    format to use when saving figures. Bitmap images are saved at 300 dpi
%    resolution. Accepted formats are: 'fig', 'tiff', 'tiffn', 'png', 'bmp',
%    'eps', 'pdf', 'svg' and 'emf'. By default, 'tiff' will save images with
%    compression. To save files without compression use format 'tiffn'. On
%    windows platforms, images can also be saved as enhanced meta files
%    ('emf') for easy import into Microsoft applications.
%
%  peak = eventer(file,TC,s,SF,...,'globvar',globvar) tells eventer how to
%    load the data. Set globvar to 0 to force eventer to load the data from
%    the specified file. Set globvar to 1 to force eventer to use the data
%    already loaded into the workspace by ephysIO.
%
%  [peak,IEI] = eventer(...) returns the interevent intervals preceding
%    each peak event. The first event from the start of the wave is assigned
%    an IEI of NaN.
%
%  See the example distributed with this function.
%
%  Dependencies: sinv, binomialf, bounce, hpfilter, lpfilter, lsqfit and
%                ephysIO.
%
%  References:
%   [1] Pernia-Andrade AJ, Goswami SP, Stickler Y, Frobe U, Schlogl A,
%    Jonas P. (2012) A deconvolution-based method with high sensitivity
%    and temporal resolution for detection of spontaneous synaptic currents
%    in vitro and in vivo. Biophys J. 103(7):1429-39.
%   [2] Jonas P, Major G and Sakmann (1993) Quantal components of unitary
%    EPSCs at the mossy fibre synapse on CA3 pyramidal cells of the rat
%    hippocampus. J Physiol. 472:615-663.
%
%  eventer v1.6 (last updated: 18/12/2020)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/


function [peak,IEI,features] = eventer(arg1,TC,s,SF,varargin)

  if nargin<2 || sum(size(TC))>3
    error('A vector defining the time constants of the template has to be specified');
  else
    if sum(sign(TC))~=2
      error('The time constants must be non-zero and non-negative');
    end
    if TC(1)>=TC(2)
      error('The first time constant (rise) must be smaller than the second time constant (decay)');
    end
  end

  if nargin<3
    error('The sign of the peaks has to be specified');
  end

  if s~='+' && s~='-'
    error('The sign of the peaks has to be specified in the + or - direction');
  end

  if nargin<4
    error('The factor of standard deviations has to be specified to set the detection threshold');
  else
    if isempty(SF)
      SF = 4;
    end
  end

  % Set additional options
  options = varargin;
  method = 1+find(strcmp('method',options));
  excl = 1+find(strcmp('exclude',options));
  rmin = 1+find(strcmp('rmin',options));
  criterion = 1+find(strcmp('criterion',options));
  hpf = 1+find(strcmp('hpf',options));
  lpf = 1+find(strcmp('lpf',options));
  taus = 1+find(strcmp('taus',options));
  base = 1+find(strcmp('baseline',options));
  win = 1+find(strcmp('win',options));
  config = 1+find(strcmp('config',options));
  lambda = 1+find(strcmp('lambda',options));
  average = 1+find(strcmp('average',options));
  export = 1+find(strcmp('export',options));
  channel = 1+find(strcmp('channel',options));
  wave = 1+find(strcmp('wave',options));
  exmode = 1+find(strcmp('exmode',options));
  figform = 1+find(strcmp('figure',options));
  globvar = 1+find(strcmp('globvar',options));
  merge = 1+find(strcmp('merge',options));
  showfig = 1+find(strcmp('showfig',options));
  absT = 1+find(strcmp('threshold',options));
  if ~isempty(method)
    try
      method = options{method};
    catch
      method = 'Deconvolution';
    end
  else
    method = 'Deconvolution';
  end
  if ~isempty(excl)
    try
      excl = options{excl};
    catch
      excl = [];
    end
  else
    excl = [];
  end
  if ~isempty(rmin)
    try
      rmin = options{rmin};
    catch
      rmin = 0.4;
    end
  else
    rmin = 0.4;
  end
  if ~isempty(criterion)
    try
      criterion = options{criterion};
    catch
      criterion = 'Pearson';
    end
  else
    criterion = 'Pearson';
  end
  if strcmpi(class(criterion),'TreeBagger')
    coeff = 'Pearson';
    rmin = -1;
  else
    coeff = criterion;
  end
  % Default value of 0.4 corresponds to the time for the event to decay by 25% of the peak
  if ~isempty(taus)
    try
      taus = options{taus};
    catch
      taus = 0.4;
    end
  else
    taus = 0.4;
  end
  if ~isempty(hpf)
    try
      hpf = options{hpf};
    catch
      hpf = 1;
    end
  else
    hpf = 1;
  end
  if ~isempty(lpf)
    try
      lpf = options{lpf};
    catch
      lpf = 200;
    end
  else
    lpf = 200;
  end
  if ~isempty(base)
    try
      base = abs(options{base});
    catch
      base = 1e-03;
    end
  else
    base = 1e-03;
  end
  if ~isempty(win)
    try
      win = options{win};
    catch
      win = [-0.01 0.04];
    end
  else
    win = [-0.01 0.04];
  end
  if ~isempty(config)
    try
      config = options{config};
    catch
      config = '';
    end
  else
    config = '';
  end
  if ~isempty(lambda)
    try
      lambda = options{lambda};
    catch
      lambda = 1;
    end
  else
    lambda = 1;
  end
  if ~isempty(average)
    try
      average = options{average};
    catch
      average = 'mean';
    end
  else
    average = 'mean';
  end
  if ~isempty(export)
    try
      export = options{export};
    catch
      export = 'none';
    end
  else
    export = 'none';
  end
  if ~isempty(channel)
    try
      channel = options{channel};
    catch
      channel = 'none';
    end
  else
    channel = 1;
  end
  if ~isempty(wave)
    try
      wave = options{wave};
    catch
      wave = 1;
    end
  else
    wave = 1;
  end
  wave = wave+1; % Add 1 since first column in data matrix is time
  if ~isempty(exmode)
    try
      exmode = options{exmode};
    catch
      exmode = 1;
    end
  else
    exmode = 1;
  end
  if ~isempty(figform)
    try
      figform = options{figform};
    catch
      figform = 'fig';
    end
  else
    figform = 'fig';
  end
  if strcmpi(figform,'eps')
    figform = 'epsc';
  end
  if strcmpi(figform,'emf')
    if ~ispc
      error('the enhanced meta file format is only supported on windows platforms');
    end
    figform = 'meta';
  end
  if ~isempty(globvar)
    try
      globvar = options{globvar};
    catch
      globvar = 0;
    end
  else
    globvar = 0;
  end
  if ~isempty(merge)
    try
      merge = options{merge};
    catch
      merge = 1;
    end
  else
    merge = 1;
  end
  if ~isempty(showfig)
    try
      showfig = options{showfig};
    catch
      showfig = 'on';
    end
  else
    showfig = 'on';
  end
  if ~isempty(absT)
    try
      absT = options{absT};
    catch
      absT = 0;
    end
  else
    absT = 0;
  end

  % Error checking
  if size(win,1)~=1 && size(win,2)~=2
      error('win must be a vector defining the limits of the event window');
  end

  % Close existing figures
  if ishandle(1)
    close(1)
  end
  if ishandle(2)
    close(2)
  end
  if ishandle(3)
    close(3)
  end
  if ishandle(4)
    close(4)
  end
  if ishandle(5)
    close(5)
  end
  if ishandle(6)
    close(6)
  end

  % Load data
  cwd = pwd;
  if isstruct(arg1)
    file = 'Data'; %#ok<*NASGU>
    pathstr = '';
    filename = 'Data';
    ext = '';
    data  = arg1.array;
    xdiff = arg1.xdiff; %#ok<*GPFST>
    xunit = arg1.xunit;
    yunit = arg1.yunit;
    names = arg1.names;
    notes = arg1.notes;
  elseif ischar(arg1)
    file = arg1;
    if ~isempty(regexpi(file(end-2:end),'.gz'))
      [pathstr,filename,ext] = fileparts(file(end-2:end));
    elseif ~isempty(regexpi(file(end-3:end),'.zip'))
      [pathstr,filename,ext] = fileparts(file(end-3:end));
    else
      [pathstr,filename,ext] = fileparts(file);
    end
    if ~isempty(pathstr)
      chdir(pathstr);
    end
    if globvar == 0
      [data,xdiff,xunit,yunit,names,notes] = ephysIO ({strcat(filename,ext),channel});
    elseif globvar == 1
       global array %#ok<*TLEV>
       data = array;
       global xdiff
       global xunit
       global yunit
       global names
       global notes
    end
  end

  if wave > size(data,2)
    error('wave number exceeds data dimensions')
  end
  filewave = sprintf('%s_ch%d_%s',filename,channel,names{wave});
  if isempty(xunit)
    warning('xunit is undefined. X dimenion units assumed to be seconds.')
    xunit = 's';
  elseif ~strcmp(xunit,'s')
    error('expected xunit to be seconds (s)')
  end
  % Overide data units with configurationn argument
  if strcmpi(config,'CC')
    yunit = 'V';
  elseif strcmp(config,'VC')
    yunit = 'A';
  elseif strcmp(config,'')
    yunit = '';
  end
  try
    if xdiff == 0
      warning('data array must have a constant sampling interval');
    end
  catch
    if strcmp(xdiff,'variable')
      warning('data array must have a constant sampling interval');
    end
  end

  % Assign variables from data
  N = size(data,1);
  RecordTime = range(data(:,1));
  sample_rate = round((N-1)/RecordTime);
  tau_rise = TC(1);
  tau_decay = TC(2);
  % Concatenate traces
  numwave = numel(wave);
  if numel(wave) > 1
    Trace = reshape(data(:,wave),numwave*N,1);
    T = (0:numwave*N-1)'*1/sample_rate;
    t = T;
  else
    Trace = data(:,wave);
    T = (0:N-1)'*1/sample_rate;
    t = data(:,1);
  end
  clear data

  % Create model template event
  Template = -exp(-T/tau_rise)+exp(-T/tau_decay);
  TemplatePeak = max(Template);
  Template = Template/TemplatePeak;
  time2peak = tau_decay*tau_rise/(tau_decay-tau_rise)*log(tau_decay/tau_rise);
  if s=='-'
    Template = Template*-1;
  end

  if strcmpi(method, 'Deconvolution')
    % Perform fourier transform-based deconvolution (pad ends to mask end effects)
    Template_FFT = fft(cat(1,Template,zeros(2*sample_rate,1)));
    Trace_FFT = fft(bounce(Trace,sample_rate));
    DEC = real(ifft(Trace_FFT./Template_FFT));
    clear Trace_FFT Template_FFT
    DEC(1:sample_rate) = [];
    DEC(end-sample_rate+1:end) = [];
  elseif strcmpi(method, 'Derivative')
    DEC = [diff(Trace);0];
  end

  % Band-pass filter the deconvoluted trace (default is 1-200 Hz)
  DEC = filter1 (DEC, t, hpf, lpf, 'median');

  % Assign NaN to deconvoluted waves for values inside user-defined exclusion zones
  % Calculate actual recording time analysed (not including exclusion zones)
  if ~isempty(excl)
    excl_idx(:,1) = dsearchn(t,excl(:,1));
    excl_idx(:,2) = dsearchn(t,excl(:,2));
  end
  for i=1:size(excl,1)
    DEC(excl_idx(i,1):excl_idx(i,2)) = NaN;
  end
  AnalysedTime = (sum(~isnan(DEC))-1)/sample_rate;

  % Create all-point histogram with optimal bin width determined by modified Silverman's
  % rule of thumb on the putative noise component based on robust statistics.
  % References:
  % -Hardle, W. SMOOTHING TECHNIQUES, with implementations in S. Springer, New York. 1991
  % -Silverman, B.W. DENSITY ESTIMATION FOR STATISTICS AND DATA ANALYSIS.
  %  Monographs on Statistics and Applied Probability, London: Chapman and Hall, 1986.
  noise = DEC(~isnan(DEC));
  lim = abs(min(noise));
  noise(noise>lim) = [];
  Q2 = median(noise);
  Q1 = median(noise(noise<Q2));
  Q3 = median(noise(noise>Q2));
  IQR = abs(Q1-Q3);
  if IQR < eps
    error('the distribution of the noise has zero standard deviation')
  end
  h = 0.79*IQR*numel(noise)^(-1/5);
  binwidth = 2*h;
  bins = ceil(range(DEC)/binwidth);
  try
    [counts,x] = hist(DEC,bins);
  catch
    error('hist returned an error. Try alterative filter settings')
  end

  % Find the center of the noise peak
  peak_count = max(counts);
  mu = mean(x(counts==max(counts)));

  % Use linear interpolation between histogram bins to find half-maximum and
  % extrapolate the FWHM from the bottom half of the noise peak. Use this to
  % approximate the standard deviation of the noise peak.
  top = x(counts>max(counts)/2);
  idx = ones(1,2)*find(x==top(1));
  idx(1) = idx(1)-1;
  LL = interp1(counts(idx),x(idx),peak_count/2,'linear','extrap');
  FWHM = abs((LL-mu))*2;
  sigma = FWHM/(2*sqrt(2*log(2)));

  % Center and scale the noise component of the deconvoluted data using initial estimates
  DEC = (DEC-mu)/sigma;
  x = (x-mu)/sigma;

  % Least-squares Gaussian fitting of the noise peak using Levenberg-Marquardt
  % algorithm with the initial estimates for the fit parameters (data is scaled
  % for best performance)
  xdata = x((x>-4)&(x<1))';
  ydata = counts((x>-4)&(x<1))'/max(counts);
  fun1 = @(p,xdata)p(1)*p(3)*sqrt(2*pi)*normpdf(xdata,p(2),p(3));
  p0 = [1;0;1];
  optimoptions = struct;
  optimoptions.Algorithm = {'levenberg-marquardt',lambda};
  optimoptions.TolX = eps;
  optimoptions.TolFun = eps;
  [p,resnorm,residual,exitflag] = lsqfit(fun1,p0,xdata,ydata,[],[],optimoptions);
  clear xdata ydata
  if exitflag<1
    error('Noise peak fitting failed to reach tolerance. Try a higher lambda value.')
  end

  % Recenter and scale the noise component of the deconvoluted data using the
  % optimization result
  DEC = (DEC-p(2))/p(3);
  x = (x-p(2))/p(3);
  
  % Absolute threshold provided. Overide scale factor setting.
  noiseSD = sigma*p(3);
  if absT > 0
    SF = absT/noiseSD;
  end
  
  % Scan superthreshold deconvoluted wave for local maxima (vectorized)
  N = numel(Trace); % Number of sample points in cropped wave
  npeaks = sum(diff(sign(diff(DEC)))==-2); % Total number of peaks without threshold
  Event = zeros(N,1);
  y = DEC.*(DEC>SF);
  Event(2:N-1) = diff(sign(diff(y)))==-2;

  % Create episodic data
  samples_pre = round(abs(win(1))*sample_rate);
  samples_post = round(win(2)*sample_rate);
  Event_idx = find(Event);
  idx_pre = Event_idx-samples_pre;
  idx_post = Event_idx+samples_post;
  n = sum(Event);
  %t_events(:,1) = win(1):1/sample_rate:win(2);
  t_events(:,1) = [-samples_pre:+samples_post]/sample_rate; % more robust
  y_events = NaN(length(t_events),n); % preallocation
  idx = cell(n,1);
  if n > 0
    for i=1:n
      if idx_pre(i) < 1
        % Event close to start of trace
        y_events(:,i) = Trace(1);
        idx{i} = transpose(1:idx_post(i));
        y_events(length(t_events)-idx_post(i)+1:end,i) = Trace(idx{i});
      elseif idx_post(i) > N
        % Event close to end of trace
        y_events(:,i) = Trace(end);
        idx{i} = transpose(idx_pre(i):N);
        y_events(1:numel(idx{i}),i) = Trace(idx{i});
      else
        idx{i} = transpose(idx_pre(i):idx_post(i));
        y_events(:,i) = Trace(idx{i});
      end
    end
  end
  clear idx
  Event_time = t(Event_idx);

  % Initialize output arguments
  peak = [];
  IEI = [];
  features = [];

  if n > 0
    % Fit template onto events using linear least squares
    SCALE = ones(1,n);
    OFFSET = zeros(1,n);
    % Time-to-peak = -log(tau_decay/tau_rise)/(1/tau_decay-1/tau_rise)
    if s=='+'
      idx = find(Template==max(Template));
    elseif s=='-'
      idx = find(Template==min(Template));
    end
    templateTime = idx+round(sample_rate*tau_decay*taus);
	if templateTime > samples_post
	  error('Event window size is to small for the events')
	end
    numelBase = round(base*sample_rate);
    A = [zeros(numelBase,1); Template(1:templateTime)];
    l = length(A);
    A = [ones(l,1),A];
    start = dsearchn(t_events,base*-1);
    t_fit = t_events(start:start+numelBase+templateTime-1);
    y_fit = y_events(start:start+numelBase+templateTime-1,:);
    y_template = NaN(size(y_fit)); % preallocation
    tpeak = NaN(n,1);  % preallocation
    t50 = NaN(n,1);    % preallocation
    tdecay = NaN(n,1); % preallocation
    auc = NaN(n,1);    % preallocation
    skew = NaN(n,1);   % preallocation
    rstdev = NaN(n,1); % preallocation
    rskew = NaN(n,1);  % preallocation
    r = NaN(n,1);      % preallocation
    before = NaN(n,1); % preallocation
    after = NaN(n,1);  % preallocation
    tzero = dsearchn(t_fit,0);
    for i=1:n
      % For fitting discard any region of overlap with next event
      if i < n
        k = min(Event_idx(i)+templateTime,Event_idx(i+1))-Event_idx(i)+numelBase;
      else
        k = min(Event_idx(i)+templateTime,N)-Event_idx(i)+numelBase;
      end
      % Perform linear least-squares fit using singular value decomposition on events
      kmin = sinv(A(1:k,:)'*A(1:k,:))*A(1:k,:)'*y_fit(1:k,i);
      OFFSET(i) = kmin(1);
      SCALE(i) = abs(kmin(2)); % Absolute value constraint: The template cannot be inverted
      y_template(:,i) = SCALE(i)*A(:,2)+OFFSET(i);
      r(i) = corr(y_template(1:k,i),y_fit(1:k,i),'type',coeff);
      residuals = y_fit(1:k,i) - y_template(1:k,i);
      [tpeak(i), t50(i), tdecay(i), auc(i), skew(i)] = shape(t_fit(1:k),y_fit(1:k,i),OFFSET(i),sample_rate,s,base);
      rstdev(i) = std(residuals);
      rskew(i) = skewness(residuals);
      after(i) = (k - numelBase)/sample_rate; % Truncated time to next event
      % Time-to-peak dead time. Discard events that are proceeded by another event
      % before the initial event reaches it's peaks
      if k < numelBase+idx
        r(i) = -inf;
      end
    end
    y_events = y_events-ones(length(t_events),1)*OFFSET;

    % Create feature matrix
    before = circshift(after,1); before(1) = templateTime / sample_rate;
    ampl = SCALE';
    features = [tpeak t50 tdecay auc skew ampl r rstdev before after];
    % Description of features
    %  1)  tpeak: time-to-peak
    %  2)  t50: duration of event above half-maximum amplitude (approximates FWHM)
    %  3)  tdecay: duration of event post-peak above half-maximum amplitude (approximates decay half-life)
    %  4)  auc: area-under-curve of the event
    %  5)  skew: skewness metric for the distribution of the event sample points around the baseline
    %  6)  ampl: event amplitude calculated from template fitting
    %  7)  r: correlation coefficient of the template fitting
    %  8)  rstdev: standard deviation of residuals from the fit
    %  9)  before: time from preceding event (truncated)
    %  10) after: time to next event (truncated)

    % Discard events that are poorly correlated with the fitted template (r < rmin)
    if strcmpi(class(criterion),'TreeBagger')
      out = predict(criterion,features);
      ridx = find(cell2mat(out)=='0');
      criterion = 'Machine Learning';
    else
      ridx = find(r<rmin);
    end
    Event(ridx) = 0; % Keep in case required for diagnostic tests
    Event_idx(ridx') = [];
    Event_time(ridx') = [];
    SCALE(ridx) = [];
    OFFSET(ridx) = []; % Keep in case required for diagnostic tests
    tpeak(ridx) = [];
    t50(ridx) = [];
    tdecay(ridx) = [];
    auc(ridx) = [];
    skew(ridx) = [];
    rstdev(ridx) = [];
    rskew(ridx) = [];
    r(ridx) = [];
    ampl(ridx) = [];
    before(ridx) = [];
    after(ridx) = [];
    y_events(:,ridx) = [];
    y_fit(:,ridx) = [];
    y_template(:,ridx) = [];
    features(ridx,:) = [];
    n = numel(Event_idx);
    y_avg = mean(y_fit,2);
    peak = SCALE';

    % Create wave of fitted templates to overlay onto the event wave
    % Note that the baseline period is not plotted
    FITtrace = nan(N,1);
    l = length(t_fit(tzero:end));
    for i=1:n
      k = min(Event_idx(i)+(l-1),N);
      FITtrace(Event_idx(i):k) = y_template(tzero:tzero+k-Event_idx(i),i);
      FITtrace(Event_idx(i)-1) = NaN;
    end
  end

  % Get screen size and set figure sizes
  set(0,'units','pixels');
  scn_sz = get(0,'screensize');
  sw = scn_sz(3);
  sh = scn_sz(4);
  fw = sw/3;
  fh = sh/2;

  % Figures
  vector = 'pdf eps epsc eps2 epsc2 meta svg ps psc ps2 psc2';
  bitmap = 'jpeg png tiff tiffn bmpmono bmp bmp16m bmp256 hdf pbm pbmraw pcxmono pcx24b pcx256 pcx16 pgm pgmraw ppm ppmraw';

  % Figure 1: All-points histogram of the deconvoluted data points
  clear h
  h1 = figure(1);set(h1,'visible',showfig);set(h1,'OuterPosition',[0 sh-fh fw fh]);
  addToolbarExplorationButtons(gcf);
  bar(x,counts,1,'EdgeColor','b','Facecolor','b');
  ax = gca;
  ax.Toolbar.Visible = 'off';
  hold on; plot(x,max(counts)*fun1([p(1);0;1],x),'r','linewidth',3);
  ylimits=ylim;
  plot(SF*ones(1,2),ylimits,'g','linewidth',2)
  hold off;
  ylim(ylimits); xlim([min(DEC),max(DEC)]);
  box('off'); grid('off');
  title('Histogram of the deconvoluted wave');
  xlabel('Deconvoluted wave (SD)');
  ylabel('Number of points');
  xlim([-5,20]);  % Set x-axis limits to -5 and +20 SD of the baseline noise for clarity

  % Figure 2: Filtered deconvoluted wave (blue) and detection threshold (green)
  h2 = figure(2);set(h2,'visible',showfig);set(h2,'OuterPosition',[fw sh-fh fw fh]);
  addToolbarExplorationButtons(gcf);
  if ~isempty(regexpi(vector,figform))
    try
      reduce_plot(t,DEC,'b','linewidth',1.25);
    catch
      plot(t,DEC,'b');
    end
  else
    plot(t,DEC,'b');
  end
  ax = gca;
  ax.Toolbar.Visible = 'off';
  hold on;
  plot([t(1),t(end)],SF*ones(1,2),'g','linewidth',2);
  p2 = get(gca);
  EventPoints = p2.YLim(2)*Event;
  if ~isempty(regexpi(vector,figform))
    try
      reduce_plot(t(Event>0),EventPoints(Event>0),'r.','markersize',5);
    catch
      plot(t(Event>0),EventPoints(Event>0),'r.','markersize',5);
    end
  else
    plot(t(Event>0),EventPoints(Event>0),'r.','markersize',5);
  end
  hold off;
  xlim([min(t),max(t)]);
  grid('off'); box('off');
  title('Deconvoluted wave (before applying event criterion)');
  ylabel('Deconvoluted wave (SD)');
  xlabel('Time (s)');

  % Print result summary and escape from the function if no events were detected
  if n == 0
    if strcmpi(class(criterion),'TreeBagger')
      criterion = 'Machine Learning';
    end

    % Print basic information
    %fprintf(fid,...
    %        ['--------------------------EVENTER---------------------------\n',...
    %         ' Automatic PSC/PSP detection using FFT-based deconvolution\n',...
    %         ' and event peak analysis by least-squares template fitting\n',...
    %         ' Version v1.0 Copyright 2014 Andrew Charles Penn                \n\n'])
    fid=fopen(filewave,'w+t');
    fprintf(fid,'Filename: %s\n',filename);
    fprintf(fid,'Channel: %d\n',channel);
    fprintf(fid,'Wave number: %d\n',wave-1);
    fprintf(fid,'Wave name: %s\n',names{wave});
    fprintf(fid,'Total number of events detected: %d\n',n);
    fprintf(fid,'Duration of recording analysed (in s): %.1f\n',AnalysedTime);
    fprintf(fid,'High-pass filter cut-off on deconvoluted wave (at -3 dB, in Hz): %.3g\n',hpf);
    fprintf(fid,'Low-pass filter cut-off on deconvoluted wave (at -3 dB, in Hz): %.3g\n',lpf);
    fprintf(fid,'Vector of model template time constants (in s): [%.3g,%.3g]\n',TC*1e3);
    fprintf(fid,'Standard deviation of the noise of the deconvoluted wave (a.u.): %.3g\n',noiseSD);
    fprintf(fid,'Scale factor of noise standard deviations for threshold setting: %.3g\n',SF);
    fprintf(fid,'Theoretical false positive detection rate before applying event criterion (in Hz): %.3g\n',(1-normcdf(SF))*sample_rate);
    fprintf(fid,'Sign of the event peaks: %s\n',s);
    fprintf(fid,'Criterion used for event screening: %s\n',criterion);
    fprintf(fid,'Maximum time after peak (expressed in time constants) used for template fit: %.3g\n',taus);
    fprintf(fid,'Dead time from event start required for template fit (ms): %.3g\n',time2peak*1e3);
    fprintf(fid,'Minimum acceptable correlation coefficient for the template fit: %.3g\n',rmin);
    fprintf(fid,'Sample rate of the recording (in kHz): %d\n',sample_rate*1e-03);
    fprintf(fid,'lsqfit exitflag for fitting the noise peak: %d\n',exitflag);
    fprintf(fid,'Exclusion zones:\n');
    for i=1:size(excl,1)
      fprintf(fid,'%.6f\t%.6f\n',[excl(i,1) excl(i,2)]);
    end
    fclose(fid);

    % Save episodic data waves to file
    if exist('eventer.output','dir')==0
      mkdir('eventer.output');
    end
    cd eventer.output
    if exist(filewave,'dir')==0
      mkdir(filewave);
    end
    chdir(filewave);

    % Save detection parameters and event amplitudes and times to text file
    dlmwrite('_parameters',[TC';noiseSD;AnalysedTime;sample_rate;npeaks],'newline','pc')
    dlmwrite('_offset',win(1),'newline','pc');
    if exist('summary.txt','file')~=0
      delete('summary.txt');
    end
    movefile(['../../' filewave],'summary.txt');
    cd ../..

    if merge == 1
      merge_data(average,s,win,export,optimoptions,cwd,figform,config,taus);
    end

    return

  end

  % Calculate axes autoscaling for figure 3
  y_autoscale = 0.1*max(range(y_events));
  y_maxlim = max(max(y_events))+y_autoscale;
  y_minlim = min(min(y_events))-y_autoscale;
  % Figure 3: Identified events (black) aligned and overlayed with the mean event (red)
  h3 = figure(3);set(h3,'visible',showfig);set(h3,'OuterPosition',[fw*2 sh-fh fw fh]);
  addToolbarExplorationButtons(gcf);
  ylimits = [y_minlim,y_maxlim]; % Encoded y-axis autoscaling
  if ~isempty(regexpi(vector,figform))
    try
      reduce_plot(t_events,y_events,'color',[0.85,0.85,0.85]);
    catch
      plot(t_events,y_events,'color',[0.85,0.85,0.85]);
    end
  else
    plot(t_events,y_events,'color',[0.85,0.85,0.85]);
  end
  ax = gca;
  ax.Toolbar.Visible = 'off';
  hold on;
  if ~isempty(regexpi(vector,figform))
    try
      reduce_plot(t_events,mean(y_events,2),'-b');
    catch
      plot(t_events,mean(y_events,2),'-b');
    end
  else
    plot(t_events,mean(y_events,2),'-b');
  end
  hold off;
  xlim(win);
  ylim([ylimits(1),ylimits(2)]);
  title('Events');
  xlabel('Time (s)');
  if strcmp(yunit,'A')
    ylabel('Current (A)');
  elseif strcmp(yunit,'V')
    ylabel('Voltage (V)');
  else
    ylabel('Amplitude');
  end
  box('off'); grid('off');

  % Figure 4: Overlay of model template (green) and mean event (red)
  h4 = figure(4);set(h4,'visible',showfig);set(h4,'OuterPosition',[0 0 fw fh]);
  addToolbarExplorationButtons(gcf);
  if ~isempty(regexpi(vector,figform))
    try
      reduce_plot(t_fit,y_avg,'-g','linewidth',3);
    catch
      plot(t_fit,y_avg,'-g','linewidth',3);
    end
  else
    plot(t_fit,y_avg,'-g','linewidth',3);
  end
  ax = gca;
  ax.Toolbar.Visible = 'off';
  hold on;
  if ~isempty(regexpi(vector,figform))
    try
      reduce_plot(t_fit,mean(y_template,2),'-r','linewidth',2);
    catch
      plot(t_fit,mean(y_template,2),'-r','linewidth',2);
    end
  else
    plot(t_fit,mean(y_template,2),'-r','linewidth',2);
  end
  hold off;
  xlim([t_fit(1) t_fit(end)]);
  ylim([min(y_avg)-range(y_avg)*0.1 max(y_avg)+range(y_avg)*0.1]);
  xlabel('Time (s)');
  if strcmp(yunit,'A')
    ylabel('Current (A)');
  elseif strcmp(yunit,'V')
    ylabel('Voltage (V)');
  else
    ylabel('Amplitude');
  end
  grid('off'); box('off');
  title('Template and Mean Event Overlay')

  % Calculate axes autoscaling for figure 5
  y_autoscale = 0.1*range(Trace);
  y_maxlim = max(Trace)+y_autoscale;
  y_minlim = min(Trace)-y_autoscale;
  % Figure 5: Event wave overlaid with template fits
  h5 = figure(5);set(h5,'visible',showfig);set(h5,'OuterPosition',[fw 0 fw fh]);
  addToolbarExplorationButtons(gcf);
  ylimits = [y_minlim y_maxlim]; % Encoded y-axis autoscaling
  if ~isempty(regexpi(vector,figform))
    try
      reduce_plot(t,Trace,'-','color',[0.85,0.85,0.85],'linewidth',1.25);
    catch
      plot(t,Trace,'-','color',[0.85,0.85,0.85]);
    end
  else
    plot(t,Trace,'-','color',[0.85,0.85,0.85]);
  end
  ax = gca;
  ax.Toolbar.Visible = 'off';
  hold on;
  if ~isempty(regexpi(vector,figform))
    try
      reduce_plot(t,FITtrace,'-r');
    catch
      plot(t,FITtrace,'-r');
    end
  else
    plot(t,FITtrace,'-r');
  end
  hold off;
  xlim([min(t),max(t)]);
  ylim([ylimits(1), ylimits(2)]);
  grid('off'); box('off');
  title('Wave and events (after applying event criterion)');
  xlabel('Time (s)');
  if strcmp(yunit,'A')
    ylabel('Current (A)');
  elseif strcmp(yunit,'V')
    ylabel('Voltage (V)');
  else
    ylabel('Amplitude');
  end

  % Calculate statistics
  ET = Event_time;
  IEI = [NaN; diff(ET)];
  if exmode==1
    % Do nothing
  elseif exmode==2
    % Convert interevent intervals corresponding to events spanning exclusion zones to NaN
    if ~isempty(excl)
      for i=1:size(excl,1)
        temp = find(ET>=excl(i,2));
        if ~isempty(temp)
          IEI(temp(1)) = NaN;
        end
      end
    end
  else
    error('exclusion method not recognised')
  end
  Amplitude = mean(SCALE);

  % Print basic information
  %fprintf(['--------------------------EVENTER---------------------------\n',...
  %         ' Automatic PSC/PSP detection using FFT-based deconvolution\n',...
  %         ' and event peak analysis by least-squares template fitting\n',...
  %         ' Version v1.0 Copyright 2014 Andrew Charles Penn                \n\n']);
  fid = fopen(filewave,'w+t');
  fprintf(fid,'Filename: %s\n',filename);
  fprintf(fid,'Channel: %d\n',channel);
  fprintf(fid,'Wave number: %d\n',wave-1);
  fprintf(fid,'Wave name: %s\n',names{wave});
  fprintf(fid,'Total number of events detected: %d\n',n);
  fprintf(fid,'Duration of recording analysed (in s): %.1f\n',AnalysedTime);
  if strcmp(yunit,'A')
    fprintf(fid,'Mean event amplitude (pA): %.3g\n',Amplitude*1e+12);
  elseif strcmp(yunit,'A')
    fprintf(fid,'Mean event amplitude (mV): %.3g\n',Amplitude*1e+3);
  else
    fprintf(fid,'Mean event amplitude: %.3g\n',Amplitude);
  end
  fprintf(fid,'Event frequency (in Hz): %.3g\n',n/AnalysedTime);
  fprintf(fid,'High-pass filter cut-off on deconvoluted wave (at -3 dB, in Hz): %.3g\n',hpf);
  fprintf(fid,'Low-pass filter cut-off on deconvoluted wave (at -3 dB, in Hz): %.3g\n',lpf);
  fprintf(fid,'Vector of model template time constants (in ms): [%.3g,%.3g]\n',TC*1e3);
  fprintf(fid,'Standard deviation of the noise of the deconvoluted wave (a.u.): %.3g\n',noiseSD);
  fprintf(fid,'Scale factor of noise standard deviations for threshold setting: %.3g\n',SF);
  fprintf(fid,'Theoretical false positive detection rate before applying event criterion (in Hz): %.3g\n',(1-normcdf(SF))*sample_rate);
  fprintf(fid,'Sign of the event peaks: %s\n',s);
  fprintf(fid,'Criterion used for event screening: %s\n',criterion);
  fprintf(fid,'Maximum time after peak (expressed in time constants) used for template fit: %.3g\n',taus);
  fprintf(fid,'Dead time from event start required for template fit (ms): %.3g\n',time2peak*1e3);
  fprintf(fid,'Minimum acceptable correlation coefficient for the template fit: %.3g\n',rmin);
  fprintf(fid,'Episodic data window limits centred around each event: [%.3g,%.3g]\n',win);
  fprintf(fid,'Sample rate of the recording (in kHz): %d\n',sample_rate*1e-03);
  fprintf(fid,'lsqfit exitflag for fitting the noise peak: %d\n',exitflag);
  fprintf(fid,'Exclusion zones:\n');
  for i=1:size(excl,1)
      fprintf(fid,'%.6f\t%.6f\n',[excl(i,1) excl(i,2)]);
  end
  fclose(fid);

  % Save episodic data waves to file
  event_data = cat(2,sample_rate^-1*[0:size(y_events,1)-1]',y_events);
  ensemble_mean = [sample_rate^-1*[0:size(y_events,1)-1]',mean(y_events,2)];
  if exist('eventer.output','dir')==0
    mkdir('eventer.output');
  end
  cd eventer.output
  if exist(filewave,'dir')==0
    mkdir(filewave);
  end
  chdir(filewave);
  %ephysIO('event_data.phy',event_data,xunit,yunit);
  %ephysIO('ensemble_mean.phy',ensemble_mean,xunit,yunit);
  %if ~strcmpi(export,'none')
  % ephysIO(sprintf('event_data.%s',export),event_data,xunit,yunit);
  % ephysIO(sprintf('ensemble_mean.%s',export),ensemble_mean,xunit,yunit);
  %end
  ephysIO(sprintf('event_data.%s',export),event_data,xunit,yunit);
  ephysIO(sprintf('ensemble_mean.%s',export),ensemble_mean,xunit,yunit);

  % Save detection parameters and event amplitudes and times to text file
  dlmwrite('_parameters',[TC';noiseSD;AnalysedTime;sample_rate;npeaks],'newline','pc')
  dlmwrite('_offset',win(1),'newline','pc');
  if exist('summary.txt','file')~=0
    delete('summary.txt');
  end
  movefile(['../../' filewave],'summary.txt');
  if exist('txt','dir')==0
    mkdir('txt');
  end
  chdir('txt');
  dlmwrite('times.txt',Event_time,'delimiter','\t','newline','pc');
  peak = SCALE';
  dlmwrite('peak.txt',SCALE','delimiter','\t','newline','pc');
  dlmwrite('IEI.txt',IEI,'delimiter','\t','newline','pc');
  dlmwrite('features.txt',features,'delimiter','\t','newline','pc');
  cd ..

  % Save figures and images
  if exist('img','dir')==0
    mkdir('img');
  end
  chdir('img');
  len = numel(figform);
  if strcmpi(figform,'fig')
    set(h1, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    set(h2, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    set(h3, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    set(h4, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    set(h5, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    savefig(h1,'output_histogram.fig','compact');
    savefig(h2,'output_decon.fig','compact');
    savefig(h3,'output_avg.fig','compact');
    savefig(h4,'output_template.fig','compact');
    savefig(h5,'output_wave.fig','compact');
  elseif ~isempty(regexpi(bitmap,figform))
    print(h1,sprintf('output_histogram.%s',figform(1:min(3,len))),sprintf('-d%s',figform),'-r300','-opengl');
    print(h2,sprintf('output_decon.%s',figform(1:min(3,len))),sprintf('-d%s',figform),'-r300','-opengl');
    print(h3,sprintf('output_avg.%s',figform(1:min(3,len))),sprintf('-d%s',figform),'-r300','-opengl');
    print(h4,sprintf('output_template.%s',figform(1:min(3,len))),sprintf('-d%s',figform),'-r300','-opengl');
    print(h5,sprintf('output_wave.%s',figform(1:min(3,len))),sprintf('-d%s',figform),'-r300','-opengl');
    if strcmpi(figform,'jpeg')
      movefile('output_histogram.jpe','output_histogram.jpg');
      movefile('output_decon.jpe','output_decon.jpg');
      movefile('output_avg.jpe','output_avg.jpg');
      movefile('output_template.jpe','output_template.jpg');
      movefile('output_wave.jpe','output_wave.jpg');
    end
  elseif ~isempty(regexpi(vector,figform))
    print(h1,sprintf('output_histogram.%s',figform(1:min(3,len))),sprintf('-d%s',figform),'-painters');
    print(h2,sprintf('output_decon.%s',figform(1:min(3,len))),sprintf('-d%s',figform),'-painters');
    print(h3,sprintf('output_avg.%s',figform(1:min(3,len))),sprintf('-d%s',figform),'-painters');
    print(h4,sprintf('output_template.%s',figform(1:min(3,len))),sprintf('-d%s',figform),'-painters');
    print(h5,sprintf('output_wave.%s',figform(1:min(3,len))),sprintf('-d%s',figform),'-painters');
    if strcmpi(figform(1:2),'ps')
      movefile('output_histogram.ps*','output_histogram.ps');
      movefile('output_decon.ps*','output_decon.ps');
      movefile('output_avg.ps*','output_avg.ps');
      movefile('output_template.ps*','output_template.ps');
      movefile('output_wave.ps*','output_wave.ps');
    elseif strcmpi(figform,'meta')
      movefile('output_histogram.met','output_histogram.emf');
      movefile('output_decon.met','output_decon.emf');
      movefile('output_avg.met','output_avg.emf');
      movefile('output_template.met','output_template.emf');
      movefile('output_wave.met','output_wave.emf');
    end
  end
  cd ../../..

  if merge == 1
    merge_data(average,s,win,export,optimoptions,cwd,figform,config,taus);
  end

end

%----------------------------------------------------------------------

function [tpeak, t50, tdecay, auc, skew] = shape (t, y, offset, Fs, s, base)

  % Calculate quick and dirty estimate of full-width half-maximum (t50)
  % and time-to-peak (tpeak) for the event

  % Calculate data size
  n = numel(y);

  % Subtract baseline
  y = y-offset;
  if s=='-'
    y = y * -1;
  end

  % Approximate FWHM by calculating the duration spent at > half-maximal amplitude
  y50 = 0.5 * max(y);
  t50 = sum(y>y50)/Fs;

  % Approximate time-to-peak by finding the point of maximum deviation from the baseline
  ypeak = max(y);
  ipeak = find(y==ypeak);
  tpeak = (sum(t<mean(t(ipeak)))/Fs) - base;

  % Approximate decay half-life by calculating the duration of the event post-peak spent at > half-maximal amplitude
  tdecay = sum((y>y50) & (t>mean(t(ipeak))))/Fs;

  % Calculate the area-under-curve
  auc = trapz(t,y);

  % Calculate event distribution skew (with respect to the baseline i.e. 0)
  skew = (n^-1)*sum(y.^3) / (sqrt((n^-1)*sum(y.^2)))^3;

end

%----------------------------------------------------------------------

function merge_data(average,s,win,export,optimoptions,cwd,figform,config,taus)

  % Merge the event data for all the waves analysed in this experiment

  % Merge all data from eventer.output folder
  count = 0;
  cd eventer.output
  numtraces=0;
  dirlist = dir('.');
  dirlist = char(dirlist.name);
  numdir = size(dirlist,1);
  dirarray = mat2cell(dirlist,ones(numdir,1),size(dirlist,2));
  t = [];
  root = pwd;
  for i=1:numdir
    dirname{i} = strtrim(dirarray(i));
    if ~strcmp(dirname{i},'ALL_events') &&... %isdir(char(dirname{i})) && was removed
    ~strcmp(dirname{i},'.') &&...
    ~strcmp(dirname{i},'..') &&...
    ~strcmp(dirname{i},'.DS_Store')
      count = count+1;
	  if i==3
	    chdir(char(dirname{i}));
	  elseif i>3
        %chdir(['../',dirname{i}]);
	    location = strcat(root,filesep,dirname{i});
	    chdir(char(location));
	  end
      %if exist('event_data.phy','file')
      %  [temp,xdiff,xunit,yunit] = ephysIO('event_data.phy');
      if exist(sprintf('event_data.%s',export),'file')
        [temp,xdiff,xunit,yunit] = ephysIO(sprintf('event_data.%s',export));
        temp(:,1) = temp(:,1)+win(1);
        temp(abs(temp(:,1))<eps,1)=0;
        if isempty(t)
          t = temp(:,1);
        elseif count>1
          if numel(t)==numel(temp(:,1))
             if all(t==temp(:,1))~=1
               error('Inconsistent window dimensions. Cannot merge data.');
             end
          else
            error('Inconsistent window dimensions. Cannot merge data.');
          end
        end
        temp(:,1) = [];
        data{count} = temp;
        numtraces(count,1) = size(data{count},2);
      else
        numtraces(count,1) = 0;
      end
      parameters(:,count) = load('-ascii','_parameters');
      TC(count,:) = parameters(1:2,count)';
      sigma(count,1) = parameters(3,count);
      AnalysedTime(count,1) = parameters(4,count);
      sample_rate(count,1) = parameters(5,count);
      if sample_rate(count)~=sample_rate(1)
        error('Inconsistent sample rate. Cannot merge data.')
      end
      if exist('txt','dir')
        cd txt
        IEI{i,1} = load('-ascii','IEI.txt');
        peak{i,1} = load('-ascii','peak.txt');
        if exist('features.txt','file')
          features{i,1} = load('-ascii','features.txt');
        end
      end
    else
      % The name in the directory list is not a folder suitable for the
      % merge process so do nothing
    end
	%chdir(root);
  end
  chdir(root);
  % Overide data units with configurationn argument
  if strcmpi(config,'CC')
    yunit = 'V';
  elseif strcmp(config,'VC')
    yunit = 'A';
  elseif strcmp(config,'')
    yunit = '';
  end
  numEvents = sum(numtraces);
  if exist('ALL_events','dir')==0
    mkdir('ALL_events');
  end
  cd('ALL_events')
  if numEvents > 0
    tau_rise = sum(TC(:,1).*numtraces/sum(numtraces));
    tau_decay = sum(TC(:,2).*numtraces/sum(numtraces));
    y = cell2mat(data);
    IEI = cell2mat(IEI);
    peak = cell2mat(peak);
    if exist('features','var')
      features = cell2mat(features);
    end
    %nanidx = isnan(IEI);
    %peak(nanidx) = [];
    %IEI(nanidx) = [];
    %y(:,nanidx) = [];
    freq = 1/median(IEI);
    numEvents = numel(peak);
    events = [sample_rate(1)^-1*[0:size(y,1)-1]',y];
    if strcmp(average,'mean')
      y_avg = mean(y,2);
    elseif strcmp(average,'median')
      y_avg = median(y,2);
    end
    ensemble_average = [sample_rate(1)^-1*[0:size(y_avg,1)-1]',y_avg];

    % Fit a sum of exponentials function to the ensemble average event
    p0 = [1,tau_rise,tau_decay];
    tpeak0 = p0(3)*p0(2)/(p0(3)-p0(2))*log(p0(3)/p0(2));
    tlim = p0(3)*taus;
    idx = t>=0 & t<=tpeak0+tlim;
    tdata = t(idx);
    fun2 = @(p,tdata)-exp(-tdata/p(1))+exp(-tdata/p(2));
    ypeak0 = fun2([p0(2),p0(3)],tpeak0);
    if s=='-'
      NF = ypeak0/min(y_avg(idx));  % Normalization factor
      ydata = y_avg(idx)*NF;
    elseif s=='+'
      NF = ypeak0/max(y_avg(idx));  % Normalization factor
      ydata = y_avg(idx)*NF;
    end
    fun3 = @(p,tdata)p(1)*(-exp(-tdata/p(2))+exp(-tdata/p(3)));
    try
      [p,resnorm,residual,exitflag] = lsqfit(fun3,p0,tdata,ydata,[],[],optimoptions);
      tpeak = p(3)*p(2)/(p(3)-p(2))*log(p(3)/p(2));
      fitAmplitude = abs((fun3(p,tpeak))/NF);
      fitIntegral = abs(p(1)*(p(3)-p(2))/NF);
      fit = fun3(p,tdata)/NF;
      residuals = y_avg(idx)-fit;
      errflag = 0;
    catch
      % do nothing
      errflag = 1;
    end

    % Get screen size and set figure sizes
    set(0,'units','pixels');
    scn_sz = get(0,'screensize');
    sw = scn_sz(3);
    sh = scn_sz(4);
    fw = sw/3;
    fh = sh/2;

    % Figures
    vector = 'pdf eps epsc eps2 epsc2 meta svg ps psc ps2 psc2';
    bitmap = 'jpeg png tiff tiffn bmpmono bmp bmp16m bmp256 hdf pbm pbmraw pcxmono pcx24b pcx256 pcx16 pgm pgmraw ppm ppmraw';

    % Calculate axes autoscaling for figure 6
    y_autoscale = 0.1*max(range(y));
    y_maxlim = max(max(y))+y_autoscale;
    y_minlim = min(min(y))-y_autoscale;
    % Figure 6: Identified events aligned and overlayed with the ensemble average event
    h6 = figure(6);set(h6,'OuterPosition',[fw*2 0 fw fh]);
    addToolbarExplorationButtons(gcf);
    ylimits = [y_minlim y_maxlim]; % Encoded y-axis autoscaling
    if ~isempty(regexpi(vector,figform))
      try
        reduce_plot(t,y,'color',[0.85,0.85,0.85]);
      catch
        plot(t,y,'color',[0.85,0.85,0.85]);
      end
    else
      plot(t,y,'color',[0.85,0.85,0.85]);
    end
    ax = gca;
    ax.Toolbar.Visible = 'off';
    hold on;
    if ~isempty(regexpi(vector,figform))
      try
        reduce_plot(t,y_avg,'-b');
        if errflag < 1
          reduce_plot(tdata,fit,'r-','linewidth',2);
        end
      catch
        plot(t,y_avg,'-b');
        if errflag < 1
          plot(tdata,fit,'r-','linewidth',2);
        end
      end
    else
      plot(t,y_avg,'-b');
      if errflag < 1
        plot(tdata,fit,'r-','linewidth',2);
      end
    end
    xlim(win);
    ylim([ylimits(1), ylimits(2)]);
    if strcmp(average,'mean')
      title('Ensemble mean (blue) and fit (red)');
    elseif strcmp(average,'median')
      title('Ensemble median (blue) and fit (red)');
    end
    xlabel('Time (s)');
    if strcmp(yunit,'A')
      ylabel('Current (A)');
    elseif strcmp(yunit,'V')
      ylabel('Voltage (V)');
    else
      ylabel('Amplitude');
    end
    box('off'); grid('off');
    hold off;
    if errflag > 0
       error('Fit to ensemble average event failed');
    end

    % Save data
    save('ensemble_average.txt','ensemble_average','-ascii','-tabs');
    save('fit.txt','fit','-ascii','-tabs');
    save('residuals.txt','residuals','-ascii','-tabs');
    %fid = fopen('event_counts.txt','w');
    %fprintf(fid,'%d\n',numtraces);
    %fclose(fid);
    %ephysIO('event_data.phy',events,xunit,yunit);
    %ephysIO('ensemble_average.phy',ensemble_average,xunit,yunit);
    %if ~strcmpi(export,'none')
    %  ephysIO(sprintf('event_data.%s',export),events,xunit,yunit);
    %  ephysIO(sprintf('ensemble_average.%s',export),ensemble_average,xunit,yunit);
    %end
    ephysIO(sprintf('event_data.%s',export),events,xunit,yunit);
    ephysIO(sprintf('ensemble_average.%s',export),ensemble_average,xunit,yunit);   
    dlmwrite('_parameters',[tau_rise;tau_decay;sigma;AnalysedTime],'newline','pc');
    dlmwrite('_offset',win(1),'newline','pc');
    if exist('img','dir')==0
      mkdir('img');
    end
    cd('img')
    len = numel(figform);
    if strcmpi(figform,'fig')
      set(h6, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
      savefig(h6,'output_avg.fig','compact')
    elseif ~isempty(regexpi(bitmap,figform))
      print(h6,sprintf('output_avg.%s',figform(1:min(3,len))),sprintf('-d%s',figform),'-r300','-opengl');
      if strcmpi(figform,'jpeg')
        movefile('output_avg.jpe','output_avg.jpg');
      end
    elseif ~isempty(regexpi(vector,figform))
      print(h6,sprintf('output_avg.%s',figform(1:min(3,len))),sprintf('-d%s',figform),'-painters');
      if strcmpi(figform(1:2),'ps')
        movefile('output_avg.ps*','output_avg.ps');
      end
      if strcmpi(figform,'meta')
        movefile('output_avg.met','output_avg.emf');
      end
    end
    cd ..
    if exist('txt','dir')==0
      mkdir('txt');
    end
    cd('txt')
    dlmwrite('IEI.txt',IEI,'delimiter','\t','newline','pc');
    dlmwrite('peak.txt',peak,'delimiter','\t','newline','pc');
    if exist('features','var')
      dlmwrite('features.txt',features,'delimiter','\t','newline','pc');
    end
    cd ..
  end
  sigma = sum(sigma.*AnalysedTime/sum(AnalysedTime));
  AnalysedTime = sum(AnalysedTime);

  % Print basic information
  fid=fopen('ALL_events','w+t');
  fprintf(fid,'Number of waves analysed: %d\n',count);
  fprintf(fid,'Total recording time analysed (in s): %.1f\n',AnalysedTime);
  fprintf(fid,'Total number of events: %d\n',numEvents);
  fprintf(fid,'Event frequency (in Hz): %.3g\n',numEvents/AnalysedTime);
  if numEvents > 0
    fprintf(fid,'Ensemble average: %s\n',average);
    if strcmp(yunit,'A')
      fprintf(fid,'Amplitude of the model PSC fit (pA): %.3g\n',fitAmplitude*1e+12);
      fprintf(fid,'Integral (charge) of the model PSC fit (fC): %.4g\n',fitIntegral*1e+15);
      fprintf(fid,'Rise time constant of the model PSC fit (ms): %.3g\n',p(2)*1e+03);
      fprintf(fid,'Decay time constant of the model PSC fit (ms): %.3g\n',p(3)*1e+03);
    elseif strcmp(yunit,'V')
      fprintf(fid,'Amplitude of the model PSP fit (mV): %.3g\n',fitAmplitude*1e+03);
      fprintf(fid,'Integral of the model PSP fit (mV.ms): %.4g\n',fitIntegral*1e+06);
      fprintf(fid,'Rise time constant of the model PSP fit (ms): %.3g\n',p(2)*1e+03);
      fprintf(fid,'Decay time constant of the model PSP fit (ms): %.3g\n',p(3)*1e+03);
    else
      fprintf(fid,'Amplitude of the model event fit: %.3g\n',fitAmplitude);
      fprintf(fid,'Integral of the model event fit: %.4g\n',fitIntegral);
      fprintf(fid,'Rise time constant of the model event fit (ms): %.3g\n',p(2)*1e+03);
      fprintf(fid,'Decay time constant of the model event fit (ms): %.3g\n',p(3)*1e+03);
    end
      fprintf(fid,'lsqfit exitflag for fitting the ensemble average event: %d\n',exitflag);
  else
    fprintf(fid,'Note: Not enough events for analysis\n');
  end
  fprintf(fid,'Standard deviation of the noise of the deconvoluted waves (a.u.): %.3g\n',sigma);
  fclose(fid);
  movefile(['./ALL_events'],'summary.txt');
  cd ../..
  chdir(cwd);

end


%--------------------------------------------------------------------------

function [dydx, y, x] = ndiff (y, x)

if nargin ~= 2
 error('Invalid number of input arguments');
end

if all(size(x) == 1) || ~any(size(x) == 1) || length(x) ~= length(y)
 error('x and y must be vectors of the same size');
end

% Assess sampling characteristics of input with precision of 10e-9
isDiscrete=~any(round(diff(x)*10e9)-mean(round(diff(x)*10e9)));
 if isDiscrete == 0
  warning('non-discrete','Input must consist of data sampled at evenly spaced points');
 end

% Set all input vectors as column vectors where applicable
x=x(:); y=y(:);

% Create matrices for numerical differentiation
% These have column length l-2 since complete numerical differentiation cannot be calulated for values within 1 point from the start and end
l=length(x);
my=ones(l-2,3);mx=ones(l-2,3);
 for i=1:3
  my(:,i)=y(i:l-(3-i));
  mx(:,i)=x(i:l-(3-i));
 end

% Calculate first derivative
dx=mx(:,3)-mx(:,1);
dy=my(:,3)-my(:,1);
dydx=dy./dx;
x=mx(:,2);
y=my(:,2);

end

%--------------------------------------------------------------------------

function [yf, x] = medianf (y, x, r)
   

if nargin ~= 3
 error('Invalid number of input arguments');
end

if all(size(x) == 1) || ~any(size(x) == 1) || any(size(x)~=size(y))
 error('x and y must be vectors of the same size');
end

if isinf(r) || ~all(size(r) == 1) || r<=0 || r~=round(r)
 error('r must be a nonnegative integer');
end

% Set all input vectors as column vectors and calculate vector size
x=x(:); y=y(:);
m = length(y);

% Calculate total number of y-points to average within sliding box
P=2*r+1;

% Treatment of the ends of the data
y = bounce(y,r);

% Implement median smoothing filter
if r < 250

  % Fast algorithm for small filter rank
  l = length(y);
  Y = zeros(l-(P-1),P);
  for i=1:P
   Y(:,i)=y(i:l-(P-i));
  end
  Y=sort(Y,2);
  yf=Y(:,r+1);

else

  % Fast algorithm for large filter rank
  % Calculate running mean and variance. Algorithm from:
  %  J.E. Hadstate (2008) Efficient Moving Average and Moving Variance Calculations
  %  https://www.dsprelated.com/showthread/comp.dsp/97276-1.php
  M = zeros(m,1);
  mu = M;
  V = zeros(m,1);
  v  = V;
  SX1 = sum(y(1:P));
  SX2 = sum(y(1:P).^2);
  X1 = 0;
  X2 = 0;
  Y1 = 0;
  Y2 = 0;
  M(1) = SX1/P;
  V(1) = (P*SX2-(SX1*SX1))/(P*(P-1));
  for k=2:m
    Y1 = y(k-1);
    Y2 = Y1^2;
    X1 = y(P+k-1);
    X2 = X1^2;
    SX1 = SX1 + X1 - Y1;
    SX2 = SX2 + X2 - Y2;
    M(k) = SX1/P;
    V(k) = (P*SX2-(SX1*SX1))/(P*(P-1));
  end

  % Use the mean and standard deviation of the window centered around first point to
  % choose initial bin edges then bin the data from that window. Note that the ends
  % of the data have already been padded with r points using the bounce function
  B = 1000;
  edges = [-inf,linspace(min(M(1)-sqrt(V(1))),max(M(1)+sqrt(V(1))),B),inf];
  [N] =  histc(y(1:P),edges);

  % Compute the running median based on use of Tibshirani's binapprox
  % algorithm with the update problem. This algorithm scales extremely
  % well with increasing filter rank.
  yf = zeros(m,1);
  a = sum(edges <= y(1));
  j = 0;
  for i=1:m
    exitflag = 0;
    while exitflag < 1
      % Find bin containing the median
      j = 0;
      n = 0;
      while n <= r+1
        j = j+1;
        n = n + N(j);
      end
      % If the median lies outside of the existing bins, calculate new bin edges
      if any(j==[1 B+1])
        edges = [-inf,linspace(min(M(i)-sqrt(V(i))),max(M(i)+sqrt(V(i))),B),inf];
        [N] =  histc(y(i:P+i-1),edges);
        exitflag = 0;
      else
        exitflag = 1;
      end
    end
    % Allocate the bin center to yf
    yf(i) = (edges(j)+edges(j+1))/2;
    % Update bin counts after sliding the window
    N(a) = N(a) - 1;
    a = sum(edges <= y(i+1));
    if i < m
      b = sum(edges <= y(P+i));
      N(b) = N(b) + 1;
    end
  end

end


end

%--------------------------------------------------------------------------

% readMeta.m from Luka Campagnola
function f = readMeta(file)
info = hdf5info(file);
f = readMetaRecursive(info.GroupHierarchy.Groups(1));
end


function f = readMetaRecursive(root)
typ = 0;
for i = 1:length(root.Attributes)
    if strcmp(root.Attributes(i).Shortname, '_metaType_')
        typ = root.Attributes(i).Value.Data;
        break
    end
end
if typ == 0
    printf('group has no _metaType_')
    typ = 'dict';
end

list = 0;
if strcmp(typ, 'list') || strcmp(typ, 'tuple')
    data = {};
    list = 1;
elseif strcmp(typ, 'dict')
    data = struct();
else
    printf('Unrecognized meta type %s', typ);
    data = struct();
end

for i = 1:length(root.Attributes)
    name = root.Attributes(i).Shortname;
    if strcmp(name, '_metaType_')
        continue
    end
    val = root.Attributes(i).Value;
    if isa(val, 'hdf5.h5string')
        val = val.Data;
    end
    if list
        ind = str2num(name)+1;
        data{ind} = val;
    else
        data.(name) = val;
    end
end

for i = 1:length(root.Datasets)
    fullName = root.Datasets(i).Name;
    name = stripName(fullName);
    file = root.Datasets(i).Filename;
    data2 = hdf5read(file, fullName);
    if list
        ind = str2num(name)+1;
        data{ind} = data2;
    else
        data.(name) = data2;
    end
end

for i = 1:length(root.Groups)
    name = stripName(root.Groups(i).Name);
    data2 = readMetaRecursive(root.Groups(i));
    if list
        ind = str2num(name)+1;
        data{ind} = data2;
    else
        data.(name) = data2;
    end
end
f = data;
return;
end


function f = stripName(str)
inds = strfind(str, '/');
if isempty(inds)
    f = str;
else
    f = str(inds(length(inds))+1:length(str));
end
end

%--------------------------------------------------------------------------

%     Function File: [y, k] = bounce (y, n)
%
%     Extend each end of the y-vector by n number of points for
%     'bounce' (or mirror) end-effect correction. The output
%     includes the extended y-vector and the number of bounces (k).
%
%     bounce v1.0 (last updated: 16/09/2011)
%     Author: Andrew Charles Penn
%     https://www.researchgate.net/profile/Andrew_Penn/


function [y] = bounce (y, n)

if nargin ~= 2
 error('Invalid number of input arguments');
end

if all(size(y) == 1) || ~any(size(y) == 1)
 error('y must be a vector');
end

if isinf(n) || ~all(size(n) == 1) || n<=0 || n~=round(n)
 error('n must be a nonnegative integer');
end


% Set input vector as column vector if applicable
y=y(:); y_ref=y;

% Extend each end by n number of points for 'bounce' end-effect correction
N=numel(y);
k=0;
while (n > N)
 if (k == 0) || (k/2 == round(k/2))
  y=cat(1,flipud(y_ref),y,flipud(y_ref));
 else
  y=cat(1,y_ref,y,y_ref);
 end
 k=k+1;
 n=n-N;
end
if (k == 0) || (k/2 == round(k/2))
 y=cat(1,flipud(y_ref(1:n)),y,flipud(y_ref(end-n+1:end)));
else
 y=cat(1,y_ref(end-n+1:end),y,flipud(y_ref(1:n)));
end

end
%--------------------------------------------------------------------------
%     Function File: filter1
%
%         Usage: YF = filter1 (Y, t, HPF, LPF, method)
%
%            OR: filter1 (filename, ending, HPF, LPF, method, '-file')
%
%     This function combines high- and/or low-pass 1-D filtering at the
%     set -3 dB cut-off frequencies, which are defined in the parameters
%     HPF and LPF (units Hertz). To switch off the high-pass filter, enter
%     an HPF value of 0. To switch off the low-pass filter, enter an LPF
%     value of Inf. This function avoids end effects by using a bounce
%     algorithm.
%
%     Low pass filtering is achieved using a digital Gaussian filter at the
%     specified cut-off frequency LPF. 

%     High pass filtering is achieved using a digital binomial filter
%     (default method='binomial'). If the method is specified as 'median', 
%     then a median filtering method is used at the specified HPF cut-off
%     frequency and the resulting trace is subtracted from the input.
%     The filter rank is estimated from the desired cut-off value for the
%     analagous linear filter (boxcar). The median filter does not cause
%     the edge effects that linear filters are prone to, but is slower
%     Note that in the case of the median filter, the HPF value reported
%     is an estimate of the -3 dB cut-off for the initial low-pass 
%     filtering prior to subtraction and not the -3 dB cutoff of the 
%     resulting high-pass filter.
%
%     In the first example of the function usage, filtering is by default
%     applied to each data column in Y as a function of time (t) defined
%     by the first and second input arguments respectively. The output
%     arguments are the filtered data values (YF).
%
%     The second example of the function usage is defined by setting a
%     sixth input argument to '-file'. In this mode, the function instead
%     loads the data from the text file named in the first input argument
%     and saves the processed data with the filename appended with the
%     ending given in the second input argument. If the ending option is
%     left empty ('[]'), the function will overwrite the original file.
%     The used cut-off values are saved in a separate file with the new
%     filename appended with '_Fc'. The data is saved in the ephysIO hdf5-
%     based MATLAB format with the .mat filename extension.
%
%     When used to band-pass filter, this function performs operations in
%     the order: 1) High-pass, then 2) Low-pass.
%
%     This function requires the following functions and their dependencies:
%     'ephysIO', 'hpfilter', 'lpfilter', 'binomialf' and 'medianf'.
%
%     Bibliography:
%     Marchand & Marmet (1983) Rev Sci Instrum 54(8): 1034-1041
%     Moore & Jorgenson (1993) Anal Chem 65: 188-191
%
%     filter1 v1.3 (last updated: 21/02/2017)
%     Author: Andrew Charles Penn
%     https://www.researchgate.net/profile/Andrew_Penn/


function YF = filter1 (argin1, argin2, HPF, LPF, method, option)

  if nargin<4
    error('Invalid number of input arguments');
  end

  if isinf(HPF) || ~all(size(HPF) == 1) || HPF<0
    error('If non-zero, the HPF cut-off must be a nonnegative, finite value in unit Hz');
  end

  if ~all(size(LPF) == 1) || LPF<=0
    error('If finite, the LPF cut-off must be a nonnegative, non-zero value in unit Hz');
  end

  if nargin>4
    if isempty(method)
      method='binomial';
    end
  else
    method='binomial';
  end

  if nargin>5
    if ~strcmp(option,'-file')
      error('If option is specified, it must be set to '-file'');
    end
  else
    option=[];
  end

  % Load data
  if strcmp(option,'-file')
    cwd = pwd;
    if ~isempty(regexpi(argin1(end-2:end),'.gz'))
      [pathstr,filename,ext] = fileparts(argin1(end-2:end));
    elseif ~isempty(regexpi(argin1(end-3:end),'.zip'))
      [pathstr,filename,ext] = fileparts(argin1(end-3:end));
    else
      [pathstr,filename,ext] = fileparts(argin1);
    end
    if ~isempty(pathstr)
      chdir(pathstr);
    end
    [data,xdiff,xunit,yunit,names,notes] = ephysIO (strcat(filename,ext));
    t=data(:,1);
    Y=data; Y(:,1)=[];
  else
    Y=argin1;
    t=argin2;
    if any(diff(diff(t)) > 1.192093e-07)
      % Variable sampling interval
      xdiff = 0;
    else
      xdiff = t(2)-t(1);
    end
  end
  l=length(t);

  % Perform linear interpolation on non-evenly spaced datasets
  if xdiff == 0
    disp(sprintf('Input must consist of data sampled at evenly spaced time points.\nData will undergo linear interpolation.'));
    sample_rate=(l-1)/(max(t)-min(t));
    tl=linspace(min(t),max(t),l);
    tl=tl(:);
    for i=1:size(Y,2)
      state = warning('query');
      warning off %#ok<WNOFF>
      Yl(:,i)=interp1q(t,Y(:,i),tl);
      warning(state)
    end
  elseif xdiff > 0
    sample_rate=(l-1)/(max(t)-min(t));
    tl=t(:);
    Yl=Y;
  end

% Filter data traces
if strcmp(option,'-file')
 figure(2)
 clf
end
for i=1:size(Yl,2)
  yf = Yl(:,i);
  if ~isinf(LPF)
    yf = gaussianf (yf, tl, LPF,'on');
  end
  yref = yf;
  if HPF > 0
    if strcmp(method,'median')
     % Analogy to a linear, boxcar filter:
     % -3 dB cut-off: Fc  ~ 0.443 / p * sample_rate,
     % where the number of points in the sliding window:
     % p = 2 * r + 1 and Fc is HPF and r is the filter rank
     p = ((0.443*sample_rate)/HPF);
     r=round((p-1)/2);
     HPF = 0.443 / (2*r+1) * sample_rate;
     arch = computer('arch');
     try
       if strcmpi(arch,'maci64')
         % Code to run on Mac 64-bit platform
         [ybase, tbase] = medianf_mex_maci64 (yref, tl, r);  %#ok<*ASGLU>
       elseif strcmpi(arch,'glnxa64')
         % Code to run on Linux 64-bit platform
         [ybase, tbase] = medianf_mex_glnxa64 (yref, tl, r);
       elseif strcmpi(arch,'win64')
         % Code to run on Windows 64-bit platform
         [ybase, tbase] = medianf_mex_win64 (yref, tl, r);
       end
     catch
       warning(sprintf(['A suitable MEX file for medianf is not available or failed to execute.\n',...
                        'Falling back to Matlab file']));
       [ybase, tbase] = medianf (yref, tl, r);
     end
     yf=yref-ybase;
     if strcmp(option,'-file')
      y_autoscale=0.05*(max(yref)-min(yref)); y_maxlim=max(yref)+y_autoscale; y_minlim=min(yref)-y_autoscale; % Encoded y-axis autoscaling
      figure(2); hold on; plot(tl,yref,'-','color',[0.75,0.75,0.75]); plot(tl,ybase,'k'); hold off; xlim([min(tl),max(tl)]); ylim([y_minlim y_maxlim]); box('off');
     end
    elseif strcmp(method, 'binomial')
     yf = hpfilter (yref, tl, HPF);
     if strcmp(option,'-file')
      y_autoscale=0.05*(max(yref)-min(yref)); y_maxlim=max(yref)+y_autoscale; y_minlim=min(yref)-y_autoscale; % Encoded y-axis autoscaling
      figure(2); hold on; plot(tl,yref,'-','color',[0.75,0.75,0.75]); plot(tl,yref-yf,'k-'); hold off; xlim([min(tl),max(tl)]); ylim([y_minlim y_maxlim]); box('off');
     end
   end
  elseif HPF == 0
   figure(2)
   close(2)
  end
 YF(:,i)=yf;
end

% Output for filter1 function in '-file' mode
if strcmp(option,'-file')
 diary('on');
 cutoffs=strcat(argin1,argin2,'_Fc.txt');
 if exist(cutoffs,'file') ~= 0
  delete(cutoffs);
 end
 method %#ok<NOPRT>
 format short g
 if strcmp(method,'binomial')
  HPF %#ok<NOPRT>
 elseif strcmp(method,'median')
  HPF %#ok<NOPRT>
  r %#ok<NOPRT>
 end
 LPF %#ok<NOPRT>
 diary('off');
 figure(1);
 y_autoscale=0.05*(max(max(YF))-min(min(YF))); y_maxlim=max(max(YF))+y_autoscale; y_minlim=min(min(YF))-y_autoscale; % Encoded y-axis autoscaling
 plot(tl,YF,'k-'); xlim([min(tl),max(tl)]); ylim([y_minlim y_maxlim]); box('off');
 ylabel(strcat('y-axis (',yunit,')'));
 xlabel(strcat('x-axis (',xunit,')'));
 output_data=cat(2,tl,YF); %#ok<NASGU>
 newfilename=strcat(filename,argin2,'.mat');
 ephysIO(newfilename,output_data,xunit,yunit,names,notes);
 if exist('filter1.output','dir') == 0
  mkdir('filter1.output');
 end
 cd filter1.output
 newfilename=strtok(newfilename,'.');
 if exist(newfilename,'dir') == 0
  mkdir(newfilename);
 end
 cd(newfilename);
 print(1,'output.png','-dpng');
 print(1,'output.eps','-depsc');
 if HPF > 0
  print(2,'baseline.png','-dpng');
  print(2,'baseline.eps','-depsc');
 end
 movefile('../../diary',cutoffs);
 cd ../..
 clear YF tl
end
end

%-------

%     Function File: [yf, tf, F] = gaussianf (y, t, Fc, correction)
%
%     Application of a Gaussian smoothing filter to the y-vector.
%     The cutoff of the low pass filter is at frequency Fc (in Hz)
%     Requires the Matlab Signal Processing toolbox.
%
%     Output vectors of the smoothed y-values and of the corresponding
%     time (t) values are returned.
%
%     gaussianf v1.0 (last updated: 16/09/2011)
%     Author: Andrew Charles Penn
%     https://www.researchgate.net/profile/Andrew_Penn/


function [yf, tf, F] = gaussianf (y, t, Fc, correction)

if nargin ~= 4
 error('Invalid number of input arguments');
end

if all(size(t) == 1) || ~any(size(t) == 1) || any(size(t)~=size(y))
 error('t and y must be vectors of the same size');
end

if isinf(Fc) || isnan(Fc) || ~all(size(Fc) == 1) || Fc<=0
 error('Fc must be a finite nonnegative number');
end

% Set all input vectors as column vectors
t=t(:); y=y(:);

% Calculate sampling frequency
N = numel(t);
Fs = (N-1)/(max(t)-min(t));

% Calculate the standard deviation of the Gaussian envelope in points
sigma = Fs*0.132505/Fc;

% Calculate Filter coefficients
if sigma < 0.62
  n = 1;
  b = sigma^2/2;
  F(1) = b;
  F(2) = 1-2*b;
  F(3) = b;
else
  p = 2*round(4*sigma)+1;
  n = (p-1)/2;
  b = -1/(2*sigma^2);
  i = (p-1)/2+1;
  F = zeros(1,p);
  for j=1:p
    F(j)=exp(b*(j-i)^2); 
  end
end

% Normalize coefficients
F = F/trapz(F);

% Bounce end-effect correction
if strcmp(correction,'on') == 1
 [y] = bounce(y,n);
 tf = t;
elseif strcmp(correction,'off') == 1
 tf = t(1+n:N-n);
end

% Gaussian smoothing filter
yf = filter(F,1,y);
yf(1:2*n)=[]; % The assymmetric point deletion compensates for the group delay
end

%----

function y = inrange(X,R,varargin)
%INRANGE tests if values are within a specified range (interval).
%   INRANGE(X,RANGE) tests if the values in X are within the range specified 
%   by RANGE.  X can be a vector or matrix.  
%
%   RANGE is a range in the form [LOW HIGH] against which each value in X will 
%   be tested.  RANGE can also be a two-column vector whose ith row is of form 
%   RANGE(i,:) = [LOW HIGH].  In this form, input X must be a vector with the
%   same length as RANGE, and each element of X is tested against the
%   range in the corresponding row of RANGE.
%
%   INRANGE(X,RANGE,BOUNDARY) specifies whether the endpoints of the specified 
%   range should be included or excluded from the interval.  The options for 
%   BOUNDARY are:
%
%       'includeboth' : Both end points included in interval (default)
%       'includeleft' : Left end point only included in interval
%       'includeright': Right end point only included in interval
%       'excludeboth' : Neither end point included in interval
%
%   If the LOW and HIGH values for RANGE are equal, that single value will be 
%   found in range only under the 'includeboth' option (default). Otherwise, 
%   for this case, no values will be found in range.
%
%   Examples:  
%      X = 1:10
%      X =
%          1     2     3     4     5     6     7     8     9    10
%
%      Y = inrange(X,[5 7.2],'includeboth')
%      Y =
%          0     0     0     0     1     1     1     0     0     0
%
%      Y = inrange(X,[5 7.2],'excludeboth')
%      Y =
%          0     0     0     0     0     1     1     0     0     0
%
%   See also ISVALID.
%Version 2: Per online comment from "Jos x" (thanks!), fixed the following:
% 1) no need for find, use logical indexing;  DONE
% 2) return a logical array (false/true) instead of a zero/one array; DONE
% 3) for a 1x2 input R you do no need the repmat, as you split it into two
%    scalars; DONE
% 4) reduce overhead: first get the size of X, than use X = X(:), and
%    finally reshape Y using the stored size of X; DONE
% 5) why not make a default value for boundary DONE
% 6) it may return some unexpected results when LOW==HIGH;  I DONT SEE IT,
%    but will add a comment.
% 7) Take a look at submission #9428 by John D'Errico how to implement
%    "boundary" effectively.  SORRY, TOO LAZY.
Narg = nargin;
error(nargchk(2,3,Narg,'struct'))
if Narg==2
    boundary = 'includeboth';
else
    boundary = varargin{1};
end
XoriginalDim = size(X);
X = X(:);
if numel(R) ~= 2,
    if ~isequal(size(R,2),2),
        error('RANGE input has too many columns.')
    end
    if ~isequal(size(R,1),numel(X))
        error('If RANGE is a matrix, X must be a vector of same length.')
    end
end
Rrep = R;
leftBound = R(:,1);
rightBound = R(:,2);
if any(leftBound > rightBound),
    error('Rows of RANGE must have form [LOW HIGH].')
end
switch boundary
    case 'includeboth',
        inRangeIX = (X >= leftBound) & (X <= rightBound);
    case 'includeleft',
        inRangeIX = (X >= leftBound) & (X < rightBound);
    case 'includeright',
        inRangeIX = (X > leftBound) & (X <= rightBound);
    case 'excludeboth',
        inRangeIX = (X > leftBound) & (X < rightBound);
    otherwise
        error('Valid options for third input are ''includeboth'', ''includeleft'', ''includeright'', ''excludeboth''.')
end
y = reshape(inRangeIX,XoriginalDim);
end