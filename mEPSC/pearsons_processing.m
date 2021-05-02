% Select the folder of the eventer output from a pearsons detection
folderpath = uigetdir();

% Load the summary.txt into a table
sumtable = readtable(fullfile(folderpath,"eventer.output/ALL_events/summary.txt"));
sumtable = table(sumtable.Var2,'RowNames',sumtable.Var1);

% Extract the rise and decay time constants of the pearsons fit
rise = table2array(sumtable({'Rise time constant of the model PSC fit (ms)'},:));
decay = table2array(sumtable({'Decay time constant of the model PSC fit (ms)'},:));
    
% change directory to the eventer.output folder
cd(fullfile(folderpath,"eventer.output"))

% generate directory listing
d=dir;
d=d(~ismember({d.name},{'.','..'}));
dirFlags = [d.isdir];
d = d(dirFlags);
waves = d(2:end);

% generate a list of file paths and extract the data
wave = zeros(size(waves,1),1);
noise = zeros(size(waves,1),1);
thresh = zeros(size(waves,1),1);
for i = 1:size(waves,1)
    wave_path(i,1) = string(fullfile(waves(i).folder,waves(i).name));
    file_path(i,1) = fullfile(wave_path(i),'summary.txt');
    [wave(i),noise(i),thresh(i),~] = pullnoise(file_path(i));
end

% outputs
[val,ind] = max(noise);
Wave_Number = wave(ind);
Noise_Level = val;
noise_levels = table(wave,noise);
[~,~,~,summary] = pullnoise(file_path(ind));
figure
plot(wave,noise)
title('Wave by noise')
ylabel('Standard deviation of the noise of the deconvoluted wave (a.u.)')
xlabel('Wave number')
box off
set(gca,'linewidth',2)
set(gcf,'color','w')
hold on
plot(wave,movmean(noise,7))
legend('Noise','Average','linewidth',1)


% create output
noise_output.filepath = file_path;
noise_output.max_noise_val = Noise_Level;
noise_output.max_noise_ind = Wave_Number;
noise_output.noise_levels = noise_levels;
noise_output.summary = summary;
noise_output.wave_path = wave_path;
noise_output.mean_noise = trimmean(noise_levels.noise,33,'floor');       

% tidy up
clear file path d waves wave noise thresh ans dirFlags i ind val Wave_Number 
clear file_path Noise_Level noise_levels summary wave_path

% cd back to same dir as raw recordings
cd ..

% Select recording summary came from
%title_str = "2. Select the recording summary.txt came from";
%if ~ispc; menu(title_str,'OK'); end
%[file,~,~] = uigetfile('*.phy');
%clear title_str ans

% save as name of file in dir
%a = split(file,'.');
%file = append(char(a(1)),'_noise.mat');
%save(file,'noise_output')

% save figure
%a = split(file,'.');
%file = append(char(a(1)),'.pdf');
%saveas(gcf,file)
%% Define Functions

function [wave_num, sd_noise, sd_thresh, summary] = pullnoise(filepath)
% simple function to pull the wave number, noise level and adjusted
% threshold from the summary.txt files generated by Eventer

% Output Arguments
%   wave_num: wave number to associate noise level with
%   sd_noise: sd of the noise (a.u.)
%   sd_thresh: scaled threshold sd

% Input Arguments:
%   filepath: character string of the fullfile where the summary.txt is
%             located

% ------------------------------------------------------------------- % 
    % Setup the Import Options
    opts = delimitedTextImportOptions("NumVariables", 2);

    % Specify range and delimiter
    opts.DataLines = [2, Inf];
    opts.Delimiter = ":";

    % Specify column names and types
    opts.VariableNames = ["Variable", "Data"];
    opts.VariableTypes = ["string", "double"];
    opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
    opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Import the data
    summary = readtable(filepath, opts);
    summary = table(summary.Data,'RowNames',summary.Variable);
    
    % Extract the variables of interest
    wave_num = table2array(summary({'Wave number'},:)); % wave number to associate noise level with
    sd_noise = table2array(summary({'Standard deviation of the noise of the deconvoluted wave (a.u.)'},:)); % scaled threshold sd
    sd_thresh = table2array(summary({'Scale factor of noise standard deviations for threshold setting'},:)); % scaled threshold sd
    %wave_num = summary.Data(2); % wave number to associate noise level with
    %sd_noise = summary.Data(11); % sd of the noise (a.u.)
    %sd_thresh = summary.Data(12); % scaled threshold sd
    
end

function [notes1] = notesimport(filename) 
% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 1);

% Specify range and delimiter
opts.DataLines = [1, 4];
opts.Delimiter = "";

% Specify column names and types
opts.VariableNames = "UsableForCompoundEventAnalysis";
opts.VariableTypes = "char";

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "UsableForCompoundEventAnalysis", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "UsableForCompoundEventAnalysis", "EmptyFieldRule", "auto");

% Import the data
notes1 = readtable(filename, opts);

% Convert to output type
notes1 = table2cell(notes1);
numIdx = cellfun(@(x) ~isnan(str2double(x)), notes1);
notes1(numIdx) = cellfun(@(x) {str2double(x)}, notes1(numIdx));

% Clear temporary variables
clear opts
end