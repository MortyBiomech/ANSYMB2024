clc
clear


%% Add necessary paths
main_project_path = 'D:\Morteza\MyProjects\ANSYMB2024\';

addpath(genpath([main_project_path, 'Code']));
addpath(genpath([main_project_path, 'data\7_STUDY\Epoched_data']));

data_path         = [main_project_path, 'data\'];
Code_path         = [main_project_path, 'Code\Matlab\data_processing\'];
all_STUDY_PATH    = [data_path, '7_STUDY\Epoched_data\', ...
                        'multiple_clustering\'];

icatimef_path     = [data_path, '5_single-subject-EEG-analysis\', ...
                        'timewarp_test\Epoched_data'];
epoched_data_path = [data_path, '6_Trials_Info_and_Epoched_data\'];
ersp_data_path    = [data_path, '7_STUDY\Epoched_data\Final_figures', ...
                        '\ERSP\Three Pressure Conditions\', ...
                        'p 0.01 ersp results\'];
Subject_ICs_in_clusters_path = [Code_path, ...
    'Group_Level_PostProcessing\Final_paper_plot_generation\', ...
    'Detailed_Analysis_on_TF_regions\', ...
    'extracting Subjects and ICs in the brain clusters'];

current_path = ['D:\Morteza\MyProjects\ANSYMB2024\Code\Matlab', ...
    '\data_processing\Group_Level_PostProcessing\', ...
    'Final_paper_plot_generation\', ...
    'Cross_correlation_Tracking_Error_Brain_Power'];

P1_color = [1, 115, 178]/255;
P3_color = [222, 143, 5]/255;
P6_color = [148, 73, 92]/255; %[148, 73, 92]/255;
colors = [P1_color; P3_color; P6_color];


%% Load necessary data structures
cd(current_path)
load("Subjects_Tracking_Error_warped_final.mat")
load("Subjects_Tracking_Error_final.mat")
load("Subjects_Epochs_toRemove_in_icatimef.mat")



%% Main Loop on Left-Parieto-Occipital Cluster Cross-Correlation

% load the subject-IC pairs for our STUDY
cd(Subject_ICs_in_clusters_path)
load("Subjects_ICs_in_clusters.mat")
cd(current_path)

idx_cluster = find(cellfun(@(x) strcmp(x, 'Left_Parieto_Occipital'), ...
    SUBJECTS_ICS(:, 1)));
Subjects = SUBJECTS_ICS{idx_cluster, 2}.Subjects + 4;
[Subjects_sorted, idx_subject_sort] = sort(Subjects, 2, "ascend");


Subject_list           = 5:18;
idx_member             = ismember(Subject_list, Subjects_sorted);

ICs                    = SUBJECTS_ICS{idx_cluster, 2}.ICs;
ICs_sorted_temp        = ICs(idx_subject_sort);
ICs_sorted             = nan(size(Subject_list));
ICs_sorted(idx_member) = ICs_sorted_temp;


data_structure = struct('EEG', [], 'Err', [], ...
    'Pressure', [], 'Score', [], 'Trial', []);
Main_data_structure = cell(length(Subject_list)+2, 3);
Main_data_structure(1:length(Subject_list), 3) = ...
    repmat({data_structure}, length(Subject_list), 1);
Main_data_structure{length(Subject_list)+1, 2} = 'EEG_time';
Main_data_structure{length(Subject_list)+2, 2} = 'EEG_freq';



cd(icatimef_path)
fileExt = '.icatimef';
for sub = 1:length(Subject_list)

    if idx_member(sub) == 1

        Main_data_structure{sub, 1} = ['Sub ', num2str(Subject_list(sub))];
        Main_data_structure{sub, 2} = ['IC ', num2str(ICs_sorted(sub))];
    
        % load icatime (just load the com_X, times, freqs)
        fileBaseName = ['S', num2str(Subject_list(sub))];
        chanList = ['comp', num2str(ICs_sorted(sub))];
        disp(['Loading S', num2str(Subject_list(sub)),'.icatimef ...']);
        icatimef = load('-mat', [ fileBaseName fileExt ], chanList, ...
            'times', 'freqs', 'trialinfo', 'parameters');
        trialinfo = icatimef.trialinfo;
        ic = icatimef.(['comp', num2str(ICs_sorted(sub))]);
        ic = ic.*conj(ic);
        % crop the ic data based on the timewarpms 
        idx = find(strcmp(icatimef.parameters, 'timewarpms'));
        timewarpms = icatimef.parameters{1, idx+1};
        times = icatimef.times;
        idx_to_keep = times < timewarpms(end);
        ic = ic(:, idx_to_keep, :);
        new_times = 100*(times(idx_to_keep)/timewarpms(end));
        
        % remove epochs which were marked in
        % main_cross_correlation_data_preparation.m 
        if ~isempty(Subject_Epochs_toRemove_in_icatimef{sub, 1})
            e = Subject_Epochs_toRemove_in_icatimef{sub, 1};
            ic(:, :, e) = [];
            trialinfo(e) = [];
        end
        trials = {trialinfo.trial}.';
        trials = cellfun(@(x) str2num(x), trials);

        % interpolate the Err to match the timing (percentage) of EEG
        Err = Subject_Tracking_Error_warped_final{sub, 1};
        Err_interp = interp1(1:size(Err, 2), Err', new_times, "linear");
        Err_interp = Err_interp';
        
        Main_data_structure{sub, 3}.EEG = ic;
        Main_data_structure{sub, 3}.Err = Err_interp;
        Main_data_structure{sub, 3}.Pressure = ...
            Subject_Tracking_Error_final{sub, 1}.pressure;
        Main_data_structure{sub, 3}.Score = ...
            Subject_Tracking_Error_final{sub, 1}.score;
        Main_data_structure{sub, 3}.Trial = trials;
        
        if sub == length(Subject_list)
            Main_data_structure{length(Subject_list)+1, 3} = ...
                icatimef.times;
            Main_data_structure{length(Subject_list)+2, 3} = ...
                icatimef.freqs;
        end
        clear icatimef
    else
        Main_data_structure{sub, 1} = ['Sub ', num2str(Subject_list(sub))];
        Main_data_structure{sub, 2} = 'No IC';
        Main_data_structure{sub, 3} = [];
    end

end



%% Moving from epochs space to trials space
cd(current_path)
Main_data_structure_trial = Main_data_structure;
for sub = 1:length(Subject_list)

    if isempty(Main_data_structure{sub, 3}), continue, end;
    
    EEG      = Main_data_structure{sub, 3}.EEG;
    Err      = Main_data_structure{sub, 3}.Err;
    Pressure = Main_data_structure{sub, 3}.Pressure;
    Score    = Main_data_structure{sub, 3}.Score;
    Trials   = Main_data_structure{sub, 3}.Trial;
    
    EEG_trial      = [];
    Err_trial      = [];
    pressure_trial = [];
    score_trial    = [];

    unique_trials = unique(Trials);
    for trial = 1:length(unique_trials)
        idx_trial      = find(Trials == unique_trials(trial));
        EEG_trial      = cat(3, EEG_trial, mean(EEG(:, :, idx_trial), 3));
        Err_trial      = cat(1, Err_trial, mean(Err(idx_trial, :), 1));
        pressure_trial = cat(1, pressure_trial, Pressure(idx_trial(1)));
        score_trial    = cat(1, score_trial, Score(idx_trial(1)));
    end

    idx_score_0 = find(score_trial == 0);
    if ~isempty(idx_score_0)
        EEG_trial(:, :, idx_score_0) = [];
        Err_trial(idx_score_0, :)    = [];
        pressure_trial(idx_score_0)  = [];
        score_trial(idx_score_0)     = [];
        unique_trials(idx_score_0)   = [];
    end

    Main_data_structure_trial{sub, 3}.EEG = EEG_trial;
    Main_data_structure_trial{sub, 3}.Err = Err_trial;
    Main_data_structure_trial{sub, 3}.Pressure = pressure_trial;
    Main_data_structure_trial{sub, 3}.Score = score_trial;
    Main_data_structure_trial{sub, 3}.Trial = unique_trials;

end

clear Main_data_structure;


%% Lets do the cross-correlation 

% Option 1: 
%   Using raw EEG power (conjucated time-frequency decompostions)
%   with no baseline correction.
% 
% Option 1A: 
%   Cross-Correlation on Residuals

% ---------------------------------- %
% Summary of discussions with ChatGPT:
% ---------------------------------- %
% Goal:
% You don't just want to see that error and θ/α have a similar average 
% shape (which is largely driven by task phase).
% You want to know:
% When error is higher than usual on a given trial at a given phase, 
% does θ/α power also change? And does it change before or after the 
% error change?
% To do that, you remove the phase-locked average and analyse 
% cycle-to-cycle fluctuations (residuals).

% Option 1A: 
% Use raw EEG power (conjucated time-frequency decompostions)
% with no baseline correction.

freqs = Main_data_structure_trial{end, 3};
idx_theta = (freqs > 4) & (freqs < 8);
idx_alpha = (freqs > 8) & (freqs < 14);

nTime = length(new_times);
timeIdx = 1:nTime;              % 1..nTime = 0..100% cycle
maxLagPct    = 50;              % max lag in % of cycle (e.g. ±10%)
maxLagSamples = round(maxLagPct/100 * nTime);
lagVec       = -maxLagSamples:maxLagSamples;
nLags        = numel(lagVec);



nCond = 3;
Conditions = [1, 3, 6];
R_theta = cell(length(Subject_list), nCond);
R_alpha = cell(length(Subject_list), nCond);
for sub = 1:length(Subject_list)
    
    if isempty(Main_data_structure_trial{sub, 3}), continue, end;

    % Tracking Error
    Err = Main_data_structure_trial{sub, 3}.Err;
    % Theta and Alpha Power
    EEG_theta = mean(Main_data_structure_trial{sub, 3}.EEG(idx_theta, :, :), 1);
    EEG_theta = squeeze(EEG_theta)';
    EEG_alpha = mean(Main_data_structure_trial{sub, 3}.EEG(idx_alpha, :, :), 1);
    EEG_alpha = squeeze(EEG_alpha)';

    pressure = Main_data_structure_trial{sub, 3}.Pressure;
    for c = 1:nCond
        
        pressure_idx = pressure == Conditions(c);
        curErr = Err(pressure_idx, :);               % [nTrials_sc x nTime]
        curPow_theta = EEG_theta(pressure_idx, :);   % [nTrials_sc x nTime]
        curPow_alpha = EEG_alpha(pressure_idx, :);   % [nTrials_sc x nTime]
        

        % --- Step 1: mean waveform per subject & condition ---
        meanErr = mean(curErr, 1);                            % [1 x nTime]
        meanPow_theta = mean(curPow_theta, 1);                % [1 x nTime]
        meanPow_alpha = mean(curPow_alpha, 1);                % [1 x nTime]
        
        % --- Step 2: residuals = epoch - mean waveform ---
        errRes = curErr - meanErr;                      % [nTrials x nTime]
        powRes_theta = curPow_theta - meanPow_theta;    % [nTrials x nTime]
        powRes_alpha = curPow_alpha - meanPow_alpha;    % [nTrials x nTime]
        
        nTrials_sc = size(errRes, 1);
        nTime      = size(errRes, 2);
        
       
        rLags_trials_theta = nan(nTrials_sc, nLags);
        for k = 1:nTrials_sc
            e = errRes(k, :);       % [1 x nTime]
            p_theta = powRes_theta(k, :);
            p_alpha = powRes_alpha(k, :);

            % e = e - mean(e, 'omitnan');
            % p = p - mean(p, 'omitnan');

            rLags_trials_theta(k, :) = xcorr(e, p_theta, maxLagSamples, "normalized");
            rLags_trials_alpha(k, :) = xcorr(e, p_alpha, maxLagSamples, "normalized");
        
        end



        % e = mean(errRes, 1);
        % p = mean(powRes_theta, 1);
        % rLags_trials = xcorr(e, p, maxLagSamples, "normalized");

        R_theta{sub, c} = rLags_trials_theta;
        R_alpha{sub, c} = rLags_trials_alpha;
    end
end



%% plot the Cross-Correlation results

% Get monitor information
monitors = get(0, 'MonitorPositions');

fig = figure('name', ['Cross-Correlation of EEG power vs Tracking Error Residuals'], ...
    'InvertHardcopy', 'off', 'PaperType', 'a2', ...
    'PaperOrientation', 'landscape', ...
    'Resize', 'off');

% For second monitor (row 2), add drawnow before setting position
drawnow;  % Let MATLAB finish drawing on primary monitor first
pause(0.1);  % Short pause helps

% Now set position on second monitor
% set(gcf, 'Position', monitors(1, :));  % Use entire second monitor
set(fig, 'Position', [monitors(1,1)-200, monitors(1,2)+600, 1200, 800]);  % Use entire second monitor


% 13 built-in marker types
markerList = {'o','s','d','^','v','>','<','p','h','s','^','v','p','o'};
% colors_markers = lines(length(Subject_list));   % 14 distinct colors
colors_markers = [
    0.00 0.45 0.74  % blue
    0.85 0.33 0.10  % reddish orange
    0.93 0.69 0.13  % yellow-ish
    0.49 0.18 0.56  % purple
    0.47 0.67 0.19  % green
    0.30 0.75 0.93  % cyan
    0.64 0.08 0.18  % dark red
    0.13 0.35 0.60  % dark blue
    0.76 0.49 0.27  % brown-ish
    0.55 0.57 0.00  % olive
    0.65 0.70 0.96  % light bluish
    0.49 0.39 0.64  % desaturated purple
    0.70 0.67 0.18  % mustard
    0.40 0.40 0.40  % gray
];


% generate the %cycle corresponding Lag
maxLagPct1    = 10;              
maxLagSamples1 = round(maxLagPct1/100 * nTime);
maxLagPct2    = 20;              
maxLagSamples2 = round(maxLagPct2/100 * nTime);
maxLagPct3    = 30;              
maxLagSamples3 = round(maxLagPct3/100 * nTime);
maxLagPct4    = 40;              
maxLagSamples4 = round(maxLagPct4/100 * nTime);
maxLagPct5    = 50;              
maxLagSamples5 = round(maxLagPct5/100 * nTime);

LagSamplesnegativ = -1*[maxLagSamples5, maxLagSamples4, maxLagSamples3, ...
    maxLagSamples2,maxLagSamples1];
LagSamplespositive = -1*fliplr(LagSamplesnegativ);
xticks = [LagSamplesnegativ 0 LagSamplespositive];
xticklabels = {'-50', '-40', '-30', '-20', '-10', '0', ...
    '+10', '+20', '+30', '+40', '+50'};

titles = {'Low', 'Medium', 'High'};

legends = cellfun(@(x) strcat('S', num2str(x)), ...
    num2cell(Subjects_sorted), 'UniformOutput', false);
    


% theta band
R_theta_abs = cellfun(@(x) abs(x), R_theta, 'UniformOutput', false);
R_theta_mean_abs = R_theta_abs;
R_theta_mean_abs(10, :) = [];
R_theta_mean_abs = cellfun(@(x) mean(x, 1), R_theta_mean_abs, 'UniformOutput', false);
R_theta_mean_abs_P1 = cell2mat(R_theta_mean_abs(:, 1));
R_theta_mean_abs_P3 = cell2mat(R_theta_mean_abs(:, 2));
R_theta_mean_abs_P6 = cell2mat(R_theta_mean_abs(:, 3));

R_theta_mean = R_theta;
R_theta_mean(10, :) = [];
R_theta_mean = cellfun(@(x) mean(x, 1), R_theta_mean, 'UniformOutput', false);
R_theta_mean_P1 = cell2mat(R_theta_mean(:, 1));
R_theta_mean_P3 = cell2mat(R_theta_mean(:, 2));
R_theta_mean_P6 = cell2mat(R_theta_mean(:, 3));

% alpha band
R_alpha_abs = cellfun(@(x) abs(x), R_alpha, 'UniformOutput', false);
R_alpha_mean_abs = R_alpha_abs;
R_alpha_mean_abs(10, :) = [];
R_alpha_mean_abs = cellfun(@(x) mean(x, 1), R_alpha_mean_abs, 'UniformOutput', false);
R_alpha_mean_abs_P1 = cell2mat(R_alpha_mean_abs(:, 1));
R_alpha_mean_abs_P3 = cell2mat(R_alpha_mean_abs(:, 2));
R_alpha_mean_abs_P6 = cell2mat(R_alpha_mean_abs(:, 3));

R_alpha_mean = R_alpha;
R_alpha_mean(10, :) = [];
R_alpha_mean = cellfun(@(x) mean(x, 1), R_alpha_mean, 'UniformOutput', false);
R_alpha_mean_P1 = cell2mat(R_alpha_mean(:, 1));
R_alpha_mean_P3 = cell2mat(R_alpha_mean(:, 2));
R_alpha_mean_P6 = cell2mat(R_alpha_mean(:, 3));

R = cat(3, R_alpha, R_theta);
R_abs = cat(3, R_alpha_abs, R_theta_abs);

% three rows for three conditions, 
% two coloums for alpha and theta
all_powers_abs = cell(3, 2); 
all_powers_abs{1, 1} = R_alpha_mean_abs_P1;
all_powers_abs{2, 1} = R_alpha_mean_abs_P3;
all_powers_abs{3, 1} = R_alpha_mean_abs_P6;
all_powers_abs{1, 2} = R_theta_mean_abs_P1;
all_powers_abs{2, 2} = R_theta_mean_abs_P3;
all_powers_abs{3, 2} = R_theta_mean_abs_P6;

all_powers = cell(3, 2); 
all_powers{1, 1} = R_alpha_mean_P1;
all_powers{2, 1} = R_alpha_mean_P3;
all_powers{3, 1} = R_alpha_mean_P6;
all_powers{1, 2} = R_theta_mean_P1;
all_powers{2, 2} = R_theta_mean_P3;
all_powers{3, 2} = R_theta_mean_P6;


X = [lagVec, fliplr(lagVec)];
ylabels = {'magnitude |r(Lag)|', 'signed r(Lag)'};
Subject_sigend_R = cell(length(Subject_list), 3, 2);
RR = {R_abs, R};
for power = 1:2 % theta and alpha

    for i = 1:3

        subplot(2, 3, (power-1)*3 + i); hold on;

        % absolute r(Lag) or signed r(Lag)?
        % signal = all_powers_abs{i, power};
        signal = all_powers{i, power};
        ylabel_s = ylabels{2};
        RRR = RR{2};



        % Y = [mean(signal, 1) + ...
        %        std(signal, [], 1)./sqrt(size(signal,1 )), ...
        %      fliplr(mean(signal, 1) - ...
        %      std(signal, [], 1)./sqrt(size(signal,1 )))];
        Y = [mean(signal, 1) - std(signal, [], 1), ...
             fliplr(mean(signal, 1) + std(signal, [], 1))];
        
        fill(X, Y, colors(i, :), 'EdgeColor', 'none', 'FaceColor', colors(i, :), ...
            'FaceAlpha', 0.2, 'HandleVisibility', 'off')
        plot(lagVec, mean(signal, 1), 'Color', colors(i, :), ...
            'LineWidth', 2, 'HandleVisibility', 'off')
        xline(0, 'LineStyle', '--', 'HandleVisibility', 'off')
        xlim([lagVec(1) lagVec(end)])
        
        
        
        
        for sub = 1:length(Subject_list)
    
            if isempty(Main_data_structure_trial{sub, 3}), continue, end;
    
            [maxR maxLagSamples] = max(RRR{sub,i, power}, [], 2);
            sigend_R = cellfun(@(x1, x2) x1(1, x2), ...
                num2cell(R{sub, i, power}, 2), num2cell(maxLagSamples));
            Subject_sigend_R{sub, i, power} = sigend_R;
            
            meanLag = mean(lagVec(maxLagSamples));
            stdLag = std(lagVec(maxLagSamples), []);
            semLag = std(lagVec(maxLagSamples), [])/sqrt(numel(maxLagSamples));
            meanR = mean(maxR);
            stdR = std(maxR, []);
            semR = std(maxR, [])/sqrt(numel(maxLagSamples));
            
           
            % errorbar(meanLag, meanR, semR, semR, semLag, semLag, ...
            %     "Color", 0.7*[1 1 1], 'HandleVisibility', 'off')
            errorbar(meanLag, meanR, stdR, stdR, stdLag, stdLag, ...
                "Color", 0.7*[1 1 1], 'HandleVisibility', 'off')
            scatter(meanLag, meanR, 40, 'Marker', markerList{sub}, ...
                'LineWidth', 2, 'MarkerFaceColor', 'w', ...
                'MarkerEdgeColor', colors_markers(sub, :));
        end
        ylim([0 1])

        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, ...
            'XTickLabelRotation', 45);
        set(gca, 'FontSize', 16)

        if power == 2
            xlh = xlabel('Lag (% Cycle)', 'FontWeight', 'bold', ...
                'FontSize', 18);
            xlh.Position(2) = -0.15;
        else
            tth = title([titles{i}, ' Pressure'], ...
                'FontSize', 18, 'FontWeight', 'bold');
            tth.Position(2) = 1.1;
        end

        if i == 1
            if power == 1
                ylh = ylabel(sprintf(['\\bf\\alpha Power vs. Tracking Error\n', ...
                              '\\bfCorrelation ' ylabel_s]), ...
                     'Interpreter','tex', ...      % <-- important
                     'FontSize',18);
                ylh.Position(1) = lagVec(1) - 30;
            else
                ylh = ylabel(sprintf(['\\bf\\theta Power vs. Tracking Error\n', ...
                              '\\bfCorrelation ' ylabel_s]), ...
                     'Interpreter','tex', ...      % <-- important
                     'FontSize',18);
                ylh.Position(1) = lagVec(1) - 30;
            end
        end

        if power == 2 && i == 3
            legend(legends, 'FontName', 'Arial', 'FontSize', 12, ...
                'Box', 'off')
        end

    end
    
    

end




%% Option 1B: 
%  Cross-Correlation on actual signal (no subtraction of mean signal)
%  EEG is raw power (no baseline correction)


freqs = Main_data_structure_trial{end, 3};
idx_theta = (freqs > 4) & (freqs < 8);
idx_alpha = (freqs > 8) & (freqs < 14);

nTime = length(new_times);
timeIdx = 1:nTime;              % 1..nTime = 0..100% cycle
maxLagPct    = 50;              % max lag in % of cycle (e.g. ±10%)
maxLagSamples = round(maxLagPct/100 * nTime);
lagVec       = -maxLagSamples:maxLagSamples;
nLags        = numel(lagVec);



nCond = 3;
Conditions = [1, 3, 6];
R_theta = cell(length(Subject_list), nCond);
R_alpha = cell(length(Subject_list), nCond);
for sub = 1:length(Subject_list)
    
    if isempty(Main_data_structure_trial{sub, 3}), continue, end;

    % Tracking Error
    Err = Main_data_structure_trial{sub, 3}.Err;
    % Theta and Alpha Power
    EEG_theta = mean(Main_data_structure_trial{sub, 3}.EEG(idx_theta, :, :), 1);
    EEG_theta = squeeze(EEG_theta)';
    EEG_alpha = mean(Main_data_structure_trial{sub, 3}.EEG(idx_alpha, :, :), 1);
    EEG_alpha = squeeze(EEG_alpha)';

    pressure = Main_data_structure_trial{sub, 3}.Pressure;
    for c = 1:nCond
        
        pressure_idx = pressure == Conditions(c);
        curErr = Err(pressure_idx, :);               % [nTrials_sc x nTime]
        curPow_theta = EEG_theta(pressure_idx, :);   % [nTrials_sc x nTime]
        curPow_alpha = EEG_alpha(pressure_idx, :);   % [nTrials_sc x nTime]
        

        nTrials_sc = size(curErr, 1);
        nTime      = size(curErr, 2);
        
       
        rLags_trials_theta = nan(nTrials_sc, nLags);
        rLags_trials_alpha = nan(nTrials_sc, nLags);
        for k = 1:nTrials_sc
            e = curErr(k, :);       % [1 x nTime]
            p_theta = curPow_theta(k, :);
            p_alpha = curPow_alpha(k, :);

            % e = e - mean(e, 'omitnan');
            % p = p - mean(p, 'omitnan');

            rLags_trials_theta(k, :) = xcorr(e, p_theta, maxLagSamples, "normalized");
            rLags_trials_alpha(k, :) = xcorr(e, p_alpha, maxLagSamples, "normalized");
        
        end



        % e = mean(errRes, 1);
        % p = mean(powRes_theta, 1);
        % rLags_trials = xcorr(e, p, maxLagSamples, "normalized");

        R_theta{sub, c} = rLags_trials_theta;
        R_alpha{sub, c} = rLags_trials_alpha;
    end
end




%% plot the Cross-Correlation results

% Get monitor information
monitors = get(0, 'MonitorPositions');

fig = figure('name', ['Cross-Correlation of EEG power vs Tracking Error actual signals'], ...
    'InvertHardcopy', 'off', 'PaperType', 'a2', ...
    'PaperOrientation', 'landscape', ...
    'Resize', 'off');

% For second monitor (row 2), add drawnow before setting position
drawnow;  % Let MATLAB finish drawing on primary monitor first
pause(0.1);  % Short pause helps

% Now set position on second monitor
% set(gcf, 'Position', monitors(1, :));  % Use entire second monitor
set(fig, 'Position', [monitors(1,1)-200, monitors(1,2)+600, 1200, 800]);  % Use entire second monitor


% 13 built-in marker types
markerList = {'o','s','d','^','v','>','<','p','h','s','^','v','p','o'};
% colors_markers = lines(length(Subject_list));   % 14 distinct colors
colors_markers = [
    0.00 0.45 0.74  % blue
    0.85 0.33 0.10  % reddish orange
    0.93 0.69 0.13  % yellow-ish
    0.49 0.18 0.56  % purple
    0.47 0.67 0.19  % green
    0.30 0.75 0.93  % cyan
    0.64 0.08 0.18  % dark red
    0.13 0.35 0.60  % dark blue
    0.76 0.49 0.27  % brown-ish
    0.55 0.57 0.00  % olive
    0.65 0.70 0.96  % light bluish
    0.49 0.39 0.64  % desaturated purple
    0.70 0.67 0.18  % mustard
    0.40 0.40 0.40  % gray
];






% theta band
R_theta_abs = cellfun(@(x) abs(x), R_theta, 'UniformOutput', false);
R_theta_mean_abs = R_theta_abs;
R_theta_mean_abs(10, :) = [];
R_theta_mean_abs = cellfun(@(x) mean(x, 1), R_theta_mean_abs, 'UniformOutput', false);
R_theta_mean_abs_P1 = cell2mat(R_theta_mean_abs(:, 1));
R_theta_mean_abs_P3 = cell2mat(R_theta_mean_abs(:, 2));
R_theta_mean_abs_P6 = cell2mat(R_theta_mean_abs(:, 3));

R_theta_mean = R_theta;
R_theta_mean(10, :) = [];
R_theta_mean = cellfun(@(x) mean(x, 1), R_theta_mean, 'UniformOutput', false);
R_theta_mean_P1 = cell2mat(R_theta_mean(:, 1));
R_theta_mean_P3 = cell2mat(R_theta_mean(:, 2));
R_theta_mean_P6 = cell2mat(R_theta_mean(:, 3));

% alpha band
R_alpha_abs = cellfun(@(x) abs(x), R_alpha, 'UniformOutput', false);
R_alpha_mean_abs = R_alpha_abs;
R_alpha_mean_abs(10, :) = [];
R_alpha_mean_abs = cellfun(@(x) mean(x, 1), R_alpha_mean_abs, 'UniformOutput', false);
R_alpha_mean_abs_P1 = cell2mat(R_alpha_mean_abs(:, 1));
R_alpha_mean_abs_P3 = cell2mat(R_alpha_mean_abs(:, 2));
R_alpha_mean_abs_P6 = cell2mat(R_alpha_mean_abs(:, 3));

R_alpha_mean = R_alpha;
R_alpha_mean(10, :) = [];
R_alpha_mean = cellfun(@(x) mean(x, 1), R_alpha_mean, 'UniformOutput', false);
R_alpha_mean_P1 = cell2mat(R_alpha_mean(:, 1));
R_alpha_mean_P3 = cell2mat(R_alpha_mean(:, 2));
R_alpha_mean_P6 = cell2mat(R_alpha_mean(:, 3));

R = cat(3, R_alpha, R_theta);
R_abs = cat(3, R_alpha_abs, R_theta_abs);

% three rows for three conditions, 
% two coloums for alpha and theta
all_powers_abs = cell(3, 2); 
all_powers_abs{1, 1} = R_alpha_mean_abs_P1;
all_powers_abs{2, 1} = R_alpha_mean_abs_P3;
all_powers_abs{3, 1} = R_alpha_mean_abs_P6;
all_powers_abs{1, 2} = R_theta_mean_abs_P1;
all_powers_abs{2, 2} = R_theta_mean_abs_P3;
all_powers_abs{3, 2} = R_theta_mean_abs_P6;

all_powers = cell(3, 2); 
all_powers{1, 1} = R_alpha_mean_P1;
all_powers{2, 1} = R_alpha_mean_P3;
all_powers{3, 1} = R_alpha_mean_P6;
all_powers{1, 2} = R_theta_mean_P1;
all_powers{2, 2} = R_theta_mean_P3;
all_powers{3, 2} = R_theta_mean_P6;


X = [lagVec, fliplr(lagVec)];
ylabels = {'magnitude |r(Lag)|', 'signed r(Lag)'};
Subject_sigend_R = cell(length(Subject_list), 3, 2);
RR = {R_abs, R};
for power = 1:2 % theta and alpha

    for i = 1:3

        subplot(2, 3, (power-1)*3 + i); hold on;

        % absolute r(Lag) or signed r(Lag)?
        % signal = all_powers_abs{i, power};
        signal = all_powers{i, power};
        ylabel_s = ylabels{2};
        RRR = RR{2};


        % Y = [mean(signal, 1) + ...
        %        std(signal, [], 1)./sqrt(size(signal,1 )), ...
        %      fliplr(mean(signal, 1) - ...
        %      std(signal, [], 1)./sqrt(size(signal,1 )))];
        Y = [mean(signal, 1) - std(signal, [], 1), ...
             fliplr(mean(signal, 1) + std(signal, [], 1))];
        
        fill(X, Y, colors(i, :), 'EdgeColor', 'none', 'FaceColor', colors(i, :), ...
            'FaceAlpha', 0.2, 'HandleVisibility', 'off')
        plot(lagVec, mean(signal, 1), 'Color', colors(i, :), ...
            'LineWidth', 2, 'HandleVisibility', 'off')
        xline(0, 'LineStyle', '--', 'HandleVisibility', 'off')
        xlim([lagVec(1) lagVec(end)])
        
        
        
        
        for sub = 1:length(Subject_list)
    
            if isempty(Main_data_structure_trial{sub, 3}), continue, end;
    
            [maxR maxLagSamples] = max(RRR{sub,i, power}, [], 2);
            sigend_R = cellfun(@(x1, x2) x1(1, x2), ...
                num2cell(R{sub, i, power}, 2), num2cell(maxLagSamples));
            Subject_sigend_R{sub, i, power} = sigend_R;
            
            meanLag = mean(lagVec(maxLagSamples));
            stdLag = std(lagVec(maxLagSamples), []);
            semLag = std(lagVec(maxLagSamples), [])/sqrt(numel(maxLagSamples));
            meanR = mean(maxR);
            stdR = std(maxR, []);
            semR = std(maxR, [])/sqrt(numel(maxLagSamples));
            
           
            % errorbar(meanLag, meanR, semR, semR, semLag, semLag, ...
            %     "Color", 0.7*[1 1 1], 'HandleVisibility', 'off')
            errorbar(meanLag, meanR, stdR, stdR, stdLag, stdLag, ...
                "Color", 0.7*[1 1 1], 'HandleVisibility', 'off')
            scatter(meanLag, meanR, 40, 'Marker', markerList{sub}, ...
                'LineWidth', 2, 'MarkerFaceColor', 'w', ...
                'MarkerEdgeColor', colors_markers(sub, :));
        end
        ylim([0 1])

        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, ...
            'XTickLabelRotation', 45);
        set(gca, 'FontSize', 16)

        if power == 2
            xlh = xlabel('Lag (% Cycle)', 'FontWeight', 'bold', ...
                'FontSize', 18);
            xlh.Position(2) = -0.15;
        else
            tth = title([titles{i}, ' Pressure'], ...
                'FontSize', 18, 'FontWeight', 'bold');
            tth.Position(2) = 1.1;
        end

        if i == 1
            if power == 1
                ylh = ylabel(sprintf(['\\bf\\alpha Power vs. Tracking Error\n', ...
                              '\\bfCorrelation ' ylabel_s]), ...
                     'Interpreter','tex', ...      % <-- important
                     'FontSize',18);
                ylh.Position(1) = lagVec(1) - 30;
            else
                ylh = ylabel(sprintf(['\\bf\\theta Power vs. Tracking Error\n', ...
                              '\\bfCorrelation ' ylabel_s]), ...
                     'Interpreter','tex', ...      % <-- important
                     'FontSize',18);
                ylh.Position(1) = lagVec(1) - 30;
            end
        end

        if power == 2 && i == 3
            legend(legends, 'FontName', 'Arial', 'FontSize', 12, ...
                'Box', 'off')
        end

    end
    
    

end





%% Option 2: 
%   Using baseline-corrected EEG power (dB); ERSP (dB)
% 
% Option 2A: 
%   Cross-Correlation on Residuals

% ---------------------------------- %
% Summary of discussions with ChatGPT:
% ---------------------------------- %
% Goal:
% You don't just want to see that error and θ/α have a similar average 
% shape (which is largely driven by task phase).
% You want to know:
% When error is higher than usual on a given trial at a given phase, 
% does θ/α power also change? And does it change before or after the 
% error change?
% To do that, you remove the phase-locked average and analyse 
% cycle-to-cycle fluctuations (residuals).

% Option 2A: 
% Use baseline-corrected EEG power  


% Calculate the baseline per subject
Subject_baseline = cell(length(Subject_list), 1);
Conditions = [1, 3, 6];
for sub = 1:length(Subject_list)
    
    if isempty(Main_data_structure_trial{sub, 3}), continue, end;
    
    pressure = Main_data_structure_trial{sub, 3}.Pressure;
    EEG = Main_data_structure_trial{sub, 3}.EEG;
    mean_EEG = zeros(size(EEG, 1), 1);
    for i = 1:3
        idx_pressure = find(pressure == Conditions(i));
        EEG_p = EEG(:, :, idx_pressure);
        mean_EEG_p = mean(EEG_p, 3);
        mean_EEG = mean_EEG + mean(mean_EEG_p, 2)/3;
    end

    Subject_baseline{sub, 1} = mean_EEG;

end


% Correct Baseline in EEG power (no log-transformed)
Main_data_structure_trial_baseln = Main_data_structure_trial;
for sub = 1:length(Subject_list)
    
    if isempty(Main_data_structure_trial{sub, 3}), continue, end;
    
    EEG = Main_data_structure_trial{sub, 3}.EEG;
    baseline = Subject_baseline{sub, 1};
    
    EEG = EEG ./ repmat(baseline, 1, size(EEG, 2), size(EEG, 3));

    Main_data_structure_trial_baseln{sub, 3} = EEG;

end


freqs = Main_data_structure_trial{end, 3};
idx_theta = (freqs > 4) & (freqs < 8);
idx_alpha = (freqs > 8) & (freqs < 14);

nTime = length(new_times);
timeIdx = 1:nTime;              % 1..nTime = 0..100% cycle
maxLagPct    = 50;              % max lag in % of cycle (e.g. ±10%)
maxLagSamples = round(maxLagPct/100 * nTime);
lagVec       = -maxLagSamples:maxLagSamples;
nLags        = numel(lagVec);



nCond = 3;
R_theta = cell(length(Subject_list), nCond);
R_alpha = cell(length(Subject_list), nCond);
for sub = 1:length(Subject_list)
    
    if isempty(Main_data_structure_trial{sub, 3}), continue, end;

    % Tracking Error
    Err = Main_data_structure_trial{sub, 3}.Err;
    % Theta and Alpha Power
    EEG_theta = mean(Main_data_structure_trial_baseln{sub, 3}(idx_theta, :, :), 1);
    EEG_theta = squeeze(EEG_theta)';
    EEG_alpha = mean(Main_data_structure_trial_baseln{sub, 3}(idx_alpha, :, :), 1);
    EEG_alpha = squeeze(EEG_alpha)';

    pressure = Main_data_structure_trial{sub, 3}.Pressure;
    for c = 1:nCond
        
        pressure_idx = pressure == Conditions(c);
        curErr = Err(pressure_idx, :);               % [nTrials_sc x nTime]
        curPow_theta = EEG_theta(pressure_idx, :);   % [nTrials_sc x nTime]
        curPow_alpha = EEG_alpha(pressure_idx, :);   % [nTrials_sc x nTime]
        

        % --- Step 1: mean waveform per subject & condition ---
        meanErr = mean(curErr, 1);                            % [1 x nTime]
        meanPow_theta = mean(curPow_theta, 1);                % [1 x nTime]
        meanPow_alpha = mean(curPow_alpha, 1);                % [1 x nTime]
        
        % --- Step 2: residuals = epoch - mean waveform ---
        errRes = curErr - meanErr;                      % [nTrials x nTime]
        powRes_theta = curPow_theta - meanPow_theta;    % [nTrials x nTime]
        powRes_alpha = curPow_alpha - meanPow_alpha;    % [nTrials x nTime]
        
        nTrials_sc = size(errRes, 1);
        nTime      = size(errRes, 2);
        
       
        rLags_trials_theta = nan(nTrials_sc, nLags);
        for k = 1:nTrials_sc
            e = errRes(k, :);       % [1 x nTime]
            p_theta = powRes_theta(k, :);
            p_alpha = powRes_alpha(k, :);

            % e = e - mean(e, 'omitnan');
            % p = p - mean(p, 'omitnan');

            rLags_trials_theta(k, :) = xcorr(e, p_theta, maxLagSamples, "normalized");
            rLags_trials_alpha(k, :) = xcorr(e, p_alpha, maxLagSamples, "normalized");
        
        end



        % e = mean(errRes, 1);
        % p = mean(powRes_theta, 1);
        % rLags_trials = xcorr(e, p, maxLagSamples, "normalized");

        R_theta{sub, c} = rLags_trials_theta;
        R_alpha{sub, c} = rLags_trials_alpha;
    end
end


%% plot the Cross-Correlation results

% Get monitor information
monitors = get(0, 'MonitorPositions');

fig = figure('name', ...
    ['Cross-Correlation of baseline-corrected EEG power', ...
    ' vs Tracking Error Residuals'], ...
    'InvertHardcopy', 'off', 'PaperType', 'a2', ...
    'PaperOrientation', 'landscape', ...
    'Resize', 'off');

% For second monitor (row 2), add drawnow before setting position
drawnow;  % Let MATLAB finish drawing on primary monitor first
pause(0.1);  % Short pause helps

% Now set position on second monitor
% set(gcf, 'Position', monitors(1, :));  % Use entire second monitor
set(fig, 'Position', [monitors(1,1)-200, monitors(1,2)+600, 1200, 800]);  % Use entire second monitor


% 13 built-in marker types
markerList = {'o','s','d','^','v','>','<','p','h','s','^','v','p','o'};
% colors_markers = lines(length(Subject_list));   % 14 distinct colors
colors_markers = [
    0.00 0.45 0.74  % blue
    0.85 0.33 0.10  % reddish orange
    0.93 0.69 0.13  % yellow-ish
    0.49 0.18 0.56  % purple
    0.47 0.67 0.19  % green
    0.30 0.75 0.93  % cyan
    0.64 0.08 0.18  % dark red
    0.13 0.35 0.60  % dark blue
    0.76 0.49 0.27  % brown-ish
    0.55 0.57 0.00  % olive
    0.65 0.70 0.96  % light bluish
    0.49 0.39 0.64  % desaturated purple
    0.70 0.67 0.18  % mustard
    0.40 0.40 0.40  % gray
];


% generate the %cycle corresponding Lag
maxLagPct1    = 10;              
maxLagSamples1 = round(maxLagPct1/100 * nTime);
maxLagPct2    = 20;              
maxLagSamples2 = round(maxLagPct2/100 * nTime);
maxLagPct3    = 30;              
maxLagSamples3 = round(maxLagPct3/100 * nTime);
maxLagPct4    = 40;              
maxLagSamples4 = round(maxLagPct4/100 * nTime);
maxLagPct5    = 50;              
maxLagSamples5 = round(maxLagPct5/100 * nTime);

LagSamplesnegativ = -1*[maxLagSamples5, maxLagSamples4, maxLagSamples3, ...
    maxLagSamples2,maxLagSamples1];
LagSamplespositive = -1*fliplr(LagSamplesnegativ);
xticks = [LagSamplesnegativ 0 LagSamplespositive];
xticklabels = {'-50', '-40', '-30', '-20', '-10', '0', ...
    '+10', '+20', '+30', '+40', '+50'};

titles = {'Low', 'Medium', 'High'};

legends = cellfun(@(x) strcat('S', num2str(x)), ...
    num2cell(Subjects_sorted), 'UniformOutput', false);
    


% theta band
R_theta_abs = cellfun(@(x) abs(x), R_theta, 'UniformOutput', false);
R_theta_mean_abs = R_theta_abs;
R_theta_mean_abs(10, :) = [];
R_theta_mean_abs = cellfun(@(x) mean(x, 1), R_theta_mean_abs, 'UniformOutput', false);
R_theta_mean_abs_P1 = cell2mat(R_theta_mean_abs(:, 1));
R_theta_mean_abs_P3 = cell2mat(R_theta_mean_abs(:, 2));
R_theta_mean_abs_P6 = cell2mat(R_theta_mean_abs(:, 3));

R_theta_mean = R_theta;
R_theta_mean(10, :) = [];
R_theta_mean = cellfun(@(x) mean(x, 1), R_theta_mean, 'UniformOutput', false);
R_theta_mean_P1 = cell2mat(R_theta_mean(:, 1));
R_theta_mean_P3 = cell2mat(R_theta_mean(:, 2));
R_theta_mean_P6 = cell2mat(R_theta_mean(:, 3));

% alpha band
R_alpha_abs = cellfun(@(x) abs(x), R_alpha, 'UniformOutput', false);
R_alpha_mean_abs = R_alpha_abs;
R_alpha_mean_abs(10, :) = [];
R_alpha_mean_abs = cellfun(@(x) mean(x, 1), R_alpha_mean_abs, 'UniformOutput', false);
R_alpha_mean_abs_P1 = cell2mat(R_alpha_mean_abs(:, 1));
R_alpha_mean_abs_P3 = cell2mat(R_alpha_mean_abs(:, 2));
R_alpha_mean_abs_P6 = cell2mat(R_alpha_mean_abs(:, 3));

R_alpha_mean = R_alpha;
R_alpha_mean(10, :) = [];
R_alpha_mean = cellfun(@(x) mean(x, 1), R_alpha_mean, 'UniformOutput', false);
R_alpha_mean_P1 = cell2mat(R_alpha_mean(:, 1));
R_alpha_mean_P3 = cell2mat(R_alpha_mean(:, 2));
R_alpha_mean_P6 = cell2mat(R_alpha_mean(:, 3));

R = cat(3, R_alpha, R_theta);
R_abs = cat(3, R_alpha_abs, R_theta_abs);

% three rows for three conditions, 
% two coloums for alpha and theta
all_powers_abs = cell(3, 2); 
all_powers_abs{1, 1} = R_alpha_mean_abs_P1;
all_powers_abs{2, 1} = R_alpha_mean_abs_P3;
all_powers_abs{3, 1} = R_alpha_mean_abs_P6;
all_powers_abs{1, 2} = R_theta_mean_abs_P1;
all_powers_abs{2, 2} = R_theta_mean_abs_P3;
all_powers_abs{3, 2} = R_theta_mean_abs_P6;

all_powers = cell(3, 2); 
all_powers{1, 1} = R_alpha_mean_P1;
all_powers{2, 1} = R_alpha_mean_P3;
all_powers{3, 1} = R_alpha_mean_P6;
all_powers{1, 2} = R_theta_mean_P1;
all_powers{2, 2} = R_theta_mean_P3;
all_powers{3, 2} = R_theta_mean_P6;


X = [lagVec, fliplr(lagVec)];
ylabels = {'magnitude |r(Lag)|', 'signed r(Lag)'};
Subject_sigend_R = cell(length(Subject_list), 3, 2);
RR = {R_abs, R};
for power = 1:2 % theta and alpha

    for i = 1:3

        subplot(2, 3, (power-1)*3 + i); hold on;

        % absolute r(Lag) or signed r(Lag)?
        % signal = all_powers_abs{i, power};
        signal = all_powers{i, power};
        ylabel_s = ylabels{2};
        RRR = RR{2};


        % Y = [mean(signal, 1) + ...
        %        std(signal, [], 1)./sqrt(size(signal,1 )), ...
        %      fliplr(mean(signal, 1) - ...
        %      std(signal, [], 1)./sqrt(size(signal,1 )))];
        Y = [mean(signal, 1) - std(signal, [], 1), ...
             fliplr(mean(signal, 1) + std(signal, [], 1))];
        
        fill(X, Y, colors(i, :), 'EdgeColor', 'none', 'FaceColor', colors(i, :), ...
            'FaceAlpha', 0.2, 'HandleVisibility', 'off')
        plot(lagVec, mean(signal, 1), 'Color', colors(i, :), ...
            'LineWidth', 2, 'HandleVisibility', 'off')
        xline(0, 'LineStyle', '--', 'HandleVisibility', 'off')
        xlim([lagVec(1) lagVec(end)])
        
        
        
        Subject_sigend_R = cell(length(Subject_list), 3, 2);
        for sub = 1:length(Subject_list)
    
            if isempty(Main_data_structure_trial{sub, 3}), continue, end;
    
            [maxR maxLagSamples] = max(RRR{sub,i, power}, [], 2);
            sigend_R = cellfun(@(x1, x2) x1(1, x2), ...
                num2cell(R{sub, i, power}, 2), num2cell(maxLagSamples));
            Subject_sigend_R{sub, i, power} = sigend_R;
            
            meanLag = mean(lagVec(maxLagSamples));
            stdLag = std(lagVec(maxLagSamples), []);
            semLag = std(lagVec(maxLagSamples), [])/sqrt(numel(maxLagSamples));
            meanR = mean(maxR);
            stdR = std(maxR, []);
            semR = std(maxR, [])/sqrt(numel(maxLagSamples));
            
           
            % errorbar(meanLag, meanR, semR, semR, semLag, semLag, ...
            %     "Color", 0.7*[1 1 1], 'HandleVisibility', 'off')
            errorbar(meanLag, meanR, stdR, stdR, stdLag, stdLag, ...
                "Color", 0.7*[1 1 1], 'HandleVisibility', 'off')
            scatter(meanLag, meanR, 40, 'Marker', markerList{sub}, ...
                'LineWidth', 2, 'MarkerFaceColor', 'w', ...
                'MarkerEdgeColor', colors_markers(sub, :));
        end
        ylim([0 1])

        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, ...
            'XTickLabelRotation', 45);
        set(gca, 'FontSize', 16)

        if power == 2
            xlh = xlabel('Lag (% Cycle)', 'FontWeight', 'bold', ...
                'FontSize', 18);
            xlh.Position(2) = -0.15;
        else
            tth = title([titles{i}, ' Pressure'], ...
                'FontSize', 18, 'FontWeight', 'bold');
            tth.Position(2) = 1.1;
        end

        if i == 1
            if power == 1
                ylh = ylabel(sprintf(['\\bf\\alpha Power vs. Tracking Error\n', ...
                              '\\bfCorrelation ' ylabel_s]), ...
                     'Interpreter','tex', ...      % <-- important
                     'FontSize',18);
                ylh.Position(1) = lagVec(1) - 30;
            else
                ylh = ylabel(sprintf(['\\bf\\theta Power vs. Tracking Error\n', ...
                              '\\bfCorrelation ' ylabel_s]), ...
                     'Interpreter','tex', ...      % <-- important
                     'FontSize',18);
                ylh.Position(1) = lagVec(1) - 30;
            end
        end

        if power == 2 && i == 3
            legend(legends, 'FontName', 'Arial', 'FontSize', 12, ...
                'Box', 'off')
        end

    end
    
    

end




%% Option 2B: 
%  Cross-Correlation on actual signal (no subtraction of mean signal)
%  EEG is baseline-corrected power 

freqs = Main_data_structure_trial{end, 3};
idx_theta = (freqs > 4) & (freqs < 8);
idx_alpha = (freqs > 8) & (freqs < 14);

nTime = length(new_times);
timeIdx = 1:nTime;              % 1..nTime = 0..100% cycle
maxLagPct    = 50;              % max lag in % of cycle (e.g. ±10%)
maxLagSamples = round(maxLagPct/100 * nTime);
lagVec       = -maxLagSamples:maxLagSamples;
nLags        = numel(lagVec);


nCond = 3;
R_theta = cell(length(Subject_list), nCond);
R_alpha = cell(length(Subject_list), nCond);
for sub = 1:length(Subject_list)
    
    if isempty(Main_data_structure_trial{sub, 3}), continue, end;

    % Tracking Error
    Err = Main_data_structure_trial{sub, 3}.Err;
    % Theta and Alpha Power
    EEG_theta = mean(Main_data_structure_trial_baseln{sub, 3}(idx_theta, :, :), 1);
    EEG_theta = squeeze(EEG_theta)';
    EEG_alpha = mean(Main_data_structure_trial_baseln{sub, 3}(idx_alpha, :, :), 1);
    EEG_alpha = squeeze(EEG_alpha)';

    pressure = Main_data_structure_trial{sub, 3}.Pressure;
    for c = 1:nCond
        
        pressure_idx = pressure == Conditions(c);
        curErr = Err(pressure_idx, :);               % [nTrials_sc x nTime]
        curPow_theta = EEG_theta(pressure_idx, :);   % [nTrials_sc x nTime]
        curPow_alpha = EEG_alpha(pressure_idx, :);   % [nTrials_sc x nTime]
        

        nTrials_sc = size(curErr, 1);
        nTime      = size(curErr, 2);
        
       
        rLags_trials_theta = nan(nTrials_sc, nLags);
        for k = 1:nTrials_sc
            e = curErr(k, :);       % [1 x nTime]
            p_theta = curPow_theta(k, :);
            p_alpha = curPow_alpha(k, :);

            % e = e - mean(e, 'omitnan');
            % p = p - mean(p, 'omitnan');

            rLags_trials_theta(k, :) = xcorr(e, p_theta, maxLagSamples, "normalized");
            rLags_trials_alpha(k, :) = xcorr(e, p_alpha, maxLagSamples, "normalized");
        
        end



        % e = mean(errRes, 1);
        % p = mean(powRes_theta, 1);
        % rLags_trials = xcorr(e, p, maxLagSamples, "normalized");

        R_theta{sub, c} = rLags_trials_theta;
        R_alpha{sub, c} = rLags_trials_alpha;
    end
end



%% plot the Cross-Correlation results

% Get monitor information
monitors = get(0, 'MonitorPositions');

fig = figure('name', ...
    ['Cross-Correlation of Baseline-Corrected EEG power', ...
    ' vs Tracking Error actual signals'], ...
    'InvertHardcopy', 'off', 'PaperType', 'a2', ...
    'PaperOrientation', 'landscape', ...
    'Resize', 'off');

% For second monitor (row 2), add drawnow before setting position
drawnow;  % Let MATLAB finish drawing on primary monitor first
pause(0.1);  % Short pause helps

% Now set position on second monitor
% set(gcf, 'Position', monitors(1, :));  % Use entire second monitor
set(fig, 'Position', [monitors(1,1)-200, monitors(1,2)+600, 1200, 800]);  % Use entire second monitor


% 13 built-in marker types
markerList = {'o','s','d','^','v','>','<','p','h','s','^','v','p','o'};
% colors_markers = lines(length(Subject_list));   % 14 distinct colors
colors_markers = [
    0.00 0.45 0.74  % blue
    0.85 0.33 0.10  % reddish orange
    0.93 0.69 0.13  % yellow-ish
    0.49 0.18 0.56  % purple
    0.47 0.67 0.19  % green
    0.30 0.75 0.93  % cyan
    0.64 0.08 0.18  % dark red
    0.13 0.35 0.60  % dark blue
    0.76 0.49 0.27  % brown-ish
    0.55 0.57 0.00  % olive
    0.65 0.70 0.96  % light bluish
    0.49 0.39 0.64  % desaturated purple
    0.70 0.67 0.18  % mustard
    0.40 0.40 0.40  % gray
];




% theta band
R_theta_abs = cellfun(@(x) abs(x), R_theta, 'UniformOutput', false);
R_theta_mean_abs = R_theta_abs;
R_theta_mean_abs(10, :) = [];
R_theta_mean_abs = cellfun(@(x) mean(x, 1), R_theta_mean_abs, 'UniformOutput', false);
R_theta_mean_abs_P1 = cell2mat(R_theta_mean_abs(:, 1));
R_theta_mean_abs_P3 = cell2mat(R_theta_mean_abs(:, 2));
R_theta_mean_abs_P6 = cell2mat(R_theta_mean_abs(:, 3));

R_theta_mean = R_theta;
R_theta_mean(10, :) = [];
R_theta_mean = cellfun(@(x) mean(x, 1), R_theta_mean, 'UniformOutput', false);
R_theta_mean_P1 = cell2mat(R_theta_mean(:, 1));
R_theta_mean_P3 = cell2mat(R_theta_mean(:, 2));
R_theta_mean_P6 = cell2mat(R_theta_mean(:, 3));

% alpha band
R_alpha_abs = cellfun(@(x) abs(x), R_alpha, 'UniformOutput', false);
R_alpha_mean_abs = R_alpha_abs;
R_alpha_mean_abs(10, :) = [];
R_alpha_mean_abs = cellfun(@(x) mean(x, 1), R_alpha_mean_abs, 'UniformOutput', false);
R_alpha_mean_abs_P1 = cell2mat(R_alpha_mean_abs(:, 1));
R_alpha_mean_abs_P3 = cell2mat(R_alpha_mean_abs(:, 2));
R_alpha_mean_abs_P6 = cell2mat(R_alpha_mean_abs(:, 3));

R_alpha_mean = R_alpha;
R_alpha_mean(10, :) = [];
R_alpha_mean = cellfun(@(x) mean(x, 1), R_alpha_mean, 'UniformOutput', false);
R_alpha_mean_P1 = cell2mat(R_alpha_mean(:, 1));
R_alpha_mean_P3 = cell2mat(R_alpha_mean(:, 2));
R_alpha_mean_P6 = cell2mat(R_alpha_mean(:, 3));

R = cat(3, R_alpha, R_theta);
R_abs = cat(3, R_alpha_abs, R_theta_abs);

% three rows for three conditions, 
% two coloums for alpha and theta
all_powers_abs = cell(3, 2); 
all_powers_abs{1, 1} = R_alpha_mean_abs_P1;
all_powers_abs{2, 1} = R_alpha_mean_abs_P3;
all_powers_abs{3, 1} = R_alpha_mean_abs_P6;
all_powers_abs{1, 2} = R_theta_mean_abs_P1;
all_powers_abs{2, 2} = R_theta_mean_abs_P3;
all_powers_abs{3, 2} = R_theta_mean_abs_P6;

all_powers = cell(3, 2); 
all_powers{1, 1} = R_alpha_mean_P1;
all_powers{2, 1} = R_alpha_mean_P3;
all_powers{3, 1} = R_alpha_mean_P6;
all_powers{1, 2} = R_theta_mean_P1;
all_powers{2, 2} = R_theta_mean_P3;
all_powers{3, 2} = R_theta_mean_P6;


X = [lagVec, fliplr(lagVec)];
ylabels = {'magnitude |r(Lag)|', 'signed r(Lag)'};
Subject_sigend_R = cell(length(Subject_list), 3, 2);
RR = {R_abs, R};
for power = 1:2 % theta and alpha

    for i = 1:3

        subplot(2, 3, (power-1)*3 + i); hold on;

        % absolute r(Lag) or signed r(Lag)?
        % signal = all_powers_abs{i, power};
        signal = all_powers{i, power};
        ylabel_s = ylabels{2};
        RRR = RR{2};


        % Y = [mean(signal, 1) + ...
        %        std(signal, [], 1)./sqrt(size(signal,1 )), ...
        %      fliplr(mean(signal, 1) - ...
        %      std(signal, [], 1)./sqrt(size(signal,1 )))];
        Y = [mean(signal, 1) - std(signal, [], 1), ...
             fliplr(mean(signal, 1) + std(signal, [], 1))];
        
        fill(X, Y, colors(i, :), 'EdgeColor', 'none', 'FaceColor', colors(i, :), ...
            'FaceAlpha', 0.2, 'HandleVisibility', 'off')
        plot(lagVec, mean(signal, 1), 'Color', colors(i, :), ...
            'LineWidth', 2, 'HandleVisibility', 'off')
        xline(0, 'LineStyle', '--', 'HandleVisibility', 'off')
        xlim([lagVec(1) lagVec(end)])
        
        
        
        Subject_sigend_R = cell(length(Subject_list), 3, 2);
        for sub = 1:length(Subject_list)
    
            if isempty(Main_data_structure_trial{sub, 3}), continue, end;
    
            [maxR maxLagSamples] = max(RRR{sub,i, power}, [], 2);
            sigend_R = cellfun(@(x1, x2) x1(1, x2), ...
                num2cell(R{sub, i, power}, 2), num2cell(maxLagSamples));
            Subject_sigend_R{sub, i, power} = sigend_R;
            
            meanLag = mean(lagVec(maxLagSamples));
            stdLag = std(lagVec(maxLagSamples), []);
            semLag = std(lagVec(maxLagSamples), [])/sqrt(numel(maxLagSamples));
            meanR = mean(maxR);
            stdR = std(maxR, []);
            semR = std(maxR, [])/sqrt(numel(maxLagSamples));
            
           
            % errorbar(meanLag, meanR, semR, semR, semLag, semLag, ...
            %     "Color", 0.7*[1 1 1], 'HandleVisibility', 'off')
            errorbar(meanLag, meanR, stdR, stdR, stdLag, stdLag, ...
                "Color", 0.7*[1 1 1], 'HandleVisibility', 'off')
            scatter(meanLag, meanR, 40, 'Marker', markerList{sub}, ...
                'LineWidth', 2, 'MarkerFaceColor', 'w', ...
                'MarkerEdgeColor', colors_markers(sub, :));
        end
        ylim([0 1])

        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, ...
            'XTickLabelRotation', 45);
        set(gca, 'FontSize', 16)

        if power == 2
            xlh = xlabel('Lag (% Cycle)', 'FontWeight', 'bold', ...
                'FontSize', 18);
            xlh.Position(2) = -0.15;
        else
            tth = title([titles{i}, ' Pressure'], ...
                'FontSize', 18, 'FontWeight', 'bold');
            tth.Position(2) = 1.1;
        end

        if i == 1
            if power == 1
                ylh = ylabel(sprintf(['\\bf\\alpha Power vs. Tracking Error\n', ...
                              '\\bfCorrelation ' ylabel_s]), ...
                     'Interpreter','tex', ...      % <-- important
                     'FontSize',18);
                ylh.Position(1) = lagVec(1) - 30;
            else
                ylh = ylabel(sprintf(['\\bf\\theta Power vs. Tracking Error\n', ...
                              '\\bfCorrelation ' ylabel_s]), ...
                     'Interpreter','tex', ...      % <-- important
                     'FontSize',18);
                ylh.Position(1) = lagVec(1) - 30;
            end
        end

        if power == 2 && i == 3
            legend(legends, 'FontName', 'Arial', 'FontSize', 12, ...
                'Box', 'off')
        end

    end
    
    

end



























