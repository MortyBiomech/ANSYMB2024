% Main EMG processing
clc
clear

%% Add and Define Necessary Paths
addpath(genpath('C:\Morteza\MyProjects\ANSYMB2024')); % main folder containing all codes and data
addpath(genpath('C:\Morteza\LSL\xdf-Matlab-master')); % required for loading the XDF files.

% Change path to the directory on your PC which raw XDF files are stored:
data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
rawdata_path = [data_path, '0_source_data\'];


%% Colors for plotting the signals
All_colours = struct('dark_blue', [0, 0.4470, 0.7410], 'light_blue', [0.5010, 0.8450, 0.9930], ...
        'dark_orange', [0.8500, 0.3250, 0.0980], 'light_orange', [0.9690, 0.7940, 0.3250], ...
        'dark_green', [0.4660, 0.6740, 0.1880], 'light_green', [0.7960, 0.9740, 0.7410]);




%% Load data
subject_id = 8;

load([data_path, '6_Trials_Info_and_Epoched_data', filesep, 'sub-', ...
    num2str(subject_id), filesep,'Trials_Info.mat']);

condition = 1;
switch condition
    case 1
        epoch_base = 'Flexion epochs';
        load([data_path, '6_Trials_Info_and_Epoched_data', filesep, 'sub-', ...
            num2str(subject_id), filesep,'Epochs_Flexion_based.mat']);
        data = Epochs_Flexion_based;
    case 2
        epoch_base = 'Extension epochs';
        load([data_path, '6_Trials_Info_and_Epoched_data', filesep, 'sub-', ...
            num2str(subject_id), filesep,'Epochs_Extension_based.mat']);
        data = Epochs_Extension_based;
end

Names = data{1,1}.EMG_stream.Names;
Names = cellfun(@(x) strrep(x, '_', ' '), Names, 'UniformOutput', false);

%% Initialize a cell for saving length-normalized EMG data
EMG_Preprocessed = repmat({struct('without_outlier_removal', struct('value', [], 'RMS', []), ...
    'with_outlier_removal', [])}, 1, length(Trials_Info));


%% Select the trials with the same Pressure
bad_trials_path = [data_path, '6_0_Trials_Info_and_Events', ...
    filesep, 'sub-', num2str(subject_id)];

b = load(fullfile(bad_trials_path, 'bad_trials_EEG_based.mat'));
bad_trials_EEG_based = b.bad_trials_EEG_based;

P1 = struct();
P1.trials = [];
P1.scores = [];

P3 = struct();
P3.trials = [];
P3.scores = [];

P6 = struct();
P6.trials = [];
P6.scores = [];

for i = 1:length(Trials_Info)
    if ~ismember(i, bad_trials_EEG_based)
        p = Trials_Info{1, i}.General.Pressure;
        switch p
            case 1
                P1.trials(1, end+1) = i;
                P1.scores(1, end+1) = Trials_Info{1, i}.General.Score;
            case 3
                P3.trials(1, end+1) = i;
                P3.scores(1, end+1) = Trials_Info{1, i}.General.Score;
            case 6
                P6.trials(1, end+1) = i;
                P6.scores(1, end+1) = Trials_Info{1, i}.General.Score;
        end
    end
end



%% Lenght Normalization - without outlier removal
%%% P1
L_P1 = zeros(1, length(P1.trials));
for i = 1:length(P1.trials)
    l = cellfun("length", data{1, P1.trials(i)}.EMG_stream.Sensors_Preprocessed);
    L_P1(1, i) = max(l);
end
[max_L_P1, ~] = max(L_P1);

for i = 1:length(P1.trials)
    Q = length(data{1, P1.trials(i)}.EMG_stream.Sensors_Preprocessed);
    for j = 1:Q
        x_old = data{1, P1.trials(i)}.EMG_stream.Times{1, j};
        x_new = linspace(x_old(1), x_old(end), max_L_P1);
        y_old = data{1, P1.trials(i)}.EMG_stream.Sensors_Preprocessed{1, j};
        y_new = interp1(x_old', y_old', x_new', "spline");
    
        EMG_Preprocessed{1, P1.trials(i)}.without_outlier_removal.value(:, :, end + 1) = y_new';
    end
    EMG_Preprocessed{1, P1.trials(i)}.without_outlier_removal.value(:, :, 1) = [];
end

%%% P3
L_P3 = zeros(1, length(P3.trials));
for i = 1:length(P3.trials)
    l = cellfun("length", data{1, P3.trials(i)}.EMG_stream.Sensors_Preprocessed);
    L_P3(1, i) = max(l);
end
[max_L_P3, ~] = max(L_P3);

for i = 1:length(P3.trials)
    Q = length(data{1, P3.trials(i)}.EMG_stream.Sensors_Preprocessed);
    for j = 1:Q
        x_old = data{1, P3.trials(i)}.EMG_stream.Times{1, j};
        x_new = linspace(x_old(1), x_old(end), max_L_P3);
        y_old = data{1, P3.trials(i)}.EMG_stream.Sensors_Preprocessed{1, j};
        y_new = interp1(x_old', y_old', x_new', "spline");
    
        EMG_Preprocessed{1, P3.trials(i)}.without_outlier_removal.value(:, :, end + 1) = y_new';
    end
    EMG_Preprocessed{1, P3.trials(i)}.without_outlier_removal.value(:, :, 1) = [];
end

%%% P6
L_P6 = zeros(1, length(P6.trials));
for i = 1:length(P6.trials)
    l = cellfun("length", data{1, P6.trials(i)}.EMG_stream.Sensors_Preprocessed);
    L_P6(1, i) = max(l);
end
[max_L_P6, ~] = max(L_P6);

for i = 1:length(P6.trials)
    Q = length(data{1, P6.trials(i)}.EMG_stream.Sensors_Preprocessed);
    for j = 1:Q
        x_old = data{1, P6.trials(i)}.EMG_stream.Times{1, j};
        x_new = linspace(x_old(1), x_old(end), max_L_P6);
        y_old = data{1, P6.trials(i)}.EMG_stream.Sensors_Preprocessed{1, j};
        y_new = interp1(x_old', y_old', x_new', "spline");
    
        EMG_Preprocessed{1, P6.trials(i)}.without_outlier_removal.value(:, :, end + 1) = y_new';
    end
    EMG_Preprocessed{1, P6.trials(i)}.without_outlier_removal.value(:, :, 1) = [];
end


%% All the EMGs for one pressure condition together
EMG_P1 = [];
EMG_P1_trial_id = [];
for i = 1:length(P1.trials)
    EMG_P1 = cat(3, EMG_P1, EMG_Preprocessed{1, P1.trials(i)}.without_outlier_removal.value);
    EMG_P1_trial_id = cat(2, EMG_P1_trial_id, ...
        repmat(P1.trials(i), 1, size(EMG_Preprocessed{1, P1.trials(i)}.without_outlier_removal.value, 3)));
end

EMG_P3 = [];
EMG_P3_trial_id = [];
for i = 1:length(P3.trials)
    EMG_P3 = cat(3, EMG_P3, EMG_Preprocessed{1, P3.trials(i)}.without_outlier_removal.value);
    EMG_P3_trial_id = cat(2, EMG_P3_trial_id, ...
        repmat(P3.trials(i), 1, size(EMG_Preprocessed{1, P3.trials(i)}.without_outlier_removal.value, 3)));
end

EMG_P6 = [];
EMG_P6_trial_id = [];
for i = 1:length(P6.trials)
    EMG_P6 = cat(3, EMG_P6, EMG_Preprocessed{1, P6.trials(i)}.without_outlier_removal.value);
    EMG_P6_trial_id = cat(2, EMG_P6_trial_id, ...
        repmat(P6.trials(i), 1, size(EMG_Preprocessed{1, P6.trials(i)}.without_outlier_removal.value, 3)));
end



%% Plot EMG signals before outlier removal
X_P1 = linspace(0, 100, size(EMG_P1, 2));
X_P3 = linspace(0, 100, size(EMG_P3, 2));
X_P6 = linspace(0, 100, size(EMG_P6, 2));

%%% P1
EMG_P1_median = median(EMG_P1, 3);
figure();
t = tiledlayout(3,4);
title(t, [epoch_base, ', Pressure 1 bar, Before outlier removal'])
for i = 1:size(EMG_P1, 1)
    nexttile; hold on;
    title(Names{i})
    xlabel('Cycle [%]')
    ylabel('EMG Activity [\muV]')
    h1 = plot(X_P1, 1e3*squeeze(EMG_P1(i, :, :))', 'Color', All_colours.light_blue, 'LineWidth', 0.5);
    h2 = plot(X_P1, 1e3*mean(EMG_P1(i, :, :), 3), 'Color', All_colours.dark_blue, 'LineWidth', 2);
    h3 = plot(X_P1, 1e3*EMG_P1_median(i, :), 'Color', [0.5, 0, 0.5], 'LineWidth', 2);
    
    hold off;
end
lgd = legend([h1(1), h2, h3], {'All epochs', 'mean', 'median'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom

%%% P3
EMG_P3_median = median(EMG_P3, 3);
figure();
t = tiledlayout(3,4);
title(t, [epoch_base, ', Pressure 3 bar, Before outlier removal'])
for i = 1:size(EMG_P1, 1)
    nexttile; hold on;
    title(Names{i})
    xlabel('Cycle [%]')
    ylabel('EMG Activity [\muV]')
    h1 = plot(X_P3, 1e3*squeeze(EMG_P3(i, :, :))', 'Color', All_colours.light_orange, 'LineWidth', 0.5);
    h2 = plot(X_P3, 1e3*mean(EMG_P3(i, :, :), 3), 'Color', All_colours.dark_orange, 'LineWidth', 2);
    h3 = plot(X_P3, 1e3*EMG_P3_median(i, :), 'Color', [0.5, 0, 0.5], 'LineWidth', 2);
    
    hold off;
end
lgd = legend([h1(1), h2, h3], {'All epochs', 'mean', 'median'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom

%%% P6
EMG_P6_median = median(EMG_P6, 3);
figure();
t = tiledlayout(3,4);
title(t, [epoch_base, ', Pressure 6 bar, Before outlier removal'])
for i = 1:size(EMG_P1, 1)
    nexttile; hold on;
    title(Names{i})
    xlabel('Cycle [%]')
    ylabel('EMG Activity [\muV]')
    h1 = plot(X_P6, 1e3*squeeze(EMG_P6(i, :, :))', 'Color', All_colours.light_green, 'LineWidth', 0.5);
    h2 = plot(X_P6, 1e3*mean(EMG_P6(i, :, :), 3), 'Color', All_colours.dark_green, 'LineWidth', 2);
    h3 = plot(X_P6, 1e3*EMG_P6_median(i, :), 'Color', [0.5, 0, 0.5], 'LineWidth', 2);

    hold off;
end
lgd = legend([h1(1), h2, h3], {'All epochs', 'mean', 'median'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom


%% Compute medians and our custom-made error to find most similar signals
%%% P1
err_P1_flx = zeros(size(EMG_P1, 1), size(EMG_P1, 3), 2); % In 3rd dimension the trial-id is stored.
for m = 1:size(EMG_P1, 1)
    for i = 1:size(err_P1_flx, 2)
        err_P1_flx(m, i, 1) = sum((1 + abs(EMG_P1(m, :, i) - EMG_P1_median(m, :)) ).^6);
        err_P1_flx(m, i, 2) = EMG_P1_trial_id(i);
    end
end

% sort error values and keep the trial-ids after sorting
err_P1_flx_sorted = zeros(size(err_P1_flx));
EMG_P1_err_sorted = zeros(size(EMG_P1));
for i = 1:size(err_P1_flx, 1) % Iterate over the first dimension (muscles)
    % Sort the values in the first slice of the third dimension
    [sorted_values, sort_idx] = sort(err_P1_flx(i,:,1), 2); 
    
    % Rearrange the first slice with the sorted values
    err_P1_flx_sorted(i,:,1) = sorted_values;
    
    % Use the same sorting indices to rearrange the second slice
    err_P1_flx_sorted(i,:,2) = err_P1_flx(i, sort_idx, 2);

    % Use the sorting indices to rearrange the EMG_P6
    EMG_P1_err_sorted(i, :, :) = EMG_P1(i, :, sort_idx);
end


%%% P3
err_P3_flx = zeros(size(EMG_P3, 1), size(EMG_P3, 3), 2); % In 3rd dimension the trial-id is stored.
for m = 1:size(EMG_P3, 1)
    for i = 1:size(err_P3_flx, 2)
        err_P3_flx(m, i, 1) = sum((1 + abs(EMG_P3(m, :, i) - EMG_P3_median(m, :)) ).^6);
        err_P3_flx(m, i, 2) = EMG_P3_trial_id(i);
    end
end

% sort error values and keep the trial-ids after sorting
err_P3_flx_sorted = zeros(size(err_P3_flx));
EMG_P3_err_sorted = zeros(size(EMG_P3));
for i = 1:size(err_P3_flx, 1) % Iterate over the first dimension (muscles)
    % Sort the values in the first slice of the third dimension
    [sorted_values, sort_idx] = sort(err_P3_flx(i,:,1), 2); 
    
    % Rearrange the first slice with the sorted values
    err_P3_flx_sorted(i,:,1) = sorted_values;
    
    % Use the same sorting indices to rearrange the second slice
    err_P3_flx_sorted(i,:,2) = err_P3_flx(i, sort_idx, 2);

    % Use the sorting indices to rearrange the EMG_P6
    EMG_P3_err_sorted(i, :, :) = EMG_P3(i, :, sort_idx);
end


%%% P6
err_P6_flx = zeros(size(EMG_P6, 1), size(EMG_P6, 3), 2); % In 3rd dimension the trial-id is stored.
for m = 1:size(EMG_P6, 1)
    for i = 1:size(err_P6_flx, 2)
        err_P6_flx(m, i, 1) = sum((1 + abs(EMG_P6(m, :, i) - EMG_P6_median(m, :)) ).^6);
        err_P6_flx(m, i, 2) = EMG_P6_trial_id(i);
    end
end

% sort error values and keep the trial-ids after sorting
err_P6_flx_sorted = zeros(size(err_P6_flx));
EMG_P6_err_sorted = zeros(size(EMG_P6));
for i = 1:size(err_P6_flx, 1) % Iterate over the first dimension (muscles)
    % Sort the values in the first slice of the third dimension
    [sorted_values, sort_idx] = sort(err_P6_flx(i,:,1), 2); 
    
    % Rearrange the first slice with the sorted values
    err_P6_flx_sorted(i,:,1) = sorted_values;
    
    % Use the same sorting indices to rearrange the second slice
    err_P6_flx_sorted(i,:,2) = err_P6_flx(i, sort_idx, 2);

    % Use the sorting indices to rearrange the EMG_P6
    EMG_P6_err_sorted(i, :, :) = EMG_P6(i, :, sort_idx);
end


%% Now let's decide how many epochs we should keep
% we keep 25, 50, 75, and 90% of epochs in each muscle with lower errors
% and then we look the effects of our selection on average signals.

P = 0.6; % keep P*100 percent of epochs

K = floor(P*size(EMG_P1, 3));
EMG_P1_selected = zeros(size(EMG_P1, 1), size(EMG_P1, 2), K);
for i = 1:size(EMG_P1, 1)
    EMG_P1_selected(i, :, :) = EMG_P1_err_sorted(i, :, 1:K);
end

K = floor(P*size(EMG_P3, 3));
EMG_P3_selected = zeros(size(EMG_P3, 1), size(EMG_P3, 2), K);
for i = 1:size(EMG_P3, 1)
    EMG_P3_selected(i, :, :) = EMG_P3_err_sorted(i, :, 1:K);
end

K = floor(P*size(EMG_P6, 3));
EMG_P6_selected = zeros(size(EMG_P6, 1), size(EMG_P6, 2), K);
for i = 1:size(EMG_P6, 1)
    EMG_P6_selected(i, :, :) = EMG_P6_err_sorted(i, :, 1:K);
end


%% Plot Selected/Removed Epochs and see the effect on the total average

figure();
tiledlayout(2, 2)
sgtitle(['Effect of removing ', num2str(100*(1-P)), ...
    '% of epochs (', num2str(size(EMG_P1, 3) - K), ' out of ', ...
    num2str(size(EMG_P1, 3)), ') with high error values - Pressure 1 bar, ', epoch_base]);
for i = 1:size(EMG_P1, 1)-6
    nexttile; hold on;

    h3 = plot(X_P1, 1e3*squeeze(EMG_P1_err_sorted(i, :, K+1:end))', 'Color', 0.5*[1, 1, 1], 'LineWidth', 0.5);

    h1 = plot(X_P1, 1e3*squeeze(EMG_P1_selected(i, :, :))', 'Color', All_colours.light_blue, 'LineWidth', 0.5, 'LineStyle', '--');

    h2 = plot(X_P1, 1e3*mean(EMG_P1_selected(i, :, :), 3), 'Color', All_colours.dark_blue, 'LineWidth', 4);

    h4 = plot(X_P1, 1e3*mean(EMG_P1(i, :, :), 3), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');

    h5 = plot(X_P1, 1e3*EMG_P1_median(i, :), 'Color', [0.5, 0, 0.5], 'LineWidth', 2, 'LineStyle', '--');
    
    hold off;

    title(Names{i})
    xlabel('Cycle [%]')
    ylabel('EMG Activity [\muV]')
end

% Create a legend for the entire tiled layout
lgd = legend([h1(1), h2, h3(1), h4, h5], ...
    {'Selected Epochs', 'Mean of Selected Epochs', 'Removed Epochs', 'Mean of All Epoch', 'Median of All Epochs'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom



figure();
tiledlayout(2, 2)
sgtitle(['Effect of removing ', num2str(100*(1-P)), ...
    '% of epochs (', num2str(size(EMG_P3, 3) - K), ' out of ', ...
    num2str(size(EMG_P3, 3)), ') with high error values - Pressure 3 bar, ', epoch_base]);
for i = 1:size(EMG_P3, 1)-6
    nexttile; hold on;

    h3 = plot(X_P3, 1e3*squeeze(EMG_P3_err_sorted(i, :, K+1:end))', 'Color', 0.5*[1, 1, 1], 'LineWidth', 0.5);

    h1 = plot(X_P3, 1e3*squeeze(EMG_P3_selected(i, :, :))', 'Color', All_colours.light_orange, 'LineWidth', 0.5, 'LineStyle', '--');

    h2 = plot(X_P3, 1e3*mean(EMG_P3_selected(i, :, :), 3), 'Color', All_colours.dark_orange, 'LineWidth', 4);

    h4 = plot(X_P3, 1e3*mean(EMG_P3(i, :, :), 3), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');

    h5 = plot(X_P3, 1e3*EMG_P3_median(i, :), 'Color', [0.5, 0, 0.5], 'LineWidth', 2, 'LineStyle', '--');
    
    hold off;

    title(Names{i})
    xlabel('Cycle [%]')
    ylabel('EMG Activity [\muV]')
end

% Create a legend for the entire tiled layout
lgd = legend([h1(1), h2, h3(1), h4, h5], ...
    {'Selected Epochs', 'Mean of Selected Epochs', 'Removed Epochs', 'Mean of All Epoch', 'Median of All Epochs'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom



figure();
tiledlayout(2, 2)
sgtitle(['Effect of removing ', num2str(100*(1-P)), ...
    '% of epochs (', num2str(size(EMG_P6, 3) - K), ' out of ', ...
    num2str(size(EMG_P6, 3)), ') with high error values - Pressure 6 bar, ', epoch_base]);
for i = 1:size(EMG_P6, 1)-6
    nexttile; hold on;

    h3 = plot(X_P6, 1e3*squeeze(EMG_P6_err_sorted(i, :, K+1:end))', 'Color', 0.5*[1, 1, 1], 'LineWidth', 0.5);

    h1 = plot(X_P6, 1e3*squeeze(EMG_P6_selected(i, :, :))', 'Color', All_colours.light_green, 'LineWidth', 0.5, 'LineStyle', '--');

    h2 = plot(X_P6, 1e3*mean(EMG_P6_selected(i, :, :), 3), 'Color', All_colours.dark_green, 'LineWidth', 4);

    h4 = plot(X_P6, 1e3*mean(EMG_P6(i, :, :), 3), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');

    h5 = plot(X_P6, 1e3*EMG_P6_median(i, :), 'Color', [0.5, 0, 0.5], 'LineWidth', 2, 'LineStyle', '--');
    
    hold off;

    title(Names{i})
    xlabel('Cycle [%]')
    ylabel('EMG Activity [\muV]')
end

% Create a legend for the entire tiled layout
lgd = legend([h1(1), h2, h3(1), h4, h5], ...
    {'Selected Epochs', 'Mean of Selected Epochs', 'Removed Epochs', 'Mean of All Epoch', 'Median of All Epochs'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom


%% look deeper into the removed epochs

P = 0.9; % keep P*100 percent of epochs

[epoch_numbers, removed_epochs_total] = ...
    removed_epochs_count(data, P, P1, P3, P6, ...
    err_P1_flx_sorted, err_P3_flx_sorted, err_P6_flx_sorted);
epoch_numbers(bad_trials_EEG_based) = [];
removed_epochs_total(:, bad_trials_EEG_based) = [];


%% plot number of removed epochs

% define the colors for the bar plot
colors = zeros(length(data), 3);
colors(P1.trials, :) = repmat(All_colours.light_blue, numel(P1.trials), 1);
colors(P3.trials, :) = repmat(All_colours.light_orange, numel(P3.trials), 1);
colors(P6.trials, :) = repmat(All_colours.light_green, numel(P6.trials), 1);
colors(bad_trials_EEG_based, :) = [];

%% main plot
h = tiledlayout(4,3);
title(h, epoch_base)
x = 1:length(epoch_numbers);
for i = 1:12
    switch i
        case 1
            muscle_m = 1;
        case 4
            muscle_m = 2;
        case 7
            muscle_m = 3;
        case 10
            muscle_m = 4;
    end

    nexttile(i);

    if mod(i, 3) == 1
        P = 0.9;
    elseif mod(i, 3) == 2
        P = 0.8;
    else
        P = 0.6;
    end

    [epoch_numbers, removed_epochs_total] = ...
        removed_epochs_count(data, P, P1, P3, P6, ...
        err_P1_flx_sorted, err_P3_flx_sorted, err_P6_flx_sorted);

    epoch_numbers(bad_trials_EEG_based) = [];
    removed_epochs_total(:, bad_trials_EEG_based) = [];

    h1 = bar(x, epoch_numbers, 'EdgeColor', 'none'); 
    xlabel('Trials')
    ylabel('Epoch Count')
    h1.FaceColor = 'flat';
    h1.CData = colors;
    
    hold on;
    
    h2 = bar(x, removed_epochs_total(muscle_m, :), 'EdgeColor', 'none');
    h2.FaceColor = 'flat';
    h2.CData = repmat(0.2*[1 1 1], numel(removed_epochs_total(1,:)), 1);
    
    xlim([1, numel(epoch_numbers)])
    ylim([0 max(epoch_numbers)+1])
    title([Names{muscle_m}, ', ', num2str(100*(1-P)), '% removal '])
end


%% Look at the removed epochs in detail in each pressure condition
% define the colors for the bar plot without removing bad-trials ids 
colors = zeros(length(data), 3);
colors(P1.trials, :) = repmat(All_colours.light_blue, numel(P1.trials), 1);
colors(P3.trials, :) = repmat(All_colours.light_orange, numel(P3.trials), 1);
colors(P6.trials, :) = repmat(All_colours.light_green, numel(P6.trials), 1);

% P1
figure();
h = tiledlayout(4,3);
title(h, [epoch_base, ', Pressure 1'])
x = 1:length(P1.trials);
for i = 1:12
    switch i
        case 1
            muscle_m = 1;
        case 4
            muscle_m = 2;
        case 7
            muscle_m = 3;
        case 10
            muscle_m = 4;
    end

    nexttile(i);

    if mod(i, 3) == 1
        P = 0.9;
    elseif mod(i, 3) == 2
        P = 0.8;
    else
        P = 0.6;
    end

    [epoch_numbers, removed_epochs_total] = ...
        removed_epochs_count(data, P, P1, P3, P6, ...
        err_P1_flx_sorted, err_P3_flx_sorted, err_P6_flx_sorted);


    h1 = bar(x, epoch_numbers(P1.trials), 'EdgeColor', 'none'); 
    xlabel('Trials')
    ylabel('Epoch Count')
    h1.FaceColor = 'flat';
    h1.CData = colors(P1.trials, :);
    
    hold on;
    
    h2 = bar(x, removed_epochs_total(muscle_m, P1.trials), 'EdgeColor', 'none');
    h2.FaceColor = 'flat';
    h2.CData = repmat(0.2*[1 1 1], numel(removed_epochs_total(1,P1.trials)), 1);
    
    xlim([1, numel(epoch_numbers(P1.trials))])
    ylim([0 max(epoch_numbers(P1.trials))+1])
    title([Names{muscle_m}, ', ', num2str(100*(1-P)), '% removal '])
end


% P3
figure();
h = tiledlayout(4,3);
title(h, [epoch_base, ', Pressure 3'])
x = 1:length(P3.trials);
for i = 1:12
    switch i
        case 1
            muscle_m = 1;
        case 4
            muscle_m = 2;
        case 7
            muscle_m = 3;
        case 10
            muscle_m = 4;
    end

    nexttile(i);

    if mod(i, 3) == 1
        P = 0.9;
    elseif mod(i, 3) == 2
        P = 0.8;
    else
        P = 0.6;
    end

    [epoch_numbers, removed_epochs_total] = ...
        removed_epochs_count(data, P, P1, P3, P6, ...
        err_P1_flx_sorted, err_P3_flx_sorted, err_P6_flx_sorted);


    h1 = bar(x, epoch_numbers(P3.trials), 'EdgeColor', 'none'); 
    xlabel('Trials')
    ylabel('Epoch Count')
    h1.FaceColor = 'flat';
    h1.CData = colors(P3.trials, :);
    
    hold on;
    
    h2 = bar(x, removed_epochs_total(muscle_m, P3.trials), 'EdgeColor', 'none');
    h2.FaceColor = 'flat';
    h2.CData = repmat(0.2*[1 1 1], numel(removed_epochs_total(1,P3.trials)), 1);
    
    xlim([1, numel(epoch_numbers(P3.trials))])
    ylim([0 max(epoch_numbers(P3.trials))+1])
    title([Names{muscle_m}, ', ', num2str(100*(1-P)), '% removal '])
end


% P6
figure();
h = tiledlayout(4,3);
title(h, [epoch_base, ', Pressure 6'])
x = 1:length(P6.trials);
for i = 1:12
    switch i
        case 1
            muscle_m = 1;
        case 4
            muscle_m = 2;
        case 7
            muscle_m = 3;
        case 10
            muscle_m = 4;
    end

    nexttile(i);

    if mod(i, 3) == 1
        P = 0.9;
    elseif mod(i, 3) == 2
        P = 0.8;
    else
        P = 0.6;
    end

    [epoch_numbers, removed_epochs_total] = ...
        removed_epochs_count(data, P, P1, P3, P6, ...
        err_P1_flx_sorted, err_P3_flx_sorted, err_P6_flx_sorted);


    h1 = bar(x, epoch_numbers(P6.trials), 'EdgeColor', 'none'); 
    xlabel('Trials')
    ylabel('Epoch Count')
    h1.FaceColor = 'flat';
    h1.CData = colors(P6.trials, :);
    
    hold on;
    
    h2 = bar(x, removed_epochs_total(muscle_m, P6.trials), 'EdgeColor', 'none');
    h2.FaceColor = 'flat';
    h2.CData = repmat(0.2*[1 1 1], numel(removed_epochs_total(1,P6.trials)), 1);
    
    xlim([1, numel(epoch_numbers(P6.trials))])
    ylim([0 max(epoch_numbers(P6.trials))+1])
    title([Names{muscle_m}, ', ', num2str(100*(1-P)), '% removal '])
end


%% Deeper look into the EMG across the course of experiment
% Calculate RMS
for i = 1:length(EMG_Preprocessed)

    if ~isempty(EMG_Preprocessed{1, i}.without_outlier_removal.value)
        EMG_Preprocessed{1, i}.without_outlier_removal.RMS = [];
        l = size(EMG_Preprocessed{1, i}.without_outlier_removal.value, 3);
        for j = 1:l
            EMG = EMG_Preprocessed{1, i}.without_outlier_removal.value(1:4, :, j);
            EMG_Preprocessed{1, i}.without_outlier_removal.RMS = cat(2, ...
                EMG_Preprocessed{1, i}.without_outlier_removal.RMS, ...
                sqrt(mean((1e3*EMG).^2, 2)));
        end
    end

end

% define a new vector to store trials id for each pressure
EMG_RMS_pressure = zeros(1, length(EMG_Preprocessed)); 
EMG_RMS_pressure(1, P1.trials) = 1;
EMG_RMS_pressure(1, P3.trials) = 3;
EMG_RMS_pressure(1, P6.trials) = 6;
% remove the bad_trials_EEG_based indexes
EMG_RMS_pressure(bad_trials_EEG_based) = [];


%% Plot the results
figure();
h = tiledlayout(2,2);
title(h, epoch_base)
ylim_max = [];
for i = 1:4
    nexttile(i); hold on;
    % P1
    x = find(EMG_RMS_pressure == 1);
    y = [];
    y_std = [];
    for j = 1:length(P1.trials)
        y = cat(2, y, mean(EMG_Preprocessed{1, P1.trials(j)}.without_outlier_removal.RMS(i, :)));
        y_std = cat(2, y_std, std(EMG_Preprocessed{1, P1.trials(j)}.without_outlier_removal.RMS(i, :)));
    end
    h1 = errorbar(x, y, y_std, 'o', 'Color', All_colours.dark_blue, 'MarkerSize', 5,...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', All_colours.dark_blue);

    % P3
    x = find(EMG_RMS_pressure == 3);
    y = [];
    y_std = [];
    for j = 1:length(P3.trials)
        y = cat(2, y, mean(EMG_Preprocessed{1, P3.trials(j)}.without_outlier_removal.RMS(i, :)));
        y_std = cat(2, y_std, std(EMG_Preprocessed{1, P3.trials(j)}.without_outlier_removal.RMS(i, :)));
    end
    h2 = errorbar(x, y, y_std, 'o', 'Color', All_colours.dark_orange, 'MarkerSize', 5,...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', All_colours.dark_orange);

    % P6
    x = find(EMG_RMS_pressure == 6);
    y = [];
    y_std = [];
    for j = 1:length(P6.trials)
        y = cat(2, y, mean(EMG_Preprocessed{1, P6.trials(j)}.without_outlier_removal.RMS(i, :)));
        y_std = cat(2, y_std, std(EMG_Preprocessed{1, P6.trials(j)}.without_outlier_removal.RMS(i, :)));
    end
    h3 = errorbar(x, y, y_std, 'o', 'Color', All_colours.dark_green, 'MarkerSize', 5,...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', All_colours.dark_green);
    
    ylimit = get(gca, 'ylim');
    ylim_max = cat(2, ylim_max, ylimit(2));

    xlim([1, length(EMG_RMS_pressure)]);
    xlabel('Trials')
    ylabel('EMG RMS')
    title(Names{i})
end

% % Set the maximum y-axis limits for all tiles
% allAxes = findall(h, 'Type', 'axes');
% set(allAxes, 'YLim', [0 max(ylim_max)]); 

% Create a legend for the entire tiled layout
lgd = legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom

































%% Compare pressure effect on EMG signals

figure();
alpha = 0.3;
tiledlayout(2,2)
sgtitle('Flexion Epochs (Mean $\pm$ STD)', 'Interpreter', 'latex');
for i = 1:size(EMG_P1, 1)-6 
    nexttile; hold on;

    upper_curve = 1e3*(mean(EMG_P1_selected(i, :, :), 3) + std(EMG_P1_selected(i, :, :), 0, 3));
    lower_curve = 1e3*(mean(EMG_P1_selected(i, :, :), 3) - std(EMG_P1_selected(i, :, :), 0, 3));

    h1 = fill([X_P1, fliplr(X_P1)], [upper_curve, fliplr(lower_curve)], All_colours.light_blue, ...
        'FaceColor', All_colours.light_blue, 'EdgeColor', 'none', 'FaceAlpha', alpha);
    h2 = plot(X_P1, 1e3*mean(EMG_P1_selected(i, :, :), 3), ...
        'Color', All_colours.dark_blue, 'LineWidth', 2);


    upper_curve = 1e3*(mean(EMG_P3_selected(i, :, :), 3) + std(EMG_P3_selected(i, :, :), 0, 3));
    lower_curve = 1e3*(mean(EMG_P3_selected(i, :, :), 3) - std(EMG_P3_selected(i, :, :), 0, 3));

    h3 = fill([X_P3, fliplr(X_P3)], [upper_curve, fliplr(lower_curve)], All_colours.light_blue, ...
        'FaceColor', All_colours.light_orange, 'EdgeColor', 'none', 'FaceAlpha', alpha);
    h4 = plot(X_P3, 1e3*mean(EMG_P3_selected(i, :, :), 3), ...
        'Color', All_colours.dark_orange, 'LineWidth', 2);


    upper_curve = 1e3*(mean(EMG_P6_selected(i, :, :), 3) + std(EMG_P6_selected(i, :, :), 0, 3));
    lower_curve = 1e3*(mean(EMG_P6_selected(i, :, :), 3) - std(EMG_P6_selected(i, :, :), 0, 3));

    h5 = fill([X_P6, fliplr(X_P6)], [upper_curve, fliplr(lower_curve)], All_colours.light_blue, ...
        'FaceColor', All_colours.light_green, 'EdgeColor', 'none', 'FaceAlpha', alpha);
    h6 = plot(X_P6, 1e3*mean(EMG_P6_selected(i, :, :), 3), ...
        'Color', All_colours.dark_green, 'LineWidth', 2);


    hold off

    title(Names{i})
    xlabel('Cycle [%]')
    ylabel('mV')
end

% Create a legend for the entire tiled layout
lgd = legend([h2, h4, h6], ...
    {'Pressure 1 bar', 'Pressure 3 bar', 'Pressure 6 bar'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom

































%% EMG preprocessing
% Author: Sonja Hanek
% This script applies the preprocessing to the EMG data and structures the
% EMG data.
% Preprocessing consists of filtering, rectification, removal of outliers,
% amplitude normalization and event extraction.
% The resulting data is structured as follows:
% EMG_data --- muscle names 
%          |-- time stamps for the raw data
%          |-- raw data of the EMG sensors
%          |-- preprocessed data of each trial  
%                   |--- trial_nr
%                   |-- pressure level of the trial
%                   |-- scoring of the trial
%                   |-- flexions            
%                   |        |-- preprocessed data of individual flexions
%                   |            within the trial
%                   |-- extensions            
%                   |        |-- preprocessed data of individual extensions
%                   |            within the trial
%                   |-- epoch packs
%                   |        |-- preprocessed data of individual epoch
%                   |        |    packs, i.e. flexion-extension events
%                   |        |-- index of the start of the flexion within
%                   |             the epoch pack
%                   |-- epoch timing
%                           |-- time stamps of beginnings of flexions of
%                           |    the trial in the timestamps of EMG data
%                           |-- index of beginnings of flexions of
%                           |    the trial within EMG data
%                           |-- time stamps of endings of flexions of
%                           |    the trial in the timestamps of EMG data
%                           |-- index of endings of flexions of
%                           |    the trial within EMG data
%                           |-- time stamps of beginnings of extensions of
%                           |    the trial in the timestamps of EMG data
%                           |-- index of beginnings of extensions of
%                           |    the trial within EMG data
%                           |-- time stamps of endings of extensions of
%                           |    the trial in the timestamps of EMG data
%                           |-- index of endings of extensions of
%                               the trial within EMG data
%                                               
% Authors: 
% Sonja Hanek (abc@gmail.com) 
% Morteza Khosrotabar (mkhosrotabar@gmail.com)

%% Load the data
% Before running this section you need to modify the paths in the 
% xdf_load_matlab.m function based on your directory
% you need to run trial_press_score first, then find_epochs.m  !

% All signals from all sessions concatenated (it takes time!)
subject_id = 7;

output = runs_concatenated(subject_id, rawdata_path);
All_EMG = output.All_EMG;
All_EMG_time = output.All_EMG_time;

%%
% % for faster loading, save data
%  load('All_EMG_raw.mat')
%  load('All_EMG_time.mat')

load([path, 'subj_', num2str(subject_id),'_epoch_timestamps.mat'])
load([path,'subj_',num2str(subject_id),'_trial_pressure_score.mat'])



%% filter data

% measurement at 2 kHz 
fs=2e3;
%first muscle is just zeros, since we used the sensors 2-11
EMG_filtered=zeros(10,length(All_EMG));
for muscle=2:11

    % get data
    Raw=All_EMG(muscle,:);



    %bandpass filter
    
    %BANDPASS FILTER 
    %A bandpass filter is used to filter the EMG.
    
    fn=fs/2;    %Hz Nyquist Frequency is 1/2 Sampling Frequency
    low_freq=20;    %Hz
    high_freq=450;  %Hz
    order=4;
    [b,a]=butter(order,([low_freq high_freq]/fn)); % determine filter coefficients
    EMG_r_filt=filtfilt(b,a,Raw); %filtfilt provides zero-lag

    % rectification
    rect=abs(EMG_r_filt);

    %LOWPASS
    fc = 4; % Cut-off frequency (Hz)
    order = 2; % Filter order
    [b,a] = butter(order,fc/fn);
    EMG_r_filt=filtfilt(b,a,rect); %filtfilt provides zero-lag

    EMG_filtered(muscle-1,:)=EMG_r_filt;
end


% plot(All_EMG_time, EMG_filtered(1, :))
% zoom on
% hold on
% plot(All_EMG_time, EMG_filtered(3, :))

%% separate events
% find and save timestamps and indices of the beginning and end of each
% flexion and extension

Emg_trial_epochs=cell(1,length(epochs));

for trial_nr=1:length(epochs)
    trial=epochs{1,trial_nr};
    nr_flexions=length(trial.flexion_start);
    emg_flexion_start_time=zeros(1,nr_flexions);
    emg_flexion_start_index=zeros(1,nr_flexions);

    emg_flexion_end_time=zeros(1,nr_flexions);
    emg_flexion_end_index=zeros(1,nr_flexions);

    nr_extensions=length(trial.extension_start);
    emg_extension_start_time=zeros(1,nr_extensions);
    emg_extension_start_index=zeros(1,nr_extensions);

    emg_extension_end_time=zeros(1,nr_extensions);
    emg_extension_end_index=zeros(1,nr_extensions);

    for flexion=1:nr_flexions
        [~,start_index]=min(abs(All_EMG_time-trial.flexion_start(flexion)));
        [~, end_index]=min(abs(All_EMG_time-trial.flexion_end(flexion)));
        emg_flexion_start_time(flexion)=All_EMG_time(start_index);
        emg_flexion_start_index(flexion)=start_index;
        emg_flexion_end_time(flexion)=All_EMG_time(end_index);
        emg_flexion_end_index(flexion)=end_index;

    end
    for extension=1:nr_extensions
        [~,start_index]=min(abs(All_EMG_time-trial.extension_start(extension)));
        [~, end_index]=min(abs(All_EMG_time-trial.extension_end(extension)));
        emg_extension_start_time(extension)=All_EMG_time(start_index);
        emg_extension_start_index(extension)=start_index;
        emg_extension_end_time(extension)=All_EMG_time(end_index);
        emg_extension_end_index(extension)=end_index;

    end
    Emg_trial_epochs{trial_nr}=struct('emg_flexion_start_time', emg_flexion_start_time, ...
        'emg_flexion_start_index',emg_flexion_start_index, ...
        'emg_flexion_end_time',emg_flexion_end_time, ...
        'emg_flexion_end_index',emg_flexion_end_index, ...
        'emg_extension_start_time',emg_extension_start_time, ...
        'emg_extension_start_index',emg_extension_start_index, ...
        'emg_extension_end_time',emg_extension_end_time, ...
        'emg_extension_end_index',emg_extension_end_index);
end



%% Remove outliers
% For each muscle and each pressure level, epoch packs are identified.
% Epoch packs consist of one flexion and one extension, i.e. they last from
% the beginning of one flexion to the beginning of the next flexion within 
% the trial.
% The epoch packs are time normalized, i.e. interpolated to the sample
% count of the longest epoch pack within the muscle-pressure-level
% category. The error to the median of the epoch packs in the category is 
% calculated:
% The sum of [absolute distances of each sample to the median plus one, to
% the power of 6]
% The 300 epoch packs with the lowest error are kept, the EMG of all other
% epochs is set to NaN.
% This ensures that the evaluated data is of epochs where the muscle
% activity expresses a somewhat typical behaviour and outliers are removed.
% This does not mean that we don't exclude epochs that are actually no
% outliers.
% 

trial_epoch_packs=cell(1,length(Emg_trial_epochs));
for muscle=1:10
    for pressure=[1,3,6]
        epoch_pack_indices=[];
        for trial_nr=1:length(Emg_trial_epochs)
           if trial_pressure_score{1,trial_nr}.pressure==pressure
                flexion_starts=Emg_trial_epochs{1,trial_nr}.emg_flexion_start_index;
                epoch_pack=cell(1,length(flexion_starts)-1);
                for flexion = 1:length(flexion_starts)-1
                    epoch_pack_indices=[epoch_pack_indices;
                        flexion_starts(flexion),flexion_starts(flexion+1),trial_nr];
                end
           end
        end
        epoch_pack_indices(:,4)=epoch_pack_indices(:,2)-epoch_pack_indices(:,1);
        time_norm_length=max(epoch_pack_indices(:,4));
        time_norm_epochs=zeros(length(epoch_pack_indices),time_norm_length);
        for pack_nr=1:length(epoch_pack_indices)
            pack=epoch_pack_indices(pack_nr,:);
            epoch_pack=EMG_filtered(muscle,pack(1):pack(2));
            epoch_time_norm=interp1(epoch_pack,(1:time_norm_length)/(time_norm_length/length(epoch_pack)));
            time_norm_epochs(pack_nr,:)=epoch_time_norm;
        end
        median_pack=median(time_norm_epochs,1,"omitnan");
        errors=abs(time_norm_epochs-median_pack);
        errors=(errors+1).^6;
        errors=sum(errors,2, "omitnan");
        outlier_detection_matrix=[epoch_pack_indices, errors,time_norm_epochs];
        outlier_detection_matrix=sortrows(outlier_detection_matrix,5);
%        % uncomment to examine plots of errors
%         figure;plot(outlier_detection_matrix(:,5))
%         title(join(['muscle: ',num2str(muscle),' pressure: ', num2str(pressure)]))
        for outs=301:length(errors) % keep only the 300 epoch packs (around half) that are most consistent
            EMG_filtered(muscle,outlier_detection_matrix(outs,1):outlier_detection_matrix(outs,2))=NaN;
        end
    end
end

%% amplitude normalization
% For each muscle, the maximum value of the muscle activity in the first 
% trial with a pressure level of 1 is taken as reference for amplitude 
% normalization

trial_press=[trial_pressure_score{:}];
trial_press=[trial_press(:).pressure];
ref_trial_nr=find(trial_press==1,1);
ref_trial=Emg_trial_epochs{1,ref_trial_nr};
trial_start_index=min(ref_trial.emg_flexion_start_index(1),ref_trial.emg_extension_start_index(1));
trial_end_index=max(ref_trial.emg_flexion_end_index(end), ref_trial.emg_extension_end_index(end));
ref_val=max(EMG_filtered(:,trial_start_index:trial_end_index),[],2,"omitnan");
for muscle=1:10
    EMG_filtered(muscle,:)=EMG_filtered(muscle,:)/ref_val(muscle);
end

%% structure data

EMG_filtered_epochs=cell(1,length(epochs));
EMG_filtered_epoch_packs=cell(1,length(epoch_pack_indices));
for trial_nr=1:length(epochs)

    flexion=cell(1,length(Emg_trial_epochs{1,trial_nr}.emg_flexion_start_index));
    for flex_nr=1:length(Emg_trial_epochs{1,trial_nr}.emg_flexion_start_index)
        flexion{1,flex_nr}=EMG_filtered(:, ...
            Emg_trial_epochs{1,trial_nr}.emg_flexion_start_index(flex_nr): ...
            Emg_trial_epochs{1,trial_nr}.emg_flexion_end_index(flex_nr));
    end

    extension=cell(1,length(Emg_trial_epochs{1,trial_nr}.emg_extension_start_index));
    
    for ext_nr=1:length(Emg_trial_epochs{1,trial_nr}.emg_extension_start_index)
        extension{1,ext_nr}=EMG_filtered(:, ...
            Emg_trial_epochs{1,trial_nr}.emg_extension_start_index(ext_nr): ...
            Emg_trial_epochs{1,trial_nr}.emg_extension_end_index(ext_nr));
    end

    flexion_starts=Emg_trial_epochs{1,trial_nr}.emg_flexion_start_index;
    extension_starts=Emg_trial_epochs{1,trial_nr}.emg_extension_start_index;
    epoch_pack=cell(1,length(flexion_starts)-1);
    for pack_nr = 1:length(flexion_starts)-1
        ext_start=extension_starts(flexion_starts(pack_nr)< extension_starts & ...
            extension_starts<flexion_starts(pack_nr+1))-flexion_starts(pack_nr);
        epoch_pack{1,pack_nr}=struct(...
            'filtered_data', EMG_filtered(:,...
            flexion_starts(pack_nr):flexion_starts(pack_nr+1)),...
            'relative_extension_start_index', ext_start);
    end

        
    EMG_filtered_epochs{trial_nr}=struct( ...
        'trial_nr', trial_nr, ...
        'pressure', epochs{1,trial_nr}.pressure, ...
        'score',epochs{1,trial_nr}.score, ...
        'flexions', {flexion}, ...    
        'extensions',{extension}, ...
        'epoch_packs', {epoch_pack},...
        'epoch_timing', Emg_trial_epochs{1,trial_nr});
end

names={'Vastus_med_R', 'Rectus_femoris_R', 'Gastrocnemius_R', 'Biceps_femoris_R', ...
    'Vastus_med_L', 'Rectus_femoris_L', 'Gastrocnemius_L', 'Biceps_femoris_L', ...
    'Trapezius_R', 'Trapezius_L'};

EMG_data=struct( ...
    'Muscle_names', {names},...
    'Timestamps', All_EMG_time, ...
    'Raw_data', All_EMG, ...
    'Epochs_preprocessed', {EMG_filtered_epochs});

%% save data
filename=sprintf([path,'subj_',num2str(subject_id),'_emg.mat']);
disp(filename)
save(filename, 'EMG_data')