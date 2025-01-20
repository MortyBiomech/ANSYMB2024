clc
clear

%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
ROIs_data_path = [data_path, '8_Classification\ROIs_features\'];


%% Load Region of Interest files
cd(ROIs_data_path)
epoch_type = 'FlextoFlex';
load(['ROIs_0_', epoch_type, '.mat'])


%% Create tables for Linear Mixed Effect Model analysis

regions_names = fieldnames(ROIs);
frequency_bands = {'delta'; 'theta'; 'alpha'; 'beta'; 'gamma'};

numFreqBands = size(frequency_bands, 1);
numRegions = size(regions_names, 1);

RMS_Freq_Region = cell(size(frequency_bands, 1), size(regions_names, 1));
for i = 1:numRegions

    region_data = ROIs.(regions_names{i});

    for j = 1:numFreqBands  

        Subject_ID = [];
        IC_ID = [];
        Condition_ID = [];
        RMS_value = [];
        for k = 1:size(region_data, 1)

            % P1
            RMS_P1 = region_data{k, 3}.P1(:, j);
            subject_id_count = repmat(region_data{k, 1}, numel(RMS_P1), 1);
            ic_id_count = repmat(region_data{k, 2}, numel(RMS_P1), 1);
            condition_id_count = 1*ones(numel(RMS_P1), 1);

            RMS_value = cat(1, RMS_value, RMS_P1);
            Subject_ID = cat(1, Subject_ID, subject_id_count);
            IC_ID = cat(1, IC_ID, ic_id_count);
            Condition_ID = cat(1, Condition_ID, condition_id_count);


            % P3
            RMS_P3 = region_data{k, 3}.P3(:, j);
            subject_id_count = repmat(region_data{k, 1}, numel(RMS_P3), 1);
            ic_id_count = repmat(region_data{k, 2}, numel(RMS_P3), 1);
            condition_id_count = 3*ones(numel(RMS_P3), 1);

            RMS_value = cat(1, RMS_value, RMS_P3);
            Subject_ID = cat(1, Subject_ID, subject_id_count);
            IC_ID = cat(1, IC_ID, ic_id_count);
            Condition_ID = cat(1, Condition_ID, condition_id_count);


            % P6
            RMS_P6 = region_data{k, 3}.P6(:, j);
            subject_id_count = repmat(region_data{k, 1}, numel(RMS_P6), 1);
            ic_id_count = repmat(region_data{k, 2}, numel(RMS_P6), 1);
            condition_id_count = 6*ones(numel(RMS_P6), 1);

            RMS_value = cat(1, RMS_value, RMS_P6);
            Subject_ID = cat(1, Subject_ID, subject_id_count);
            IC_ID = cat(1, IC_ID, ic_id_count);
            Condition_ID = cat(1, Condition_ID, condition_id_count);

        end
        T = table(Subject_ID, IC_ID, Condition_ID, RMS_value);
        T.Subject_ID = categorical(T.Subject_ID);
        T.IC_ID = categorical(T.IC_ID);
        T.Condition_ID = categorical(T.Condition_ID);

        RMS_Freq_Region{j, i} = T;

    end

end



%% Remove this later after checking the data and finding the problem. 
% For now to plot the heatmaps I will remove subject 15 at Frontal region
% on theta band. Then I will fit the LMM models.
Frontal_theta = RMS_Freq_Region{2, 3};

% Specify the Subject_ID to be removed
subjectToRemove = categorical(15);

% Create a logical index for rows where Subject_ID is not equal to the specified value
rowsToKeep = Frontal_theta.Subject_ID ~= subjectToRemove;

% Apply the logical index to retain only the rows to keep
Frontal_theta = Frontal_theta(rowsToKeep, :);

% Substitute the region data table
RMS_Freq_Region{2, 3} = Frontal_theta;


%% Fit Linear Mixed Effect Model
results = cell(size(frequency_bands, 1), size(regions_names, 1)); 

for freq = 1:numFreqBands   

    for region = 1:numRegions

        % Extract the data table for this cell
        data = RMS_Freq_Region{freq, region};

        % Define the linear mixed-effects model
        % RMS_value ~ Condition_ID + (1|Subject_ID) + (1|Subject_ID:IC_ID)
        formula = 'RMS_value ~ Condition_ID + (1|Subject_ID) + (1|Subject_ID:IC_ID)';

        % First fit: Condition 1 as reference
        data.Condition_ID = reordercats(data.Condition_ID, {'1', '3', '6'});
        lme1 = fitlme(data, formula);

        % Second fit: Condition 2 as reference
        data.Condition_ID = reordercats(data.Condition_ID, {'3', '1', '6'});
        lme2 = fitlme(data, formula);

        % Extract results (e.g., p-values for fixed effects)
        coef1 = anova(lme1); 
        coef2 = anova(lme2); 

        % Save results
        results{freq, region} = struct('Model1', lme1, 'coef1', coef1, ...
            'Model2', lme2, 'coef2', coef2);

    end

end

disp('Model fitting completed for all cells.');



%% Collect Adjusted P-Values and Estimates

% For Condition 1 vs. 2
heatmap_diff_1_2 = nan(numFreqBands, numRegions); % Effect sizes
significance_1_2 = false(numFreqBands, numRegions); % Significance markers

% For Condition 1 vs. 3
heatmap_diff_1_3 = nan(numFreqBands, numRegions);
significance_1_3 = false(numFreqBands, numRegions);

% For Condition 2 vs. 3
heatmap_diff_2_3 = nan(numFreqBands, numRegions);
significance_2_3 = false(numFreqBands, numRegions);

% Adjusted p-value arrays
adjusted_pvalues_1_2 = []; % For Condition 1 vs. 2
adjusted_pvalues_1_3 = []; % For Condition 1 vs. 3
adjusted_pvalues_2_3 = []; % For Condition 2 vs. 3

% Loop through each frequency band and region
for freq = 1:numFreqBands
    for region = 1:numRegions
        % Access Model1 (Condition 1 as reference)
        model1 = results{freq, region}.Model1; % Replace with your data structure
        [fixedEffectsEstimates1, names1, stats1] = fixedEffects(model1);

        % Extract Condition 1 vs. 2
        cond2_index = strcmp(names1.Name, 'Condition_ID_3'); % Assuming Condition_ID_3 represents Condition 2
        if any(cond2_index)
            heatmap_diff_1_2(freq, region) = fixedEffectsEstimates1(cond2_index);
            adjusted_pvalues_1_2 = [adjusted_pvalues_1_2; stats1.pValue(cond2_index)];
        end

        % Extract Condition 1 vs. 3
        cond3_index = strcmp(names1.Name, 'Condition_ID_6'); % Assuming Condition_ID_6 represents Condition 3
        if any(cond3_index)
            heatmap_diff_1_3(freq, region) = fixedEffectsEstimates1(cond3_index);
            adjusted_pvalues_1_3 = [adjusted_pvalues_1_3; stats1.pValue(cond3_index)];
        end

        % Access Model2 (Condition 2 as reference)
        model2 = results{freq, region}.Model2; % Replace with your data structure
        [fixedEffectsEstimates2, names2, stats2] = fixedEffects(model2);

        % Extract Condition 2 vs. 3
        cond2_3_index = strcmp(names2.Name, 'Condition_ID_6'); % Assuming Condition_ID_6 represents Condition 3
        if any(cond2_3_index)
            heatmap_diff_2_3(freq, region) = fixedEffectsEstimates2(cond2_3_index);
            adjusted_pvalues_2_3 = [adjusted_pvalues_2_3; stats2.pValue(cond2_3_index)];
        end
    end
end

% Apply FDR correction
adjusted_pvalues_1_2 = mafdr(adjusted_pvalues_1_2, 'BHFDR', true);
adjusted_pvalues_1_3 = mafdr(adjusted_pvalues_1_3, 'BHFDR', true);
adjusted_pvalues_2_3 = mafdr(adjusted_pvalues_2_3, 'BHFDR', true);

% Update significance markers
index = 1;
for freq = 1:numFreqBands
    for region = 1:numRegions
        if index <= length(adjusted_pvalues_1_2)
            significance_1_2(freq, region) = adjusted_pvalues_1_2(index) < 0.05;
        end
        if index <= length(adjusted_pvalues_1_3)
            significance_1_3(freq, region) = adjusted_pvalues_1_3(index) < 0.05;
        end
        if index <= length(adjusted_pvalues_2_3)
            significance_2_3(freq, region) = adjusted_pvalues_2_3(index) < 0.05;
        end
        index = index + 1;
    end
end


%% Plot heatmaps

% Define extreme colors
synch_color = [214, 40, 40]/255; % Replace with your desired RGB value for the maximum
desynch_color = [58, 134, 255]/255;  % Replace with your desired RGB value for the minimum

% Define the number of colors for the colormap
num_colors = 256;

% Create a gradient for the negative side (red to white)
neg_colors = [linspace(desynch_color(1), 1, num_colors/2)', ...
              linspace(desynch_color(2), 1, num_colors/2)', ...
              linspace(desynch_color(3), 1, num_colors/2)'];

% Create a gradient for the positive side (white to blue)
pos_colors = [linspace(1, synch_color(1), num_colors/2)', ...
              linspace(1, synch_color(2), num_colors/2)', ...
              linspace(1, synch_color(3), num_colors/2)'];

% Combine the gradients to create the full colormap
custom_cmap = [neg_colors; pos_colors];

% Calculate the global range for all heatmaps
global_min = min([heatmap_diff_1_2(:); heatmap_diff_1_3(:); heatmap_diff_2_3(:)]);
global_max = max([heatmap_diff_1_2(:); heatmap_diff_1_3(:); heatmap_diff_2_3(:)]);

% Set symmetric limits around zero for consistent colormap
global_limit = max(abs([global_min, global_max]));



% Labels for brain regions and frequency bands
brain_regions_complete = {'Left PreMot SuppMot', 'Left Paracentral Lobule', ...
    'Left Dorsal ACC', 'Left VisMotor', 'Left PrimVisual', ...
    'Right PreMot SuppMot', 'Right VisMotor', 'Right PrimVisual'}; 
brain_regions = {'LPS', 'LPL', 'F', 'LP', 'LO', 'RPS', 'RO', 'RP'}; 
freq_bands = {'\delta', '\theta', '\alpha', '\beta', '\gamma'};


% Plot for Condition 1 vs. 2
figure;
imagesc(flipud(heatmap_diff_1_2)); % Use imagesc for axes-based plotting
colormap(custom_cmap); % Set colormap
clim([-global_limit, global_limit]); % Set symmetric color limits based on data

% Add colorbar and label
cb = colorbar;
cb.Label.String = 'Effect Size'; % Label for color bar
cb.Label.FontSize = 14; % Set font size for the label
cb.Ticks = linspace(-global_limit, global_limit, 5); % Generate ticks
cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), cb.Ticks, 'UniformOutput', false); % Format ticks to 2 decimals

title('Pressure 3 bar - Pressure 1 bar', 'FontSize', 14, 'FontWeight', 'normal');
xlabel('Brain Regions', 'FontSize', 14);
ylabel('Frequency Bands', 'FontSize', 14);

% Reverse the y-axis labels and align with flipped data
yticks(1:length(freq_bands));
yticklabels(flip(freq_bands)); % Flip the labels to match the flipped data

% Customize axis ticks and labels
xticks(1:length(brain_regions));
xticklabels(brain_regions);

% Increase font size of axis labels
set(gca, 'FontSize', 14);

% Overlay significance markers
[row, col] = find(flipud(significance_1_2)); % Find significant cells
for i = 1:length(row)
    text(col(i), row(i), '*', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'Color', 'black', ...
        'FontWeight', 'bold', 'FontSize', 20);
end

% Make the figure more horizontal
set(gcf, 'Position', [100, 100, 700, 400]); % Adjust figure size


% Plot for Condition 1 vs. 3
figure;
imagesc(flipud(heatmap_diff_1_3)); % Use imagesc for axes-based plotting
colormap(custom_cmap); % Set colormap
clim([-global_limit, global_limit]); % Set symmetric color limits based on data

% Add colorbar and label
cb = colorbar;
cb.Label.String = 'Effect Size'; % Label for color bar
cb.Label.FontSize = 14; % Set font size for the label
cb.Ticks = linspace(-global_limit, global_limit, 5); % Generate ticks
cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), cb.Ticks, 'UniformOutput', false); % Format ticks to 2 decimals

title('Pressure 6 bar - Pressure 1 bar', 'FontSize', 14, 'FontWeight', 'normal');
xlabel('Brain Regions');
ylabel('Frequency Bands');

% Reverse the y-axis labels and align with flipped data
yticks(1:length(freq_bands));
yticklabels(flip(freq_bands)); % Flip the labels to match the flipped data

% Customize axis ticks and labels
xticks(1:length(brain_regions));
xticklabels(brain_regions);

% Increase font size of axis labels
set(gca, 'FontSize', 14);

% Overlay significance markers
[row, col] = find(flipud(significance_1_3)); % Find significant cells
for i = 1:length(row)
    text(col(i), row(i), '*', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'Color', 'black', ...
        'FontWeight', 'bold', 'FontSize', 20);
end

% Make the figure more horizontal
set(gcf, 'Position', [100, 100, 700, 400]); % Adjust figure size



% Plot for Condition 2 vs. 3
figure;
imagesc(flipud(heatmap_diff_2_3)); % Use imagesc for axes-based plotting
colormap(custom_cmap); % Set colormap
clim([-global_limit, global_limit]); % Set symmetric color limits based on data

% Add colorbar and label
cb = colorbar;
cb.Label.String = 'Effect Size'; % Label for color bar
cb.Label.FontSize = 14; % Set font size for the label
cb.Ticks = linspace(-global_limit, global_limit, 5); % Generate ticks
cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), cb.Ticks, 'UniformOutput', false); % Format ticks to 2 decimals

title('Pressure 6 bar - Pressure 3 bar', 'FontSize', 14, 'FontWeight', 'normal');
xlabel('Brain Regions')
ylabel('Frequency Bands');

% Reverse the y-axis labels and align with flipped data
yticks(1:length(freq_bands));
yticklabels(flip(freq_bands)); % Flip the labels to match the flipped data

% Customize axis ticks and labels
xticks(1:length(brain_regions));
xticklabels(brain_regions);

% Increase font size of axis labels
set(gca, 'FontSize', 14);

% Overlay significance markers
[row, col] = find(flipud(significance_2_3)); % Find significant cells
for i = 1:length(row)
    text(col(i), row(i), '*', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'Color', 'black', ...
        'FontWeight', 'bold', 'FontSize', 20);
end

% Make the figure more horizontal
set(gcf, 'Position', [100, 100, 700, 400]); % Adjust figure size




