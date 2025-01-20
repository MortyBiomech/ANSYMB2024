clc
clear

%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
ROIs_data_path = [data_path, '8_Classification\ROIs_features\'];

%% Load Region of Interest files
current_path = pwd;
cd(ROIs_data_path)
epoch_type = 'FlextoFlex';
load(['ROIs_0_', epoch_type, '.mat'])
cd(current_path)


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

%% Fit Linear Mixed Effect Model
results = cell(size(frequency_bands, 1), size(regions_names, 1)); 

for freq = 1:numFreqBands   

    for region = 1:numRegions

        % Extract the data table for this cell
        data = RMS_Freq_Region{freq, region};

        % Define the standard LMM formula
        formula = 'RMS_value ~ Condition_ID + (1|Subject_ID) + (1|Subject_ID:IC_ID)';

        % First fit: Condition 1 as reference
        data.Condition_ID = reordercats(data.Condition_ID, {'1', '3', '6'});
        lme1 = fitlme(data, formula);
        standard_residuals1 = residuals(lme1); % Residuals for first fit

        % Second fit: Condition 3 as reference
        data.Condition_ID = reordercats(data.Condition_ID, {'3', '1', '6'});
        lme2 = fitlme(data, formula);
        standard_residuals2 = residuals(lme2); % Residuals for second fit

        % Robust LMM: First fit (Condition 1 as reference)
        data.Condition_ID = reordercats(data.Condition_ID, {'1', '3', '6'});
        data_struct = table2struct(data, 'ToScalar', true);
        save('data_for_R.mat', 'data_struct');
        system('Rscript robust_lmm.R');
        robust_results1 = load('results.mat');

        % Robust LMM: Second fit (Condition 3 as reference)
        data.Condition_ID = reordercats(data.Condition_ID, {'3', '1', '6'});
        data_struct = table2struct(data, 'ToScalar', true);
        save('data_for_R.mat', 'data_struct');
        system('Rscript robust_lmm.R');
        robust_results2 = load('results.mat');

        % Store residuals and models
        results{freq, region} = struct('StandardModel1', lme1, ...
                                       'StandardModel2', lme2, ...
                                       'StandardResiduals1', standard_residuals1, ...
                                       'StandardResiduals2', standard_residuals2, ...
                                       'RobustFixedEffects1', robust_results1.fixedEffects, ...
                                       'RobustRandomEffects1', robust_results1.randomEffects, ...
                                       'RobustResiduals1', robust_results1.residuals, ...
                                       'RobustFixedEffects2', robust_results2.fixedEffects, ...
                                       'RobustRandomEffects2', robust_results2.randomEffects, ...
                                       'RobustResiduals2', robust_results2.residuals);
    end

end

disp('Model fitting completed for all cells.');