function [RMS_Freq_Region, heatmap_diff_1_2, heatmap_diff_1_3, heatmap_diff_2_3, ...
    significance_1_2, significance_1_3, significance_2_3] = ...
    EEG_features_heatmap_tables_per_pressure(EEG_ROIs_path, filename, subject)

    %% Load Region of Interest files
    data = load(fullfile(EEG_ROIs_path, filename));
    name = fieldnames(data);
    ROIs = data.(name{1});
    
    
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



    %% Select the rows of each subject and remove all other rows
    for freq = 1:numFreqBands   
    
        for region = 1:numRegions

            temp = RMS_Freq_Region{freq, region};
            rowsToKeep = temp.Subject_ID == categorical(subject);
            temp = temp(rowsToKeep, :);
            RMS_Freq_Region{freq, region} = temp;
        
        end

    end
    
   
    
    %% Fit Linear Mixed Effect Model
    results = cell(size(frequency_bands, 1), size(regions_names, 1)); 
    
    for freq = 1:numFreqBands   
    
        for region = 1:numRegions
            
            if isempty(RMS_Freq_Region{freq, region})
                continue
            end

            % Extract the data table for this cell
            data = RMS_Freq_Region{freq, region};
    
            % Define the linear mixed-effects model
            % RMS_value ~ Condition_ID + (1|Subject_ID) + (1|Subject_ID:IC_ID)
            if size(unique(data.IC_ID), 1) ~= 1
                formula = 'RMS_value ~ Condition_ID + (1|IC_ID)';
            elseif size(unique(data.IC_ID), 1) == 1
                formula = 'RMS_value ~ Condition_ID';
            end
    
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
    
    disp('LMM fitting completed for all cells.');
    
    
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

            if isempty(results{freq, region})
                continue
            end

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
    


    %% Update significance markers
    index = 1;
    for freq = 1:numFreqBands
        for region = 1:numRegions

            if isempty(results{freq, region})
                continue
            end

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




end