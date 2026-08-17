function [subjectTable, ICs_in_regions] = ...
    createInputTableforClassification(ROI, subject, classes)

    regions = fieldnames(ROI); % Get all region names (e.g., Left_PreMot_SuppMot)

    % Step 1: Initialize an empty table with no columns
    columnNames = {};
    data = [];
    ICs_in_regions = zeros(1, length(regions));

    for regionIdx = 1:numel(regions)
        regionName = regions{regionIdx};
        regionData = ROI.(regionName); % Extract data for the current region

        % Find rows corresponding to the current subject
        subjectRows = find(cell2mat(regionData(:, 1)) == subject);
        if isempty(subjectRows)
            continue; % Skip if no ICs for this subject in this region
        end

        % Step 2: Process each IC for the subject in the current region
        for rowIdx = subjectRows'
            
            icID = regionData{rowIdx, 2}; % IC ID (2nd column)
            featureStruct = regionData{rowIdx, 3}; % 3: nonNorm RMS, 4: Norm RMS

            if strcmp(classes, 'P1P3P6')
                % Three classes (P1, P3, P6)
                data_temp = [featureStruct.P1; featureStruct.P3; featureStruct.P6];
            elseif strcmp(classes, 'P1P6')
                % Two Classes (P1, P6)
                data_temp = [featureStruct.P1; featureStruct.P6]; 
            elseif strcmp(classes, 'P1P3')
                % Two Classes (P1, P3)
                data_temp = [featureStruct.P1; featureStruct.P3];
            elseif strcmp(classes, 'P3P6')
                % Two Classes (P3, P6)
                data_temp = [featureStruct.P3; featureStruct.P6];
            elseif strcmp(classes, 'S1S2S3')
                % Three classes (S1, S2, S3)
                data_temp = [featureStruct.S1; featureStruct.S2; featureStruct.S3];
            elseif strcmp(classes, 'S1S3')
                % Three classes (S1, S3)
                data_temp = [featureStruct.S1; featureStruct.S3];
            elseif strcmp(classes, 'S2S3')
                % Three classes (S2, S3)
                data_temp = [featureStruct.S2; featureStruct.S3];
            elseif strcmp(classes, 'S1S2')
                % Three classes (S1, S2)
                data_temp = [featureStruct.S1; featureStruct.S2];
            end

            data = cat(2, data, data_temp);
            
            % for the case we extracted the flexion and extension features
            % together. PSD was calsulated for the whole FlextoFlex epochs.
            % featureNames = strcat(regionName, '_', {'delta', 'theta', 'alpha', 'beta', 'gamma'}, ...
            %     '_IC', string(icID));
            
            % for the case we extracted the flexion and extension features
            % separately and created a 1*10 feature vector instead of 1*5.
            featureNames = strcat(regionName, '_', ...
                {'delta_Flx', 'theta_Flx', 'alpha_Flx', 'beta_Flx', 'gamma_Flx', ...
                 'delta_Ext', 'theta_Ext', 'alpha_Ext', 'beta_Ext', 'gamma_Ext'}, ...
                '_IC', string(icID));

            columnNames = cat(2, columnNames, featureNames);
            ICs_in_regions(regionIdx) = ICs_in_regions(regionIdx) + 1;

        end
    end

    % Z-score normalization
    data_zscore_norm = zscore(data);

    % Step 3: Add the 'Label' column at the end
    columnNames = cat(2, columnNames, 'Label');
    
    
    % Replace numeric classes with string labels and make them categorical
    if strcmp(classes, 'P1P3P6')
        % Three classes (P1, P3, P6)
        classLabels = [repmat({'P1'}, size(featureStruct.P1, 1), 1); ...
                       repmat({'P3'}, size(featureStruct.P3, 1), 1); ...
                       repmat({'P6'}, size(featureStruct.P6, 1), 1)];
    elseif strcmp(classes, 'P1P6')
        % Two classes (P1, P6)
        classLabels = [repmat({'P1'}, size(featureStruct.P1, 1), 1); ...
                       repmat({'P6'}, size(featureStruct.P6, 1), 1)];
    elseif strcmp(classes, 'P3P6')
        % Two classes (P3, P6)
        classLabels = [repmat({'P3'}, size(featureStruct.P3, 1), 1); ...
                       repmat({'P6'}, size(featureStruct.P6, 1), 1)];
    elseif strcmp(classes, 'P1P3')
        % Two classes (P1, P3)
        classLabels = [repmat({'P1'}, size(featureStruct.P1, 1), 1); ...
                       repmat({'P3'}, size(featureStruct.P3, 1), 1)];
    elseif strcmp(classes, 'S1S2S3')
        % Three classes (S1, S2, S3)
        classLabels = [repmat({'S1'}, size(featureStruct.S1, 1), 1); ...
                       repmat({'S2'}, size(featureStruct.S2, 1), 1); ...
                       repmat({'S3'}, size(featureStruct.S3, 1), 1)];
    elseif strcmp(classes, 'S1S3')
        % Three classes (S1, S3)
        classLabels = [repmat({'S1'}, size(featureStruct.S1, 1), 1); ...
                       repmat({'S3'}, size(featureStruct.S3, 1), 1)];
    elseif strcmp(classes, 'S2S3')
        % Three classes (S2, S3)
        classLabels = [repmat({'S2'}, size(featureStruct.S2, 1), 1); ...
                       repmat({'S3'}, size(featureStruct.S3, 1), 1)];
    elseif strcmp(classes, 'S1S2')
        % Three classes (S1, S2)
        classLabels = [repmat({'S1'}, size(featureStruct.S1, 1), 1); ...
                       repmat({'S2'}, size(featureStruct.S2, 1), 1)];
    end 
    
    classLabels = categorical(classLabels); % Convert to categorical
    
    % Create the table for the current subject
    subjectTable = array2table(data_zscore_norm, 'VariableNames', columnNames(1:end-1));
    subjectTable.Label = classLabels; % Add the 'Pressure' column with categorical labels
    

    
    

end