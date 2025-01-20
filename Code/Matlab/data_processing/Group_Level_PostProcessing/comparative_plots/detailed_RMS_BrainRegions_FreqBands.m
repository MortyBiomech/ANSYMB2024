freq_indx = 5;
region_indx = 3;

dataTable = RMS_Freq_Region{freq_indx, region_indx};

% Labels for brain regions and frequency bands
brain_regions_complete = {'Left PreMot SuppMot', 'Left Paracentral Lobule', ...
    'Left Dorsal ACC', 'Left VisMotor', 'Left PrimVisual', ...
    'Right PreMot SuppMot', 'Right VisMotor', 'Right PrimVisual'}; 
brain_regions = {'LPS', 'LPL', 'F', 'LP', 'LO', 'RPS', 'RO', 'RP'}; 
freq_bands = {'\delta', '\theta', '\alpha', '\beta', '\gamma'};

% Extract the data from the table
Condition_ID = dataTable.Condition_ID;
RMS_value = dataTable.RMS_value;
Subject_ID = dataTable.Subject_ID;
IC_ID = dataTable.IC_ID;

% Define Subject_IDs to exclude
excludeSubjects = [15]; % Replace with the Subject_IDs you want to exclude

% Get unique Subject_IDs
uniqueSubjects = unique(Subject_ID);
catnames = categories(uniqueSubjects);
uniqueSubjects = str2double(catnames(uniqueSubjects));

% Shuffle colors for better distinction
rng(0); % Set seed for reproducibility
colors = jet(length(uniqueSubjects));
colors = colors(randperm(size(colors, 1)), :); % Shuffle the rows of the colormap

% Initialize the figure
figure;
hold on;

% Define marker styles
markerStyles = {'o', 's', 'd', '^', 'v', '<', '>', 'p', 'h', '*'}; % Extend as needed

% Offset factor for x-axis to separate Subject/IC combinations
offsetFactor = 0.05; % Adjust this value as needed

% Convert Condition_ID to numeric if it's categorical
if iscategorical(Condition_ID)
    Condition_ID = double(Condition_ID);
end

multiplier = 0;
% Loop through each unique Subject_ID
for s = 1:length(uniqueSubjects)
    realSubjectID = uniqueSubjects(s); % Use the real Subject_ID
    if ismember(realSubjectID, excludeSubjects)
        continue; % Skip the excluded subjects
    end
    
    subjectMask = Subject_ID == categorical(realSubjectID);
    subjectColor = colors(mod(s-1, length(uniqueSubjects)) + 1, :); % Cycle through shuffled colors
    
    % Get the unique IC_IDs for the current Subject_ID
    subjectICs = unique(IC_ID(subjectMask));
    
    % Loop through each unique IC_ID for the current Subject_ID
    for i = 1:length(subjectICs)
        icMask = IC_ID == subjectICs(i);
        combinedMask = subjectMask & icMask;
        
        % Add an offset to the x-values for this Subject/IC combination
        % xOffset = (s - 1) * offsetFactor + (i - 1) * (offsetFactor / length(subjectICs));
        adjustedCondition = Condition_ID(combinedMask) + multiplier*offsetFactor;
        multiplier = multiplier +1;
        
        % Generate a marker for the current IC_ID
        marker = markerStyles{mod(i-1, length(markerStyles)) + 1};
        
        % Plot the data
        catnames2 = categories(subjectICs);
        catnames2 = str2double(catnames2);
        scatter(adjustedCondition, RMS_value(combinedMask), ...
            40, 'MarkerEdgeColor', subjectColor, 'Marker', marker, ...
            'DisplayName', sprintf('Subject %d - IC %d', realSubjectID, catnames2(subjectICs(i))));
    end
end

% Customize the plot
% xlabel('Condition');
ylabel('RMS Value');
legend('show', 'Location', 'bestoutside'); % Display legend
title(strjoin([brain_regions_complete(region_indx), ' - ', freq_bands(freq_indx)], ''))
set(gca, "FontSize", 12)

% Set x-axis ticks to the main three conditions
xticks([1, 2, 3] + ceil((multiplier-1)/2)*offsetFactor); % Replace these values with the actual Condition_ID values if different
xticklabels({'P1', 'P3', 'P6'}); % Customize labels if needed
xlim([1 - 0.25, 3 + (multiplier-1)*offsetFactor + 0.25])

hold off;
