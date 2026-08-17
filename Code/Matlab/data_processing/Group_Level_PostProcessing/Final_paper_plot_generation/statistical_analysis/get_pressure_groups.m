function [P1, P3, P6, isExperiment] = get_pressure_groups(subject, trial_info)
    % Get experimental trials
    if subject < 10
        isExperiment = 1:length(trial_info);
    else
        isExperiment = find(cellfun(@(x) strcmp(x.General.Description, 'Experiment'), trial_info));
    end
    
    % Extract pressure levels more efficiently
    pressures = cellfun(@(x) x.General.Pressure, trial_info);
    
    P1 = intersect(find(pressures == 1), isExperiment);
    P3 = intersect(find(pressures == 3), isExperiment);
    P6 = intersect(find(pressures == 6), isExperiment);
end