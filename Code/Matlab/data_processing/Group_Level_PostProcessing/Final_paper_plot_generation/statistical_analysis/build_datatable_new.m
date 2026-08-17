function tbl = build_datatable_new(cluster_data, fIdx, tIdx)
% Assemble trial-wise data table for a given (freq,time) sample across all ICs
% 
% INPUT:
%   cluster_data: Cell array where each row represents one IC
%                 Column 1: Subject ID
%                 Column 2: IC ID
%                 Column 3: Struct with fields P1, P3, P6
%   fIdx, tIdx:   Frequency and time indices
%
% OUTPUT:
%   tbl: Table with columns Power, Pressure, SubjectID

    power_vals = [];
    press_vals = [];
    subj_vals = [];
    ic_vals = [];
    
    nIC = size(cluster_data, 1);
    condition_names = {'P1', 'P3', 'P6'};
    condition_values = [1, 3, 6];  % Map P1->1, P3->3, P6->6
    
    for ic = 1:nIC
        subject_id = cluster_data{ic, 1};  % Column 1: Subject ID
        ic_id = cluster_data{ic, 2};       % Column 2: Subject IC id
        tf_struct = cluster_data{ic, 3};   % Column 3: TF data struct
        
        % Loop through each condition (P1, P3, P6)
        for cond_idx = 1:length(condition_names)
            cond_name = condition_names{cond_idx};
            cond_value = condition_values(cond_idx);
            
            % Check if this condition exists for this IC
            if isfield(tf_struct, cond_name) && ~isempty(tf_struct.(cond_name))
                trials = tf_struct.(cond_name);  % Cell array of trial matrices
                nTrials = length(trials);
                
                % Extract power values for this (freq, time) point across all trials
                trial_powers = zeros(nTrials, 1);
                for trial = 1:nTrials
                    trial_matrix = trials{trial};  % 183 x 1001 matrix
                    trial_powers(trial) = trial_matrix(fIdx, tIdx);
                end
                
                % Append to main arrays
                power_vals = [power_vals; trial_powers]; %#ok<AGROW>
                press_vals = [press_vals; repmat(cond_value, nTrials, 1)]; %#ok<AGROW>
                subj_vals = [subj_vals; repmat(subject_id, nTrials, 1)]; %#ok<AGROW>
                ic_vals = [ic_vals; repmat(ic_id, nTrials, 1)]; %#ok<AGROW>
            end
        end
    end
    
    % Create the table
    tbl = table(power_vals, press_vals, subj_vals, ic_vals, ...
                'VariableNames', {'Power', 'Pressure', 'SubjectID', 'IC_ID'});
end