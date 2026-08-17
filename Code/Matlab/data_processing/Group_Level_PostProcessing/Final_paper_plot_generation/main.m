clc
clear


%% Add essential paths
data_path = 'D:\Morteza\MyProjects\ANSYMB2024\data\';
study_path = 'D:\Morteza\MyProjects\ANSYMB2024\data\7_STUDY\multiple_clustering\_Final_Selected_Clusters\';
processing_path = 'D:\Morteza\MyProjects\ANSYMB2024\Code\Matlab\data_processing\';
epoched_data_path = [data_path, '6_Trials_Info_and_Epoched_data\'];

addpath(genpath('D:\Morteza\MyProjects\ANSYMB2024'))
addpath(genpath('D:\Morteza\LSL\xdf-Matlab-master'))


%% Run EEGLAB
path = pwd;
cd('D:\Morteza\Toolboxes\EEGLAB\eeglab2025.1.0')
if ~exist('ALLCOM','var')
    eeglab;
end
cd(path)


%% Load epoched data one time (takes time)
subjects_data = struct();
subjects_list = 5:18;
for subject = subjects_list

    disp(['Subject ', num2str(subject)])

    data = load([epoched_data_path, 'sub-', num2str(subject), ...
        '\Epochs_FlextoFlex_based.mat']);
    name = fieldnames(data);
    data = data.(name{1});

    subjects_data.(['sub', num2str(subject)]) = data;

end

subjects_trialsInfo = struct();
for subject = subjects_list
    
    Trials_Info = load([epoched_data_path, 'sub-', num2str(subject), ...
        '\Trials_Info.mat']);
    name = fieldnames(Trials_Info);
    Trials_Info = Trials_Info.(name{1});

    subjects_trialsInfo.(['sub', num2str(subject)]) = Trials_Info;

end



%% Load STUDY files from multiple clustering solutions

filename_pre = 'main_study_potential_brain_ICs_RV-15_';
% %%%% Left_Dorsal_ACC (20, 57, -10)
% filePath = [study_path, 'Left_Dorsal_ACC (20, 57, -10)'];
% fileName = [filename_pre, 'Left_Dorsal_ACC (20, 57, -10).study'];
% [STUDY, ALLEEG] = pop_loadstudy('filename', fileName, ...
%     'filepath', filePath);
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
% eeglab redraw;

%%%% Right_Prim_Motor (43, -27, 41)
filePath = [study_path, 'Right_Prim_Motor (43, -27, 41)'];
fileName = [filename_pre, 'Right_Prim_Motor (43, -27, 41).study'];
[STUDY, ALLEEG] = pop_loadstudy('filename', fileName, ...
    'filepath', filePath);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
eeglab redraw; 

% %%%% Left_Prim_Motor (-38, -28, 59)
% filePath = [study_path, 'Left_Prim_Motor (-38, -28, 59)'];
% fileName = [filename_pre, 'Left_Prim_Motor (-38, -28, 59).study'];
% [STUDY, ALLEEG] = pop_loadstudy('filename', fileName, ...
%     'filepath', filePath);
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
% eeglab redraw;

% %%%% Right_PreMot_SuppMot (15, 15, 60)
% filePath = [study_path, 'Right_PreMot_SuppMot (15, 15, 60)'];
% fileName = [filename_pre, 'Right_PreMot_SuppMot (15, 15, 60).study'];
% [STUDY, ALLEEG] = pop_loadstudy('filename', fileName, ...
%     'filepath', filePath);
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
% eeglab redraw;  

% %%%% Left_PreMot_SuppMot (-15, 15, 60)
% filePath = [study_path, 'Left_PreMot_SuppMot (-15, 15, 60)'];
% fileName = [filename_pre, 'Left_PreMot_SuppMot (-15, 15, 60).study'];
% [STUDY, ALLEEG] = pop_loadstudy('filename', fileName, ...
%     'filepath', filePath);
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
% eeglab redraw;  

% %%%% Right_Parieto_Occipital (35, -65, 40)
% filePath = [study_path, 'Right_Parieto_Occipital (35, -65, 40)'];
% fileName = [filename_pre, 'Right_Parieto_Occipital (35, -65, 40).study'];
% [STUDY, ALLEEG] = pop_loadstudy('filename', fileName, ...
%     'filepath', filePath);
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
% eeglab redraw;  

% %%%% Left_Parieto_Occipital (-35, -65, 40)
% filePath = [study_path, 'Left_Parieto_Occipital (-35, -65, 40)'];
% fileName = [filename_pre, 'Left_Parieto_Occipital (-35, -65, 40).study'];
% [STUDY, ALLEEG] = pop_loadstudy('filename', fileName, ...
%     'filepath', filePath);
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
% eeglab redraw;  

% %%%% Prime_Visual (0, -85, 10)
% filePath = [study_path, 'Prime_Visual (0, -85, 10)'];
% fileName = [filename_pre, 'Prime_Visual (0, -85, 10).study'];
% [STUDY, ALLEEG] = pop_loadstudy('filename', fileName, ...
%     'filepath', filePath);
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
% eeglab redraw;  

% %%%% Paracentral_Lobule (-2, -10, 60)
% filePath = [study_path, 'Paracentral_Lobule (-2, -10, 60)'];
% fileName = [filename_pre, 'Paracentral_Lobule (-2, -10, 60).study'];
% [STUDY, ALLEEG] = pop_loadstudy('filename', fileName, ...
%     'filepath', filePath);
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
% eeglab redraw;  


% filename_pre = 'main_study_potential_brain_ICs_RV-15_';

all_study_names = {'Right_Prim_Motor (43, -27, 41)', ...
                   'Right_PreMot_SuppMot (15, 15, 60)', ...
                   'Right_Parieto_Occipital (35, -65, 40)', ...
                   'Prime_Visual (0, -85, 10)', ...
                   'Paracentral_Lobule (-2, -10, 60)', ...
                   'Left_Prim_Motor (-38, -28, 59)', ...
                   'Left_PreMot_SuppMot (-15, 15, 60)', ...
                   'Left_Parieto_Occipital (-35, -65, 40)', ...
                   'Left_Dorsal_ACC (20, 57, -10)'};

sets = cell(1, length(all_study_names));
sets_names = {'Right_Prim_Motor', ...
              'Right_PreMot_SuppMot', ...
              'Right_Parieto_Occipital', ...
              'Prime_Visual', ...
              'Paracentral_Lobule', ...
              'Left_Prim_Motor', ...
              'Left_PreMot_SuppMot', ...
              'Left_Parieto_Occipital', ...
              'Left_Dorsal_ACC'};

% for i = 1:length(all_study_names)
% 
%     filePath = [study_path, all_study_names{i}];
%     fileName = [filename_pre, all_study_names{i}, '.study'];
%     [STUDY, ALLEEG] = pop_loadstudy('filename', fileName, ...
%         'filepath', filePath);
%     CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
%     eeglab redraw; 
% 
%     % extract Subjects and Components
%     best_cluster_id = STUDY.etc.bemobil.clustering.cluster_ROI_index;
%     subjects = STUDY.cluster(best_cluster_id).sets + 4;
%     components = STUDY.cluster(best_cluster_id).comps;
%     subjects_components = [subjects', components'];
% 
%     sets{i} = subjects_components;
% 
% end

%%%%%%%%%%%%%%%% sets was saved to the drive on 24.08.2025 %%%%%%%%%%%%%%%%
%%% load on new runs
load('Clusters__subjects_components_24082025.mat');


% %% Create a map to store ICs and their corresponding clusters
% point_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
% 
% % Iterate through each cluster
% for i = 1:length(sets)
%     current_set = sets{i};
%     for j = 1:size(current_set, 1)
%         point = current_set(j, :);
%         % Create a string key for the point
%         key = sprintf('%.0f,%.0f', point(1), point(2));
% 
%         if isKey(point_map, key)
%             % Point already exists, add current set to its list
%             existing_sets = point_map(key);
%             existing_sets{end+1} = sets_names{i}; %#ok<SAGROW>
%             point_map(key) = existing_sets;
%         else
%             % New point, create new entry
%             point_map(key) = {sets_names(i)};
%         end
%     end
% end
% 
% %% Find ICs that belong to more than one cluster
% point_coords = [];
% set_memberships = {};
% point_strings = {};
% keys = point_map.keys;
% 
% fprintf('Points that belong to multiple sets:\n');
% for i = 1:length(keys)
%     sets_for_point = point_map(keys{i});
%     if length(sets_for_point) > 1
%         % Parse the point coordinates back from the key
%         coords = str2double(strsplit(keys{i}, ','));
% 
%         % Store data for table
%         point_coords = cat(1, point_coords, coords'); 
%         set_memberships{end+1} = sets_for_point; %#ok<SAGROW>
%         point_strings{end+1} = sprintf('[%.0f, %.0f]', coords(1), coords(2)); %#ok<SAGROW>
% 
%         fprintf('Point [%.0f, %.0f] belongs to sets: %s\n', ...
%             coords(1), coords(2), strjoin(sets_for_point, ', '));
%     end
% end
% 
% % Create the table
% if ~isempty(point_coords)
%     common_points_table = table(point_coords(:,1), point_coords(:,2), ...
%                                set_memberships', ...
%                                'VariableNames', {'Subject', 'IC', 'Clusters'});
% 
%     % Display the table
%     disp('Table of Common Points:');
%     disp(common_points_table);
% else
%     disp('No common points found.');
% end

%%%%%%%%% common_points_table was saved to the drive on 24.08.2025 %%%%%%%%
% Conclusion after looking at the common ICs:
% (1) 
% Cluster Paracentral_Lobule has the most common ICs with other clusters
% Except common ICs in this cluster there are 5 subjects which the ICs
% don't have good quality. So the Paracentral Lobule cluster is not
% consideref for further analysis. 
% (2)
% Subject 11 IC 29 is assigned to the Left_Premotor Cluster

%%% Load common ICs table
% load('Common_ICs_among_Clusters__24082025.mat');


%% Generate the figure containing all clusters with their topoplots
studyFiles = cell(1, 8); % eight files (all except paracentral lobule)
all_study_names(5) = [];
for i = 1:length(studyFiles)
    studyFiles{i} = [study_path, all_study_names{i}, '\', filename_pre, all_study_names{i}, '.study']; 
end
% 8 colors for 8 clusters
colors = [0, 134, 103; ...
          82, 48, 192; ...
          192, 60, 66; ...
          202, 182, 0; ...
          83, 0, 143; ...
          0, 218, 234; ...
          0, 76, 145; ...
          223, 128, 0]./255;

% open plotMultiStudyDipolesWithMRI.m

%% load one study to have all EEG files
%%%% Right_Prim_Motor (43, -27, 41)
filePath = [study_path, 'Right_Prim_Motor (43, -27, 41)'];
fileName = [filename_pre, 'Right_Prim_Motor (43, -27, 41).study'];
[STUDY, ALLEEG] = pop_loadstudy('filename', fileName, ...
    'filepath', filePath);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
eeglab redraw; 



%% Optional
% topoplots of all components in the cluster
STUDY = std_topoplot(STUDY,ALLEEG,'clusters',best_cluster_id, 'design', 1, 'plotsubjects', 'on' );



%% TF_calculation_per_IC   %%%%%% has been saved on disc (18 GB)
% load("Clusters__subjects_components_24082025.mat");
% sets(5) = []; % removing the Paracentral Lobule cluster
% load("Clusters__names_24082025.mat")
% sets_names(5) = [];
% 
% All_Clusters_TF_data = struct();
% 
% for c = 2:numel(sets)
% 
%     subjects_components = sets{c};
% 
%     cluster_TF_data = cell(size(subjects_components,1), 4);
%     % Extranct TF data for each IC in the cluster
%     for i = 1:size(subjects_components, 1)
% 
%         subject = subjects_components(i, 1); 
%         ic = subjects_components(i, 2);
%         disp(['Subject ', num2str(subject), ' IC ', num2str(ic), ' ', sets_names{c}])
%         sub_id = sprintf('sub%d', subject);
%         trial_info = subjects_trialsInfo.(sub_id);
%         subject_data = subjects_data.(sub_id);
% 
%         % for baseline calculations:
%         EEGLAB_source = ALLEEG(subject-4).icaact(ic, :);
%         EEGLAB_events = ALLEEG(subject-4).event;
% 
% 
%         [TF_struct, ~, freqs_baseline] = ...
%             main_TF_calculation_per_IC(subject, ic, trial_info, subject_data, ...
%                                   EEGLAB_source, EEGLAB_events);
% 
%         cluster_TF_data{i, 1} = subject;
%         cluster_TF_data{i, 2} = ic;
%         cluster_TF_data{i, 3} = TF_struct;
%         cluster_TF_data{i, 4} = freqs_baseline;
% 
%     end
%     cluster_TF_data = cell2table(cluster_TF_data, ...
%         "VariableNames", ["Subject", "IC", "TF_Data", "Frequency"]);
% 
%     All_Clusters_TF_data.(sets_names{c}) = cluster_TF_data;
% 
% end

% All_Clusters_TF_data was saved to the disk on 27.08.2025 ~ 18 GB 
% load("Clusters__Time_Frequency_data_27082025.mat") % takes some time!




%% finding median and time-warping process

All_Clusters_TF_data_warped = struct();
flexion_extension_median_lengths = struct();
for c = 1:numel(sets)

    cluster_TF_data = All_Clusters_TF_data.(sets_names{c});

    %% finding the medain flexion and extension in the cluster for time-warping
    all_flx_lens = [];
    all_ext_lens = [];
    for i = 1:size(cluster_TF_data, 1)
    
        all_flx_lens = cat(2, all_flx_lens, ...
            [cell2mat(cluster_TF_data.TF_Data(i, 1).P1(2, :)), ...
             cell2mat(cluster_TF_data.TF_Data(i, 1).P3(2, :)), ...
             cell2mat(cluster_TF_data.TF_Data(i, 1).P6(2, :))]); 
    
        all_ext_lens = cat(2, all_ext_lens, ...
            [cell2mat(cluster_TF_data.TF_Data(i, 1).P1(3, :)), ...
             cell2mat(cluster_TF_data.TF_Data(i, 1).P3(3, :)), ...
             cell2mat(cluster_TF_data.TF_Data(i, 1).P6(3, :))]); 
    
    end
    
    median_all_flx = round(median(all_flx_lens));
    median_all_ext = round(median(all_ext_lens));
    
    flexion_extension_median_lengths.(sets_names{c}) = [median_all_flx, median_all_ext];
    
    %% Time warping of TF data in the cluster
    cluster_TF_data_warped = table2cell(cluster_TF_data);
    for i = 1:size(cluster_TF_data_warped, 1)
    
        disp(['subject ', num2str(cluster_TF_data_warped{i, 1}), ...
            ' IC ', num2str(cluster_TF_data_warped{i, 2})])
        for p = ["P1", "P3", "P6"]
    
            before_interp_time_points = cellfun(@(x) size(x, 2), cluster_TF_data_warped{i, 3}.(p)(1, :), 'UniformOutput', false);
            before_interp_time_vectors = cellfun(@(x) linspace(0, 100, x), before_interp_time_points, 'UniformOutput', false);
        
            % flexion part
            flexion_lens = cluster_TF_data_warped{i, 3}.(p)(2, :); % flexion lengths 
            after_interp_time_vector_flx = cellfun(@(x, y) linspace(x(1), x(y), median_all_flx), ...
                before_interp_time_vectors, flexion_lens, 'UniformOutput', false);
            % flexion_TF_warped = cellfun(@(x, x_flx, y, z) interp1(x(1:x_flx)', y(:, 1:x_flx)', z, "linear"), ...
            %     before_interp_time_vectors, flexion_lens, cluster_TF_data_warped{i, 3}.(p)(1, :), after_interp_time_vector_flx, ...
            %     'UniformOutput', false);
            flexion_TF_warped = cellfun(@(x, x_flx, y, z) interp1(x(1:x_flx)', pagetranspose(y(:, 1:x_flx, :)), z, "linear"), ...
                before_interp_time_vectors, flexion_lens, cluster_TF_data_warped{i, 3}.(p)(1, :), after_interp_time_vector_flx, ...
                'UniformOutput', false);
            flexion_TF_warped = cellfun(@(x) pagetranspose(x), flexion_TF_warped, 'UniformOutput', false);
        
            % extension part
            after_interp_time_vector_ext = cellfun(@(x, y) linspace(x(y+1), x(end), median_all_ext), ...
                before_interp_time_vectors, flexion_lens, 'UniformOutput', false);
            % extension_TF_warped = cellfun(@(x, x_flx, y, z) interp1(x(x_flx+1:end)', y(:, x_flx+1:end)', z, "linear"), ...
            %     before_interp_time_vectors, flexion_lens, cluster_TF_data_warped{i, 3}.(p)(1, :), after_interp_time_vector_ext, ...
            %     'UniformOutput', false);
            extension_TF_warped = cellfun(@(x, x_flx, y, z) interp1(x(x_flx+1:end)', pagetranspose(y(:, x_flx+1:end, :)), z, "linear"), ...
                before_interp_time_vectors, flexion_lens, cluster_TF_data_warped{i, 3}.(p)(1, :), after_interp_time_vector_ext, ...
                'UniformOutput', false);
            extension_TF_warped = cellfun(@(x) pagetranspose(x), extension_TF_warped, 'UniformOutput', false);
        
            % combine warped flexion and extension parts
            TF_warped = cellfun(@(x, y) [x, y], flexion_TF_warped, extension_TF_warped, 'UniformOutput', false);
        
        
            % fill the cell with the warped TF data
            cluster_TF_data_warped{i, 3}.(p) = TF_warped;
    
        end
        
    end

    All_Clusters_TF_data_warped.(sets_names{c}) = cluster_TF_data_warped;

end



%% Main Statistical Analysis

% to free up the memory!
clear All_Clusters_TF_data

% TFCE Statistical Analysis
% ANALYSIS PARAMETERS
p_value = 0.01;                 % Significance level
num_permutations = 1000;        % Number of permutations
tfce_H = 2;                     % TFCE height exponent
tfce_E = 0.5;                   % TFCE extent exponent
dh = 0.1;                       % TFCE integration step

% SEGMENTATION PARAMETERS (for the sake of calculation time!)
time_segment_percent = 0.02;    % Use 2% time segments (50 segments total)
freq_group_size = 2;            % Group every 2 frequencies together

All_Clusters_TFCE_data = struct();
for c = 6:numel(sets)
    
    cluster_TF_data_warped = All_Clusters_TF_data_warped.(sets_names{c});

    % Get original dimensions
    [orig_nFreq, orig_nTime] = size(cluster_TF_data_warped{1,3}.P1{1});
    nIC = size(cluster_TF_data_warped, 1);
    
    % Calculate segmentation parameters
    [seg_nFreq, seg_nTime, freq_groups, time_segments] = ...
        calculate_segmentation_params(orig_nFreq, orig_nTime, freq_group_size, time_segment_percent);
    
    % Prepare segmented data structure with trials and averages
    segmented_data = prepare_segmented_trial_data(cluster_TF_data_warped, freq_groups, time_segments);
    
    % Compute observed F-map on segmented data
    observed_Fmap = compute_segmented_fmap(segmented_data, seg_nFreq, seg_nTime);
    
    % Apply TFCE to observed F-map
    observed_TFCE = tfce2d(observed_Fmap, tfce_H, tfce_E, dh);
    
    % Permutation testing with trial-level shuffling
    max_TFCE_null = zeros(num_permutations, 1);
    
    % Use smaller parallel pool for memory efficiency
    current_pool = gcp('nocreate');
    if isempty(current_pool)
        parpool('threads', min(4, feature('numcores')));
    end
    

    % Create progress object (set useWaitbar=false for console prints)
    pp = ParforProgress(num_permutations, 'Progress', 5, true); % update every 10 ticks

    % DataQueue + callback that calls the separate .m file method
    dq = parallel.pool.DataQueue;
    afterEach(dq, @(~) pp.tick());  % runs on client
    
    
    parfor perm = 1:num_permutations
        % Create trial-level permuted segmented averages
        perm_seg_averages = create_trial_permuted_segmented_data(segmented_data);
        
        % Compute F-map for permuted data
        perm_Fmap = compute_segmented_fmap(perm_seg_averages, seg_nFreq, seg_nTime);
        
        % Apply TFCE
        perm_TFCE = tfce2d(perm_Fmap, tfce_H, tfce_E, dh);
        
        % Store maximum TFCE value
        max_TFCE_null(perm) = max(perm_TFCE(:));


        send(dq, 1);  % ping progress
    end
    
    % Determine significance threshold
    tfce_threshold = prctile(max_TFCE_null, 100*(1-p_value));
    
    % Generate final significance mask
    significance_mask = observed_TFCE > tfce_threshold;
    
    
    % Keep only clusters (8-connected) that contain ≥8 points
    conn = conndef(2,'maximal');        % 8-connected neighbourhood in 2-D
    CC   = bwconncomp(significance_mask, conn);
    
    final_mask = false(size(significance_mask));
    cluster_sizes = cellfun(@numel, CC.PixelIdxList);
    
    for cl = 1:CC.NumObjects
        if cluster_sizes(cl) >= 8          % cluster size threshold
            final_mask(CC.PixelIdxList{cl}) = true;
        end
    end



    % save data in a structure
    P1_TF = []; P3_TF = []; P6_TF = [];
    for i = 1:size(cluster_TF_data_warped, 1)
        TF_P1 = cat(3, cluster_TF_data_warped{i, 3}.P1{:});
        TF_P3 = cat(3, cluster_TF_data_warped{i, 3}.P3{:});
        TF_P6 = cat(3, cluster_TF_data_warped{i, 3}.P6{:});
    
        P1_TF = cat(3, P1_TF, mean(TF_P1, 3));
        P3_TF = cat(3, P3_TF, mean(TF_P3, 3));
        P6_TF = cat(3, P6_TF, mean(TF_P6, 3));
    end

    % Normalize to the extreme
    P1_TF_norm = zeros(size(P1_TF));
    P3_TF_norm = zeros(size(P3_TF));
    P6_TF_norm = zeros(size(P6_TF));
    for i = 1:size(cluster_TF_data_warped, 1)
    
        M = max(max( abs(P1_TF(:, :, i)) ));
    
        P1_TF_norm(:, :, i) = P1_TF(:, :, i) / M;
        P3_TF_norm(:, :, i) = P3_TF(:, :, i) / M;
        P6_TF_norm(:, :, i) = P6_TF(:, :, i) / M;
    
    end

    P1_TF = mean(P1_TF_norm, 3);
    P3_TF = mean(P3_TF_norm, 3);
    P6_TF = mean(P6_TF_norm, 3);

    freqs = cluster_TF_data_warped{1, 4};
    freqs_stat = cellfun(@(x) mean(freqs(x)), freq_groups);

    innerstruct = struct();
    innerstruct.P1_TF = P1_TF;
    innerstruct.P3_TF = P3_TF;
    innerstruct.P6_TF = P6_TF;
    innerstruct.F_map = observed_Fmap;
    innerstruct.TFCE_map = observed_TFCE;
    innerstruct.max_TFCE_null = max_TFCE_null;
    innerstruct.p_value = p_value;
    innerstruct.final_mask = final_mask;
    innerstruct.flexion_median_length = median_all_flx;
    innerstruct.extension_median_length = median_all_ext;
    innerstruct.freqs_TF = freqs;
    innerstruct.freqs_stat = freqs_stat;

    
    All_Clusters_TFCE_data.(sets_names{c}) = innerstruct;
    
    
end



%% plot the TF ERS/ERD data

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


% TF_P1 = []; TF_P3 = []; TF_P6 = [];
% for i = 1:size(cluster_TF_data_warped, 1)
% 
%     matrix_3D = cat(3, cluster_TF_data_warped{i, 3}.P1{:});
%     TF_P1 = cat(3, TF_P1, matrix_3D);
% 
%     matrix_3D = cat(3, cluster_TF_data_warped{i, 3}.P3{:});
%     TF_P3 = cat(3, TF_P3, matrix_3D);
% 
%     matrix_3D = cat(3, cluster_TF_data_warped{i, 3}.P6{:});
%     TF_P6 = cat(3, TF_P6, matrix_3D);
% 
% end
% 
% mean_TF_P1 = mean(TF_P1, 3);
% mean_TF_P3 = mean(TF_P3, 3);
% mean_TF_P6 = mean(TF_P6, 3);
% 
% 
% T = linspace(0, 100, median_all_flx + median_all_ext);
% 
% 
% 
% %%
% figure()
% tiledlayout(1,3, "TileSpacing", "compact", "Padding", "compact")
% 
% % Calculate the global range for all heatmaps
% global_min = min(min([mean_TF_P1, mean_TF_P3, mean_TF_P6]));
% global_max = max(max([mean_TF_P1, mean_TF_P3, mean_TF_P6]));
% % Set symmetric limits around zero for consistent colormap
% global_limit = max(abs([global_min, global_max]));
% 
% nexttile
% 
% pcolor(T, freqs, mean_TF_P1); hold on
% shading interp; 
% xline(T(median_all_flx))
% colormap(custom_cmap)
% clim([-global_limit, global_limit])
% set(gca, 'YScale', 'log');  
% set(gca, 'FontSize', 14)
% yticks([2 4 8 14 30 50])
% title('P1', 'FontSize', 16, 'FontWeight', 'normal')
% xlabel('Cycle(%)')
% ylabel('Frequency (Hz)')
% 
% nexttile
% 
% pcolor(T, freqs, mean_TF_P3); hold on
% shading interp; 
% xline(T(median_all_flx))
% colormap(custom_cmap)
% clim([-global_limit, global_limit])
% set(gca, 'YScale', 'log');  
% set(gca, 'FontSize', 14)
% yticks([2 4 8 14 30 50])
% title('P3', 'FontSize', 16, 'FontWeight', 'normal')
% xlabel('Cycle(%)')
% 
% nexttile
% 
% pcolor(T, freqs, mean_TF_P6); hold on
% shading interp; 
% xline(T(median_all_flx))
% colormap(custom_cmap)
% clim([-global_limit, global_limit])
% set(gca, 'YScale', 'log');  
% set(gca, 'FontSize', 14)
% yticks([2 4 8 14 30 50])
% title('P6', 'FontSize', 16, 'FontWeight', 'normal')
% xlabel('Cycle(%)')
% 
% cb = colorbar('eastoutside');  
% cb.Label.String = 'z-scored ERS/ERD';  
% cb.Label.FontSize = 14;
% 
% set(gcf, "Position", [100 315 1200 400])


%% plot TF ERS/ERD data per subject
figure()
H = height(cluster_TF_data_warped);
maintile = tiledlayout(3, H, "Padding", "compact", "TileSpacing", "compact");

T = linspace(0, 100, median_all_flx + median_all_ext);
for i = 1:size(cluster_TF_data_warped, 1)
    
    TF_P1 = cat(3, cluster_TF_data_warped{i, 3}.P1{:});
    TF_P3 = cat(3, cluster_TF_data_warped{i, 3}.P3{:});
    TF_P6 = cat(3, cluster_TF_data_warped{i, 3}.P6{:});

    % Calculate the global range for all heatmaps
    global_min = min(min([mean(TF_P1, 3), mean(TF_P3, 3), mean(TF_P6, 3)]));
    global_max = max(max([mean(TF_P1, 3), mean(TF_P3, 3), mean(TF_P6, 3)]));
    % Set symmetric limits around zero for consistent colormap
    global_limit = max(abs([global_min, global_max]));
    
    climit = [-global_limit, global_limit];

    % subject = cluster_TF_data_warped{i, 1};
    % if subject < 10 
    %     climit = [0, global_max];
    % else
    %     climit = [-global_limit, global_limit];
    % end


    freq = cluster_TF_data_warped{i, 4};


    nexttile(i)

    pcolor(T, freq, mean(TF_P1, 3)); hold on;
    shading flat; 
    xline(T(median_all_flx))
    colormap('turbo')
    % colormap(custom_cmap)
    clim(climit)
    % clim(climit)
    set(gca, 'YScale', 'log');  
    set(gca, 'FontSize', 12)
    yticks([2 4 8 14 30 50])
    title(['P1 - ', 'S', num2str(cluster_TF_data_warped{i, 1}), ...
        ' IC', num2str(cluster_TF_data_warped{i, 2})], ...
        'FontSize', 12, 'FontWeight', 'normal')
    

    nexttile(i + H)
    
    pcolor(T, freq, mean(TF_P3, 3)); hold on
    shading interp; 
    xline(T(median_all_flx))
    colormap('turbo')
    % colormap(custom_cmap)
    clim(climit)
    % clim(climit)
    set(gca, 'YScale', 'log');  
    set(gca, 'FontSize', 12)
    yticks([2 4 8 14 30 50])
    title(['P3 - ', 'S', num2str(cluster_TF_data_warped{i, 1}), ...
        ' IC', num2str(cluster_TF_data_warped{i, 2})], ...
        'FontSize', 12, 'FontWeight', 'normal')
    

    nexttile(i + 2*H)
    
    pcolor(T, freq, mean(TF_P6, 3)); hold on
    shading interp; 
    xline(T(median_all_flx))
    colormap('turbo')
    % colormap(custom_cmap)
    clim(climit)
    % clim(climit)
    set(gca, 'YScale', 'log');  
    set(gca, 'FontSize', 12)
    yticks([2 4 8 14 30 50])
    xlabel('Cycle (%)')
    title(['P6 - ', 'S', num2str(cluster_TF_data_warped{i, 1}), ...
        ' IC', num2str(cluster_TF_data_warped{i, 2})], ...
        'FontSize', 12, 'FontWeight', 'normal')
    
    cb = colorbar('southoutside');  
    ylabel(cb, '$10\log_{10}(\mu V^2/\mathrm{Hz})$', 'Interpreter', 'latex');
    % cb.Label.String = '$10\log_{10}(\mu V^2)$';  % Add label
    
    % cb.Label.String = 'z-scored ERS/ERD';  % Add label
    % cb.Label.String = 'Power (\muV^2)';  % Add label
    
    % if subject < 10
    %     cb.Label.String = 'Power (\muV^2)';  % Add label
    % else
    %     cb.Label.String = '10$log$(\muV^2/Hz)';  % Add label
    % end

    cb.Label.FontSize = 12;

    if i == 1
        nexttile(i)
        ylabel('Frequency (Hz)')

        nexttile(i+H)
        ylabel('Frequency (Hz)')
        
        nexttile(i+2*H)
        ylabel('Frequency (Hz)')
    end
    

    
end

% title(maintile, 'Left Dorsal ACC', 'FontSize', 16)
title(maintile, 'Left Paracentral Lobule', 'FontSize', 16)
set(gcf, "Position", [1750 -80 2450 1250])


%% plot the mean TF of the cluster (average on subject-ic)
c = 2;
cluster_TF_data_warped = All_Clusters_TF_data_warped.(sets_names{c});

% making the colormap for the stat plot
middle_c = [152, 30, 30]/255;
lower_c = [255, 239, 232]/255;
higher_c = [47, 0, 0]/255;
num_colors = 256;
low_colors = [linspace(lower_c(1), middle_c(1), num_colors/2)', ...
              linspace(lower_c(2), middle_c(2), num_colors/2)', ...
              linspace(lower_c(3), middle_c(3), num_colors/2)'];
up_colors = [linspace(middle_c(1), higher_c(1), num_colors/2)', ...
              linspace(middle_c(2), higher_c(2), num_colors/2)', ...
              linspace(middle_c(3), higher_c(3), num_colors/2)'];
% Combine the gradients to create the full colormap
cmap_stat = [low_colors; up_colors];


% making the colormap for the TF plots
% Define RGB values (0–255)
rgb_values = [...
     66,  85, 255;   % RISD Blue
    102, 173, 255;   % Ruddy Blue (#66ADFF)
    112, 210, 255;   % Pale Azure (#70D2FF)
     122, 246, 255;   % Electric Blue (#7AF6FF)
    152, 251, 192;   % Aquamarine (#98FBC0)
    181, 255, 128;   % Screamin' Green (#B5FF80)
    208, 255, 125;   % Mindaro (#D0FF7D)
    235, 255, 122;   % Mindaro (#EBFF7A)
    255, 249, 133;   % Icterine (#FFF985)
    255, 234, 128;   % Jasmine (#FFEA80)
    255, 200, 133;   % Sunset (#FFC885)
    255, 146, 110;   % Atomic Tangerine (#FF926E)
    255,  92,  87];  % Bittersweet (#FF5C57)
% rgb_values = [...
%       0, 120, 129;
%      12, 157, 168;
%      40, 203, 215;
%      55, 223, 236;
%      70, 229, 241;   % Electric Blue (#46E5F1)
%      90, 235, 246;   % Electric Blue (#5AEBF6)
%     110, 239, 249;   % Electric Blue (#6EEFF9)
%     125, 242, 251;   % Electric Blue (#7DF2FB)
%     140, 245, 253;   % Electric Blue (#8CF5FD)
%     166, 248, 254;   % Celeste (#A6F8FE)
%     192, 251, 255;   % Celeste (#C0FBFF)
%     240, 254, 255;   % Azure (Web) (#F0FEFF)
%     255, 255, 255;   % White (#FFFFFF)
%     255, 240, 240;   % Lavender Blush (#FFF0F0)
%     255, 217, 216;   % Misty Rose (#FFD9D8)
%     255, 190, 188;   % Melon (#FFBEBC)
%     255, 163, 160;   % Melon (#FFA3A0)
%     255,  86,  80;   % Bittersweet (#FF5650)
%     250,  47,  40;   % Vermilion (#FA2F28)
%     244,   7,   0;   % Off Red (RGB) (#F40700)
%     232,   6,   0;   % Off Red (RGB) (#E80600)
%     215,   5,   0;   % Engineering Orange (#D70500)
%     207,   4,   0;   % Engineering Orange (#CF0400)
%     198,   3,   0;   % Engineering Orange (#C60300)
%     162,   0,   0];  % Penn Red (#A20000)
% Normalize to 0–1 range
cmap = rgb_values ./ 255;
% % Number of colors in final colormap (e.g. 256 for smooth)
% n = 256;
% % Interpolate between the 5 colors
% cmap_cont = interp1(...
%     linspace(1, n, size(cmap,1)), cmap, ...
%     1:n, 'linear');
cmap_cont = makeColormap(rgb_values, 256, 'sigmoid', 7);



T = linspace(0, 100, sum(flexion_extension_median_lengths.(sets_names{c})));
T_stat = linspace(0, 100, numel(time_segments));

freqs = cluster_TF_data_warped{1, 4};
freqs_stat = cellfun(@(x) mean(freqs(x)), freq_groups);

TF_P1_cluster = [];
TF_P3_cluster = [];
TF_P6_cluster = [];
% selected_ICs = [1, 2, 3, 4, 6, 8, 12, 15, 16, 17, 18, 20]; % Left Prieto-Occipital
% selected_ICs = [4, 5, 6, 11, 12, 13, 17, 18]; % Left Paracentral Lobule
% cluster_TF_data_warped = cluster_TF_data_warped(selected_ICs, :);
for i = 1:size(cluster_TF_data_warped, 1)
    
    TF_P1 = cat(3, cluster_TF_data_warped{i, 3}.P1{:});
    TF_P3 = cat(3, cluster_TF_data_warped{i, 3}.P3{:});
    TF_P6 = cat(3, cluster_TF_data_warped{i, 3}.P6{:});

    TF_P1_cluster = cat(3, TF_P1_cluster, mean(TF_P1, 3));
    TF_P3_cluster = cat(3, TF_P3_cluster, mean(TF_P3, 3));
    TF_P6_cluster = cat(3, TF_P6_cluster, mean(TF_P6, 3));

end

% % Calculate the global range for all heatmaps
% global_min = min(min([mean(TF_P1_cluster, 3), mean(TF_P3_cluster, 3), mean(TF_P6_cluster, 3)]));
% global_max = max(max([mean(TF_P1_cluster, 3), mean(TF_P3_cluster, 3), mean(TF_P6_cluster, 3)]));
% % Set symmetric limits around zero for consistent colormap
% global_limit = max(abs([global_min, global_max]));


% Normalize to the extreme
TF_P1_cluster_norm = [];
TF_P3_cluster_norm = [];
TF_P6_cluster_norm = [];
for i = 1:size(cluster_TF_data_warped, 1)

    M = max(max( abs(TF_P1_cluster(:, :, i)) ));

    TF_P1_cluster_norm(:, :, i) = TF_P1_cluster(:, :, i) / M;
    TF_P3_cluster_norm(:, :, i) = TF_P3_cluster(:, :, i) / M;
    TF_P6_cluster_norm(:, :, i) = TF_P6_cluster(:, :, i) / M;

end


% Calculate the global range for all heatmaps
global_min = min(min([mean(TF_P1_cluster_norm, 3), mean(TF_P3_cluster_norm, 3), mean(TF_P6_cluster_norm, 3)]));
global_max = max(max([mean(TF_P1_cluster_norm, 3), mean(TF_P3_cluster_norm, 3), mean(TF_P6_cluster_norm, 3)]));
% Set symmetric limits around zero for consistent colormap
global_limit = max(abs([global_min, global_max]));


figure()
tiledlayout(1,4, "TileSpacing", "compact", "Padding", "compact")


ax1 = nexttile;

pcolor(T, freqs, mean(TF_P1_cluster_norm, 3)); hold on;
shading flat; 
xline(T(median_all_flx))
% colormap("turbo")
colormap(ax1, cmap_cont)
clim([-global_limit, global_limit])
set(gca, 'YScale', 'log');  
set(gca, 'FontSize', 14)
yticks([2 4 8 14 30 50])
title('P1', 'FontSize', 16, 'FontWeight', 'normal')
ylabel('Frequency (Hz)')
xlabel('Cycle (%)')

ax2 = nexttile;

pcolor(T, freqs, mean(TF_P3_cluster_norm, 3)); hold on
shading interp; 
xline(T(median_all_flx))
% colormap("turbo")
colormap(ax2, cmap_cont)
clim([-global_limit, global_limit])
set(gca, 'YScale', 'log');  
set(gca, 'FontSize', 14)
yticks([2 4 8 14 30 50])
title('P3', 'FontSize', 16, 'FontWeight', 'normal')
xlabel('Cycle (%)')

ax3 = nexttile;

pcolor(T, freqs, mean(TF_P6_cluster_norm, 3)); hold on
shading interp; 
xline(T(median_all_flx))
% colormap("turbo")
colormap(ax3, cmap_cont)
clim([-global_limit, global_limit])
set(gca, 'YScale', 'log');  
set(gca, 'FontSize', 14)
yticks([2 4 8 14 30 50])
title('P6', 'FontSize', 16, 'FontWeight', 'normal')
xlabel('Cycle (%)')

% cb = colorbar('eastoutside');  
cb = colorbar('southoutside');  
cb.Label.String = 'Normalized Power (dB)';  % Add label
cb.Label.FontSize = 14;


ax4 = nexttile;

pcolor(T_stat, freqs_stat, All_Clusters_TFCE_data.(sets_names{c}).TFCE_map); hold on
shading interp; 
colormap(ax4, cmap_stat);

xline(100*median_all_flx/(median_all_flx + median_all_ext))

% Create contour of significant regions
final_mask = All_Clusters_TFCE_data.(sets_names{c}).final_mask;
[time_grid, freq_grid] = meshgrid(T_stat, freqs_stat);
contour(time_grid, freq_grid, double(final_mask), [0.5 0.5], ...
       'LineWidth', 2, 'Color', 'white');

set(gca, 'YScale', 'log');  
set(gca, 'FontSize', 14)
yticks([2 4 8 14 30 50])

en = char(8211);  % – en dash U+2013
em = char(8212);  % — em dash U+2014
title(ax4, {['RM ANOVA ', en, ' Condition Effect'], '(Contours: p_{FWE} < 0.01)'}, ...
      'Interpreter','tex','FontWeight', 'normal', 'FontSize', 14);

xlabel('Cycle (%)')

% cb = colorbar('eastoutside');  
cb = colorbar('southoutside');  
cb.Label.String = ['TFCE', en, 'Enhanced F', en, 'statistic'];  % Add label
cb.Label.FontSize = 14;
% cb.Ruler.Exponent = 3;              % puts “×10^3” on the colorbar
% cb.Ruler.TickLabelFormat = '%.0f';  % ticks show 0,1,2,... (scaled by 10^3)
% cb.Label.HorizontalAlignment = 'right';

% if isprop(cb.Ruler,'SecondaryLabel')
%     cb.Ruler.SecondaryLabel.String  = '\times10^{3}';
%     cb.Ruler.SecondaryLabel.Visible = 'on';
% end

set(gcf, "Position", [2400 350 1300 500])





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%                       %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%    HELPER FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%                       %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%% Power Spectral Density Calculation

All_Clusters_PSD_data = struct();

for c = 1:numel(sets)

    subjects_components = sets{c};

    cluster_PSD_data = cell(size(subjects_components,1), 4);
    % Extranct TF data for each IC in the cluster
    for i = 1:size(subjects_components, 1)

        subject = subjects_components(i, 1); 
        ic = subjects_components(i, 2);
        disp(['Subject ', num2str(subject), ' IC ', num2str(ic), ' ', sets_names{c}])
        sub_id = sprintf('sub%d', subject);
        trial_info = subjects_trialsInfo.(sub_id);
        subject_data = subjects_data.(sub_id);

        % for baseline calculations:
        EEGLAB_source = ALLEEG(subject-4).icaact(ic, :);
        EEGLAB_events = ALLEEG(subject-4).event;


        [PSD_struct, ~, freqs_baseline] = ...
            main_PSD_calculation_per_IC(subject, ic, trial_info, subject_data, ...
                                  EEGLAB_source, EEGLAB_events);

        cluster_PSD_data{i, 1} = subject;
        cluster_PSD_data{i, 2} = ic;
        cluster_PSD_data{i, 3} = PSD_struct;
        cluster_PSD_data{i, 4} = freqs_baseline;

    end
    cluster_PSD_data = cell2table(cluster_PSD_data, ...
        "VariableNames", ["Subject", "IC", "PSD_Data", "Frequency"]);

    All_Clusters_PSD_data.(sets_names{c}) = cluster_PSD_data;

end


% All_Clusters_PSD_data was saved on Disc on 29.08.2025
% load()



%% plot absolute PSD (dB) data per subject
color_P1 = [0.2, 0.4, 0.8];  % Blue
color_P3 = [0.8, 0.4, 0.2];  % Orange/Red
color_P6 = [0.2, 0.7, 0.3];  % Green

figure()
c = 8;
cluster_PSD_data = All_Clusters_PSD_data.(sets_names{c});
H = height(cluster_PSD_data);
maintile = tiledlayout(3, ceil(H/3)); %, "Padding", "compact", "TileSpacing", "compact");

F = cluster_PSD_data{1, 4}{:};

for i = 1:size(cluster_PSD_data, 1)
    

    PSD_P1 = cat(2, cluster_PSD_data{i, 3}.P1{:});
    PSD_P3 = cat(2, cluster_PSD_data{i, 3}.P3{:});
    PSD_P6 = cat(2, cluster_PSD_data{i, 3}.P6{:});

    baseline_psd = cluster_PSD_data.PSD_Data(i, 1).baseline.P1;  

    
    nexttile()

    plot(F, mean(PSD_P1, 2), 'Color', color_P1, 'LineWidth', 2, 'DisplayName', 'P1'); hold on;
    plot(F, mean(PSD_P3, 2), 'Color', color_P3, 'LineWidth', 2, 'DisplayName', 'P3'); 
    plot(F, mean(PSD_P6, 2), 'Color', color_P6, 'LineWidth', 2, 'DisplayName', 'P6'); 
    plot(F, 10*log10(baseline_psd), 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2, 'DisplayName', 'Baseline');
    xlim([F(1) F(end)])
    
    % set(gca, 'XScale', 'log');  
    set(gca, 'FontSize', 12)
    % xticks([2 4 8 14 30 50])
    title(['S', num2str(cluster_PSD_data{i, 1}), ...
        ' IC', num2str(cluster_PSD_data{i, 2})], ...
        'FontSize', 12, 'FontWeight', 'normal')
    
    if i == 1
        ylabel('PSD (dB)');
    end
    xlabel('Frequency (Hz)')
    
    grid on
    grid minor
    set(gca, 'MinorGridLineStyle', ':')

    if i == size(cluster_PSD_data, 1)
        legend('Location', 'northeast', 'FontSize', 10);
    end
    
end

% title(maintile, 'Left Dorsal ACC', 'FontSize', 16)
title(maintile, strrep(sets_names{c}, '_', ' '), 'FontSize', 16)
% set(gcf, "Position", [1750 -80 2450 1250])






%% plot relative PSD change (dB) data per IC / Baseline: mean log PSD of all conditions
% Per subject:
% rel_PSD(dB) = PSD(dB) - mean_over_all_conditions_and_trials(PSD(dB))

color_P1 = [0.2, 0.4, 0.8];  % Blue
color_P3 = [0.8, 0.4, 0.2];  % Orange/Red
color_P6 = [0.2, 0.7, 0.3];  % Green

figure()
c = 2;
cluster_PSD_data = All_Clusters_PSD_data.(sets_names{c});
H = height(cluster_PSD_data);
maintile = tiledlayout(4, ceil(H/4)); %, "Padding", "compact", "TileSpacing", "compact");

avg_PSD = zeros(1, size(cluster_PSD_data, 1));
for i = 1:size(cluster_PSD_data, 1)
    avg_PSD(i) = mean([cluster_PSD_data.PSD_Data(i,1).P1{:}, ...
        cluster_PSD_data.PSD_Data(i, 1).P3{:}, ...
        cluster_PSD_data.PSD_Data(i, 1).P6{:}], [1, 2]);
end

F = cluster_PSD_data{1, 4}{:};

for i = 1:size(cluster_PSD_data, 1)
    

    PSD_P1 = cat(2, cluster_PSD_data{i, 3}.P1{:});
    PSD_P3 = cat(2, cluster_PSD_data{i, 3}.P3{:});
    PSD_P6 = cat(2, cluster_PSD_data{i, 3}.P6{:});


    PSD_P1 = PSD_P1 - repmat(avg_PSD(i), size(PSD_P1));
    PSD_P3 = PSD_P3 - repmat(avg_PSD(i), size(PSD_P3));
    PSD_P6 = PSD_P6 - repmat(avg_PSD(i), size(PSD_P6));

    
    nexttile()

    plot(F, mean(PSD_P1, 2), 'Color', color_P1, 'LineWidth', 2, 'DisplayName', 'P1'); hold on;
    plot(F, mean(PSD_P3, 2), 'Color', color_P3, 'LineWidth', 2, 'DisplayName', 'P3'); 
    plot(F, mean(PSD_P6, 2), 'Color', color_P6, 'LineWidth', 2, 'DisplayName', 'P6'); 
    xlim([F(1) F(end)])
    
    % set(gca, 'XScale', 'log');  
    set(gca, 'FontSize', 12)
    % xticks([2 4 8 14 30 50])
    title(['S', num2str(cluster_PSD_data{i, 1}), ...
        ' IC', num2str(cluster_PSD_data{i, 2})], ...
        'FontSize', 12, 'FontWeight', 'normal')
    
    if i == 1
        ylabel('Normalized PSD (dB)');
    end
    xlabel('Frequency (Hz)')
    
    grid on
    grid minor
    set(gca, 'MinorGridLineStyle', ':')

    if i == size(cluster_PSD_data, 1)
        legend('Location', 'northeast', 'FontSize', 10);
    end
    
end

% title(maintile, 'Left Dorsal ACC', 'FontSize', 16)
title(maintile, strrep(sets_names{c}, '_', ' '), 'FontSize', 16)
set(gcf, "Position", [1750 -80 2450 1250])



%% plot relative PSD (dB) data per cluster / Baseline: mean log PSD of all conditions

color_P1 = [0.2, 0.4, 0.8];  % Blue
color_P3 = [0.8, 0.4, 0.2];  % Orange/Red
color_P6 = [0.2, 0.7, 0.3];  % Green

c = 2;
cluster_PSD_data = All_Clusters_PSD_data.(sets_names{c});


F = cluster_PSD_data{1, 4}{:};

for i = 1:size(cluster_PSD_data, 1)
    

    PSD_P1 = cat(2, cluster_PSD_data{i, 3}.P1{:});
    PSD_P3 = cat(2, cluster_PSD_data{i, 3}.P3{:});
    PSD_P6 = cat(2, cluster_PSD_data{i, 3}.P6{:});

    if i == 1
        rel_PSD = zeros(size(PSD_P1, 1), 3, size(cluster_PSD_data, 1));
    end


    PSD_P1 = PSD_P1 - repmat(avg_PSD(i), size(PSD_P1));
    PSD_P3 = PSD_P3 - repmat(avg_PSD(i), size(PSD_P3));
    PSD_P6 = PSD_P6 - repmat(avg_PSD(i), size(PSD_P6));

    rel_PSD(:, 1, i) = mean(PSD_P1, 2);
    rel_PSD(:, 2, i) = mean(PSD_P3, 2);
    rel_PSD(:, 3, i) = mean(PSD_P6, 2);
    
    
end




figure()
% H = height(cluster_PSD_data);
% maintile = tiledlayout(4, ceil(H/4)); %, "Padding", "compact", "TileSpacing", "compact");


% Mean
h1 = plot(F, mean(rel_PSD(:, 1, :), 3), 'Color', color_P1, 'LineWidth', 2); hold on;
h2 = plot(F, mean(rel_PSD(:, 2, :), 3), 'Color', color_P3, 'LineWidth', 2); 
h3 = plot(F, mean(rel_PSD(:, 3, :), 3), 'Color', color_P6, 'LineWidth', 2); 

% Upper STD
plot(F, mean(rel_PSD(:, 1, :), 3) + std(rel_PSD(:, 1, :), 0, 3), 'Color', color_P1, 'LineWidth', 1, 'LineStyle', '--'); 
plot(F, mean(rel_PSD(:, 2, :), 3) + std(rel_PSD(:, 2, :), 0, 3), 'Color', color_P3, 'LineWidth', 1, 'LineStyle', '--'); 
plot(F, mean(rel_PSD(:, 3, :), 3) + std(rel_PSD(:, 3, :), 0, 3), 'Color', color_P6, 'LineWidth', 1, 'LineStyle', '--');
% Lower STD
plot(F, mean(rel_PSD(:, 1, :), 3) - std(rel_PSD(:, 1, :), 0, 3), 'Color', color_P1, 'LineWidth', 1, 'LineStyle', '--'); 
plot(F, mean(rel_PSD(:, 2, :), 3) - std(rel_PSD(:, 2, :), 0, 3), 'Color', color_P3, 'LineWidth', 1, 'LineStyle', '--'); 
plot(F, mean(rel_PSD(:, 3, :), 3) - std(rel_PSD(:, 3, :), 0, 3), 'Color', color_P6, 'LineWidth', 1, 'LineStyle', '--');



% xlim([F(1) F(end)])
xlim([2 40])

% set(gca, 'XScale', 'log');  
set(gca, 'FontSize', 12)
xticks([4 8 14 30])
title(strrep(sets_names{c}, '_', ' '), 'FontSize', 16, 'FontWeight', 'normal')

ylabel('Normalized PSD (dB)');
xlabel('Frequency (Hz)')

grid on
grid minor
set(gca, 'MinorGridLineStyle', ':')

legend([h1 h2 h3], {'P1 (Mean $\pm$ Std)', 'P3 (Mean $\pm$ Std)', 'P6 (Mean $\pm$ Std)'}, ...
    'Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex');


% set(gcf, "Position", [1750 -80 2450 1250])