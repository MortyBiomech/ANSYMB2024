clc
clear

%% Add and Define Necessary Paths
eeglab_path = 'C:\Morteza\Toolboxes\EEGLAB\eeglab2024.2.1';
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
processed_data_path = [data_path, '5_single-subject-EEG-analysis\'];
study_folder = [data_path, '7_STUDY'];


%% Load ALLEEG 
if ~exist('ALLCOM','var')
    eeglab;
end

% Define subjects
subject_list = 5:18;  % List of subject IDs

for i = 1:length(subject_list)
    file_name = ['sub-', num2str(subject_list(i)), '_cleaned_with_ICA.set'];
    dataset_path = fullfile([processed_data_path, 'sub-', ...
        num2str(subject_list(i)), filesep]);
    EEG = pop_loadset('filename', file_name, 'filepath', dataset_path);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, i);
end


%% Load STUDY files
cd(study_folder)
load("all_STUDY_names.mat")
load("all_STUDY_files.mat")


%% Plot dipoles in 3D space with MRI slices
% Vibrant Colors
colors = [
    [1.0, 0.0, 0.0];   % Red
    [0.0, 0.0, 1.0];   % Blue
    [0.0, 1.0, 0.0];   % Green
    [1.0, 1.0, 0.0];   % Yellow
    [1.0, 0.5, 0.0];   % Orange
    [0.5, 0.0, 0.5];   % Purple
    [0.0, 1.0, 1.0];   % Cyan
    [1.0, 0.0, 1.0];   % Magenta
];

% Specify a high-resolution MRI template
high_res_mri = fullfile(eeglab_path, 'plugins', 'dipfit', 'standard_BEM', 'standard_mri.mat');
load(high_res_mri);

% Initialize dipole structure
dipole_struct = [];
color = [];
dipolesize = [];

row = 1;

for i = 1:length(all_STUDY_files)
    
    disp(all_STUDY_names{i})
    best_cluster_id = all_STUDY_files{1, i}.etc.bemobil.clustering.cluster_ROI_index;
    subjects = all_STUDY_files{1, i}.cluster(best_cluster_id).sets;
    ICs = all_STUDY_files{1, i}.cluster(best_cluster_id).comps;
    
    for j = 1:length(subjects)
        
        source = ALLEEG(subjects(j)).dipfit.model(ICs(j));
        dipole_struct(row).posxyz = source.posxyz;
        dipole_struct(row).momxyz = source.momxyz;
        dipole_struct(row).rv = source.rv;

        color(row, :) = colors(i, :);
        dipolesize(row, :) = 30;

        row = row + 1;
        
    end

end


% Convert HEX to RGB
hex2rgb = @(hex) sscanf(hex(2:end), '%2x%2x%2x', [1 3]) / 255;

% Darker vibrant colors in RGB
dark_colors = [
    hex2rgb('#B20000'); % Dark Red
    hex2rgb('#0000B2'); % Dark Blue
    hex2rgb('#008000'); % Dark Green
    hex2rgb('#B2B200'); % Dark Yellow
    hex2rgb('#B25900'); % Dark Orange
    hex2rgb('#4D004D'); % Dark Purple
    hex2rgb('#008080'); % Dark Cyan
    hex2rgb('#B200B2'); % Dark Magenta
];

for i = 1:length(all_STUDY_files)
    
    best_cluster_id = all_STUDY_files{1, i}.etc.bemobil.clustering.cluster_ROI_index;
    
    source = all_STUDY_files{1, i}.cluster(best_cluster_id).dipole;
    dipole_struct(row).posxyz = source.posxyz;
    dipole_struct(row).momxyz = source.momxyz;
    dipole_struct(row).rv = source.rv;

    color(row, :) = dark_colors(i, :);
    dipolesize(row, :) = 50;

    row = row + 1;
        
end

color = mat2cell(color, ones(size(color,1), 1), 3);
%'color', colors(i, :), ...  % Dipole color
            
%%
figure();
dipplot(dipole_struct, ...
    'coordformat', 'MNI', ...   % Use MNI coordinates
    'mri', high_res_mri, ...    % Matlab file containing an MRI volume
    'image', 'mri', ...         % Use default MRI template
    'view', [0, -1, 0], ...      % sagittal view
    'gui', 'off', ...           % Display controls
    'projlines', 'off', ...     % Project lines to slices
    'projimg', 'off', ...       % Project images to slices
    'dipolesize', dipolesize, ...       % Size of the dipole sphere(s)
    'dipolelength', 0, ...      % Length of the dipole bar
    'color', color, ...
    'spheres', 'on');           % Show dipoles as spheres