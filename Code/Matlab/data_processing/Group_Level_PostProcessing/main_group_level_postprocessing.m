clc
clear

%% Add and Define Necessary Paths
addpath(genpath('C:\Morteza\MyProjects\ANSYMB2024')); % main folder containing all codes and data

% Change path to the directory on your PC which raw XDF files are stored:
data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
processed_data_path = [data_path, '5_single-subject-EEG-analysis\'];


%% Run EEGLAB
if ~exist('ALLCOM','var')
	eeglab;
end


%% Create a new STUDY set

% Define subjects
subject_list = {'5', '6', '7', '8', '9', '10'};  % List of subject IDs

% Load datasets into ALLEEG
for i = 1:length(subject_list)
    file_name = ['sub-', subject_list{i}, '_cleaned_with_ICA.set'];
    dataset_path = fullfile([processed_data_path, 'sub-', ...
        subject_list{i}, filesep]);
    EEG = pop_loadset('filename', file_name, 'filepath', dataset_path);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, i);
end


% Create STUDY
study_name = 'main_study.study';
study_folder = [data_path, '7_STUDY'];

% Create STUDY from loaded datasets using commands
commands = cell(1, length(subject_list));
for i = 1:length(subject_list)
    commands{i} = {'index', i, 'subject', subject_list{i}, ...
        'inbrain', 'on', 'dipselect', 0.15};
end

[STUDY, ALLEEG] = std_editset([], ALLEEG, ...
    'name', study_name, ...
    'filepath', study_folder, ...
    'commands', commands, ...
    'updatedat', 'on', 'savedat', 'off');

% Save STUDY
[STUDY, ALLEEG] = pop_savestudy(STUDY, ALLEEG, 'filename', study_name, 'filepath', study_folder);


%% Compute the pre-cluster
clustering_weights = struct();
clustering_weights.dipoles = 3;
clustering_weights.scalp = 1;

clustering_weights.spec = 0;
freqrange = [1 50];
timewindow = [];

out_filename = 'main_study_preclustered';
out_filepath = study_folder;

[STUDY, ALLEEG, EEG] = ...
    bemobil_precluster(STUDY, ALLEEG, EEG, clustering_weights, ...
    freqrange, timewindow, out_filename, out_filepath);


%% Implementing the repeated clustering

% select the region of interest (ROI) - MNI Coordinates
Left_PrimMotor = [-36, -19, 48];

cluster_ROI_name  = 'Left_PrimMotor';
cluster_ROI_MNI.x = Left_PrimMotor(1);
cluster_ROI_MNI.y = Left_PrimMotor(2);
cluster_ROI_MNI.z = Left_PrimMotor(3);


% standard deviations boundary for outlier detection
outlier_sigma = 3; 

% number of clusters to be created
n_clust = 14;

% number of iterations to be performed for repeated clustering
n_iterations = 1000;

% quality_measure_weights: vector of weights for quality measures (6 entries): 
%                          - subjects, 
%                          - ICs/subjects, 
%                          - normalized spread, 
%                          - mean RV, 
%                          - distance from ROI, 
%                          - mahalanobis distance from median of multivariate distribution
%                            (put this very high to get the most "normal" solution)
quality_measure_weights = [3, -1, -1, -1, -2, -1];

% whether or not the clustering should be done (it takes a lot of time).
do_clustering = 1;

% whether or not the creation of the multivariate dataset should be done.
do_multivariate_data = 1;

% filepath where the new STUDY should be saved
filepath_STUDY = [study_folder, '\multiple_clustering\', cluster_ROI_name];
% filename of the new STUDY
filename_STUDY = ['main_study_', cluster_ROI_name];

% filepath where the clustering solutions should be saved
filepath_clustering_solutions = [filepath_STUDY, '\clustering_solutions'];
% filename of the repeated clustering solutions' data set
filename_clustering_solutions = [cluster_ROI_name, '_clustering_solutions'];

% filepath where the multivariate data should be saved
filepath_multivariate_data = [filepath_STUDY, '\multivariate_data'];
% filename of the multivariate data set
filename_multivariate_data = [cluster_ROI_name, '_multivariate_data'];


[STUDY, ALLEEG, EEG] = ...
    bemobil_repeated_clustering_and_evaluation(STUDY, ALLEEG, EEG, ...
    outlier_sigma, n_clust, n_iterations, cluster_ROI_MNI, ...
    quality_measure_weights, do_clustering, do_multivariate_data, ...
    filepath_STUDY, filename_STUDY, ...
    filepath_clustering_solutions, filename_clustering_solutions, ...
    filepath_multivariate_data, filename_multivariate_data);




