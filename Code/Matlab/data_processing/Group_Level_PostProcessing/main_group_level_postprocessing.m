clc
clear

%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
processed_data_path = [data_path, '5_single-subject-EEG-analysis\'];
study_folder = [data_path, '7_STUDY'];


%% Run EEGLAB
if ~exist('ALLCOM','var')
    eeglab;
end


%% Create and save main STUDY Files. 
% If Studies are created before just "Load existing study" with EEGLAB GUI
% (1) main_study_all_ICs_RV-15.study
% (2) main_study_brain_ICs_RV-15.study

% Define subjects
subject_list = [5, 6, 7, 8, 9, 10];  % List of subject IDs

%% calling function to create and save the main study files
% create_and_save_main_study_files(subject_list, data_path, processed_data_path, ALLEEG)


%% Compute the pre-cluster
clustering_weights = struct();
clustering_weights.dipoles = 3;
clustering_weights.scalp = 1;

clustering_weights.spec = 0;
freqrange = [1 50];
timewindow = [];

out_filename = [STUDY.name(1:end-6), '_preclustered'];
out_filepath = study_folder;

[STUDY, ALLEEG, EEG] = ...
    bemobil_precluster(STUDY, ALLEEG, EEG, clustering_weights, ...
    freqrange, timewindow, out_filename, out_filepath);


%% Implementing the repeated clustering

% select the region of interest (ROI) - MNI Coordinates
% Left_PrimMotor = [-36, -19, 48];
% Left_PreMot_SuppMot = [-28, -2, 52]; % ref: https://bioimagesuiteweb.github.io/webapp/mni2tal.html
% Left_Paracentral_Lobule = [0, -20, 62]; % ref: https://pmc.ncbi.nlm.nih.gov/articles/PMC5663902/table/Tab2/
% Left_Dorsal_ACC = [-5, 39, 20]; % ACC: Anterior Cingulate Cortex - ref: https://bioimagesuiteweb.github.io/webapp/mni2tal.html
% Left_VisMotor = [-18, -67, 40]; % ref: https://bioimagesuiteweb.github.io/webapp/mni2tal.html
% Left_PrimVisual = [-11, -81, 7]; % ref: https://bioimagesuiteweb.github.io/webapp/mni2tal.html
% Right_PreMot_SuppMot = [28, -1, 51]; % ref: https://bioimagesuiteweb.github.io/webapp/mni2tal.html
% Right_VisMotor = [23, -60, 61]; % ref: https://bioimagesuiteweb.github.io/webapp/mni2tal.html
Right_PrimVisual = [11, -78, 9]; % ref: https://bioimagesuiteweb.github.io/webapp/mni2tal.html

ROI = Right_PrimVisual;
cluster_ROI_name  = 'Right_PrimVisual';
cluster_ROI_MNI.x = ROI(1);
cluster_ROI_MNI.y = ROI(2);
cluster_ROI_MNI.z = ROI(3);


% standard deviations boundary for outlier detection
outlier_sigma = 3; 

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
filename_STUDY = [STUDY.name(1:end-6), '_', cluster_ROI_name];

% filepath where the clustering solutions should be saved
filepath_clustering_solutions = [filepath_STUDY, '\clustering_solutions'];
% filename of the repeated clustering solutions' data set
filename_clustering_solutions = [cluster_ROI_name, '_clustering_solutions'];

% filepath where the multivariate data should be saved
filepath_multivariate_data = [filepath_STUDY, '\multivariate_data'];
% filename of the multivariate data set
filename_multivariate_data = [cluster_ROI_name, '_multivariate_data'];


N = floor(size(STUDY.etc.preclust.preclustdata, 1) / size(subject_list, 2));
number_of_clusters = [N-2:N+2];
score_NCluster = zeros(length(number_of_clusters),2);

for i = 1:numel(number_of_clusters)

    n_clust = number_of_clusters(i);
    
    [STUDY, ALLEEG, EEG] = ...
        repeated_clustering_and_evaluation_custom(STUDY, ALLEEG, EEG, ...
        outlier_sigma, n_clust, n_iterations, cluster_ROI_MNI, ...
        quality_measure_weights, do_clustering, do_multivariate_data, ...
        filepath_STUDY, filename_STUDY, ...
        filepath_clustering_solutions, filename_clustering_solutions, ...
        filepath_multivariate_data, filename_multivariate_data);
    
    
    cd([study_folder, '\multiple_clustering\', cluster_ROI_name, '\clustering_solutions'])
    load([cluster_ROI_name, '_multivariate_data'])
    
    [ranked_scores, ranked_solutions] = clustering_rank_solutions_custom(cluster_multivariate_data,quality_measure_weights);
    
    score_NCluster(i, :) = [ranked_scores(1), n_clust];

end


% n_clust = 
% 
[STUDY, ALLEEG, EEG] = ...
    bemobil_repeated_clustering_and_evaluation(STUDY, ALLEEG, EEG, ...
    outlier_sigma, score_NCluster(1, 2), n_iterations, cluster_ROI_MNI, ...
    quality_measure_weights, do_clustering, do_multivariate_data, ...
    filepath_STUDY, filename_STUDY, ...
    filepath_clustering_solutions, filename_clustering_solutions, ...
    filepath_multivariate_data, filename_multivariate_data);

%% End of code - Remove the folder and its subfolders from the MATLAB path
rmpath(main_project_folder);