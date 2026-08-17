clc
clear


%% Add paths
addpath(genpath('D:\Morteza\MyProjects\ANSYMB2024\Code'))
addpath(genpath('D:\Morteza\MyProjects\ANSYMB2024\data\7_STUDY\Epoched_data'))
data_path = 'D:\Morteza\MyProjects\ANSYMB2024\data\';
ersp_data_path = [data_path, '7_STUDY\Epoched_data\Final_figures', ...
    '\ERSP\Three Pressure Conditions\p 0.01 ersp results\'];
study_data_path = [data_path, '7_STUDY\Epoched_data\multiple_clustering'];
eeglabpath = 'D:\Morteza\Toolboxes\EEGLAB\eeglab2026.0.0';


%% Study names
% studyNames = {'Left_Dorsal_ACC', 'Left_Parieto_Occipital', ...
%     'Left_PreMot_SuppMot', 'Left_Prim_Motor', 'Prime_Visual', ...
%     'Right_Parieto_Occipital', 'Right_PreMot_SuppMot', 'Right_Prim_Motor'};

studyNames = {'Prime_Visual'};


%% run eeglab
current_path = pwd;
if ~exist("ALLEEG", "var")
    cd(eeglabpath)
    eeglab
end
cd(current_path)


%% Load Study files and store the Subjects and ICs
inner_struct = struct('Subjects', [], 'ICs', []);
SUBJECTS_ICS = cell(length(studyNames), 2);
SUBJECTS_ICS(:, 1) = studyNames;
SUBJECTS_ICS(:, 2) = repmat({inner_struct}, length(studyNames), 1);

for s = 1:length(studyNames)

    fileName = [studyNames{s}, '.study'];
    filePath = [study_data_path, filesep, studyNames{s}];
    
    [STUDY ALLEEG] = pop_loadstudy('filename', fileName, ...
        'filepath', filePath);
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
    eeglab redraw;
    
    STUDY = std_makedesign(STUDY, ALLEEG, 1, ...
        'name', '3-condition design', ...
        'delfiles', 'off', 'defaultdesign', 'off', ...
        'variable1', 'cond', 'values1', {'1','3','6'}, ...
        'vartype1','categorical', ...
        'pairing','on'); % within-subject since every subject did all conds
    
    STUDY = add_anatomical_labels(STUDY, eeglabpath);
    
    % adjust clusters so there's one subject per cluster using 
    % lowest IC number (highest explained standard deviation)
    STUDY = oneSubPerCluster(STUDY); 
    
    CL = STUDY.etc.bemobil.clustering.cluster_ROI_index;
    subjects_in_cluster = STUDY.cluster(CL).sets;
    ICs_in_cluster = STUDY.cluster(CL).comps;

    SUBJECTS_ICS{s, 2}.Subjects = subjects_in_cluster;
    SUBJECTS_ICS{s, 2}.ICs = ICs_in_cluster;

    
end
