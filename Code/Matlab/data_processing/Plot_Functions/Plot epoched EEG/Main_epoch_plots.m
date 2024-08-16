function Main_epoch_plots(subject_id)
    
    % Inputs:
    % subjet_id  = subject id number
    % std_sem    = option for plot. Must be either 'std' standard deviation
    %              or 'sem' (standard error of mean). This input will be
    %              specified in the App by the user

    %% Add and Define Necessary Paths
    addpath(genpath('C:\Morteza\MyProjects\ANSYMB2024')); % main folder containing all codes and data
    
    % Change data_path according to your system:
    data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
    epochs_path = [data_path, '6_Trials_Info_and_Epoched_data', ...
        filesep, 'sub-', num2str(subject_id)];

    % Create a UI figure for loading screen
    loadingFig = uifigure('Units', 'normalized', 'Position', [0.4, 0.4, 0.4, 0.3], 'Name', '');

    % Create a progress bar
    progressBar = uiprogressdlg(loadingFig, 'Title', 'Loading Epoched Datasets', ...
        'Message', 'Loading ...', 'Indeterminate', 'off', 'Value', 0);

    pause(1); % Simulate loading time
    progressBar.Message = sprintf('Loading Epochs_Trial_based.mat');
    a1 = load(fullfile(epochs_path, 'Epochs_Trial_based.mat'));
    Epochs_Trial_based = a1.Epochs_Trial_based;
    progressBar.Value = 1 / 4;

    progressBar.Message = sprintf('Loading Epochs_Flexion_based.mat');
    a1 = load(fullfile(epochs_path, 'Epochs_Flexion_based.mat'));
    Epochs_Flexion_based = a1.Epochs_Flexion_based;
    progressBar.Value = 2 / 4;

    progressBar.Message = sprintf('Loading Epochs_Extension_based.mat');
    a1 = load(fullfile(epochs_path, 'Epochs_Extension_based.mat'));
    Epochs_Extension_based = a1.Epochs_Extension_based;
    progressBar.Value = 3 / 4;

    pause(1); % Simulate loading time
    progressBar.Message = sprintf('Loading Epochs_FlextoFlex_based.mat');    
    a1 = load(fullfile(epochs_path, 'Epochs_FlextoFlex_based.mat'));
    Epochs_FlextoFlex_based = a1.Epochs_FlextoFlex_based;
    progressBar.Value = 4 / 4;

    % Close the loading figure
    close(loadingFig);


    %% Select the trials with the same Pressure
    a = load(fullfile(epochs_path, 'Trials_Info.mat'));
    Trials_Info = a.Trials_Info;

    P1 = struct();
    P1.trials = [];
    P1.scores = [];
    
    P3 = struct();
    P3.trials = [];
    P3.scores = [];
    
    P6 = struct();
    P6.trials = [];
    P6.scores = [];
    
    for i = 1:length(Trials_Info)
        p = Trials_Info{1, i}.General.Pressure;
        switch p
            case 1
                P1.trials(1, end+1) = i;
                P1.scores(1, end+1) = Trials_Info{1, i}.General.Score;
            case 3
                P3.trials(1, end+1) = i;
                P3.scores(1, end+1) = Trials_Info{1, i}.General.Score;
            case 6
                P6.trials(1, end+1) = i;
                P6.scores(1, end+1) = Trials_Info{1, i}.General.Score;
        end
                
    end


    % cd('C:\Morteza\Analysis\ANSYMB2024\Code\data_processing\Plot functions\Plot epoched EEG')



    %% Create the main UI figure for selection plot options
    fig = uifigure('Units', 'normalized','Position', [0.2 0.4 0.5 0.32], 'Name', 'Plot EEG Epoched data');

    % Create radio buttons for options
    bg = uibuttongroup(fig, 'Units', 'normalized', 'Position', [0.032, 0.18, 0.936, 0.75], 'Title', 'Select to plot:');

    d = [0 -70 0 0];
    uiradiobutton(bg, 'Position', [10, 240, 300, 22]+d, 'Text', 'EEG entire trials - Time Domain - Channels');
    uiradiobutton(bg, 'Position', [10, 220, 300, 22]+d, 'Text', 'EEG entire trials - Time Domain - ICs');
    uiradiobutton(bg, 'Position', [10, 200, 300, 22]+d, 'Text', 'EEG entire trials - Frequency Domain - Channels');
    uiradiobutton(bg, 'Position', [10, 180, 300, 22]+d, 'Text', 'EEG entire trials - Frequency Domain - ICs');

    uiradiobutton(bg, 'Position', [10, 140, 300, 22]+d, 'Text', 'EEG flexion epochs - Time Domain - Channels');
    uiradiobutton(bg, 'Position', [10, 120, 300, 22]+d, 'Text', 'EEG flexion epochs - Time Domain - ICs');
    uiradiobutton(bg, 'Position', [10, 100, 300, 22]+d, 'Text', 'EEG flexion epochs - Frequency Domain - Channels');
    uiradiobutton(bg, 'Position', [10, 80, 300, 22]+d, 'Text', 'EEG flexion epochs - Frequency Domain - ICs');

    uiradiobutton(bg, 'Position', [350, 240, 400, 22]+d, 'Text', 'EEG extension epochs - Time Domain - Channels');
    uiradiobutton(bg, 'Position', [350, 220, 400, 22]+d, 'Text', 'EEG extension epochs - Time Domain - ICs');
    uiradiobutton(bg, 'Position', [350, 200, 400, 22]+d, 'Text', 'EEG extension epochs - Frequency Domain - Channels');
    uiradiobutton(bg, 'Position', [350, 180, 400, 22]+d, 'Text', 'EEG extension epochs - Frequency Domain - ICs');

    uiradiobutton(bg, 'Position', [350, 140, 400, 22]+d, 'Text', 'EEG flex2flex epochs - Time Domain - Channels');
    uiradiobutton(bg, 'Position', [350, 120, 400, 22]+d, 'Text', 'EEG flex2flex epochs - Time Domain - ICs');
    uiradiobutton(bg, 'Position', [350, 100, 400, 22]+d, 'Text', 'EEG flex2flex epochs - Frequency Domain - Channels');
    uiradiobutton(bg, 'Position', [350, 80, 400, 22]+d, 'Text', 'EEG flex2flex epochs - Frequency Domain - ICs');


    % Create std/sem slider
    uilabel(fig, 'Position', [35, 20, 260, 22], 'Text', 'Mean\pm', 'Interpreter', 'latex', 'VerticalAlignment', 'center');
    dd = uidropdown(fig, 'Position', [90, 20, 220, 22], 'Items', ...
        {'SEM (Standard Error of the Mean)', 'STD (Standard Deviation)'}, ...
        'Value', 'SEM (Standard Error of the Mean)');


    % Create Plot button
    uibutton(fig, 'Position', [600, 20, 100, 22], 'Text', 'Plot', ...
        'ButtonPushedFcn', @(btn, event) plotCallback(bg, dd.Value(1:3)));

    % Callback function
    function plotCallback(bg, statMethod)
        selectedOption = bg.SelectedObject.Text;
        switch selectedOption

            % Entire Trials
            case 'EEG entire trials - Time Domain - Channels'
                disp(['EEG entire trials - Time Domain - Channels (with ', statMethod, ')']);
                data = Epochs_Trial_based;
                EEG_entire_trials_Time_Domain_channels(data, statMethod, P1, P3, P6);
            case 'EEG entire trials - Time Domain - ICs'
                disp(['EEG entire trials - Time Domain - ICs (with ', statMethod, ')']);
                data = Epochs_Trial_based;
                EEG_entire_trials_Time_Domain_sources(data, statMethod, P1, P3, P6);
            case 'EEG entire trials - Frequency Domain - Channels'
                disp(['EEG entire trials - Frequency Domain - Channels (with ', statMethod, ')']);
                data = Epoches_Trial_based;
                EEG_entire_trials_Frequency_Domain_channels(data, statMethod, P1, P3, P6);
            case 'EEG entire trials - Frequency Domain - ICs'
                disp(['EEG entire trials - Frequency Domain - ICs (with ', statMethod, ')']);
                data = Epochs_Trial_based;
                EEG_entire_trials_Frequency_Domain_sources(data, statMethod, P1, P3, P6);

            % Flexion Epochs
            case 'EEG flexion epochs - Time Domain - Channels'
                disp(['EEG flexion epochs - Time Domain - Channels (with ', statMethod, ')']);
                data = Epochs_Flexion_based;
                EEG_flexion_epochs_Time_Domain_channels(data, statMethod, P1, P3, P6);
            case 'EEG flexion epochs - Time Domain - ICs'
                disp(['EEG flexion epochs - Time Domain - ICs (with ', statMethod, ')']);
                data = Epochs_Flexion_based;
                EEG_flexion_epochs_Time_Domain_sources(data, statMethod, P1, P3, P6);
            case 'EEG flexion epochs - Frequency Domain - Channels'
                disp(['EEG flexion epochs - Frequency Domain - Channels (with ', statMethod, ')']);
                data = Epochs_Flexion_based;
                EEG_flexion_epochs_Frequency_Domain_channels(data, statMethod, P1, P3, P6);
            case 'EEG flexion epochs - Frequency Domain - ICs'
                disp(['EEG flexion epochs - Frequency Domain - ICs (with ', statMethod, ')']);
                data = Epochs_Flexion_based;
                EEG_flexion_epochs_Frequency_Domain_sources(data, statMethod, P1, P3, P6);

            % Extension Epochs
            case 'EEG extension epochs - Time Domain - Channels'
                disp(['EEG extension epochs - Time Domain - Channels (with ', statMethod, ')']);
                data = Epochs_Extension_based;
                EEG_extension_epochs_Time_Domain_channels(data, statMethod, P1, P3, P6);
            case 'EEG extension epochs - Time Domain - ICs'
                disp(['EEG extension epochs - Time Domain - ICs (with ', statMethod, ')']);
                data = Epochs_Extension_based;
                EEG_extension_epochs_Time_Domain_sources(data, statMethod, P1, P3, P6);
            case 'EEG extension epochs - Frequency Domain - Channels'
                disp(['EEG extension epochs - Frequency Domain - Channels (with ', statMethod, ')']);
                data = Epochs_Extension_based;
                EEG_extension_epochs_Frequency_Domain_channels(data, statMethod, P1, P3, P6);
            case 'EEG extension epochs - Frequency Domain - ICs'
                disp(['EEG extension epochs - Frequency Domain - ICs (with ', statMethod, ')']);
                data = Epochs_Extension_based;
                EEG_extension_epochs_Frequency_Domain_sources(data, statMethod, P1, P3, P6);

            % Flexion-to-Flexion Epochs
            case 'EEG flex2flex epochs - Time Domain - Channels'
                disp(['EEG flex2flex epochs - Time Domain - Channels (with ', statMethod, ')']);
                data = Epochs_FlextoFlex_based;
                EEG_flextoflex_epochs_Time_Domain_channels(data, statMethod, P1, P3, P6);
            case 'EEG flex2flex epochs - Time Domain - ICs'
                disp(['EEG flex2flex epochs - Time Domain - ICs (with ', statMethod, ')']);
                data = Epochs_FlextoFlex_based;
                EEG_flextoflex_epochs_Time_Domain_sources(data, statMethod, P1, P3, P6);
            case 'EEG flex2flex epochs - Frequency Domain - Channels'
                disp(['EEG flex2flex epochs - Frequency Domain - Channels (with ', statMethod, ')']);
                data = Epochs_FlextoFlex_based;
                EEG_flextoflex_epochs_Frequency_Domain_channels(data, statMethod, P1, P3, P6);
            case 'EEG flex2flex epochs - Frequency Domain - ICs'
                disp(['EEG flex2flex epochs - Frequency Domain - ICs (with ', statMethod, ')']);
                data = Epochs_FlextoFlex_based;
                EEG_flextoflex_epochs_Frequency_Domain_sources(data, statMethod, P1, P3, P6);

            otherwise
                disp('Unknown option');
        end

    end
end