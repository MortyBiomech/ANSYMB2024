function Main_epoch_plots(subject)
    
    % Inputs:
    % std_sem    = option for plot. Must be either 'std' standard deviation
    %              or 'sem' (standard error of mean).
    
    filepath = ['C:\Morteza\Analysis\ANSYMB2024\data\', ...
        '5_single-subject-EEG-analysis_with_epochs\', ...
        'sub-', num2str(subject)];
    cd(filepath)

    % Create a UI figure for loading screen
    loadingFig = uifigure('Units', 'normalized', 'Position', [0.4, 0.4, 0.4, 0.3], 'Name', '');

    % Create a progress bar
    progressBar = uiprogressdlg(loadingFig, 'Title', 'Loading Data', ...
        'Message', 'Loading datasets...', 'Indeterminate', 'off', 'Value', 0);

    pause(1); % Simulate loading time
    progressBar.Message = sprintf('Loading dataset: EEG_epoched_trials.mat');
    EEG_epoched_trials = load('EEG_epoched_trials.mat');
    EEG_epoched_trials = EEG_epoched_trials.EEG_epoched_trials;
    progressBar.Value = 1 / 4;

    pause(1); % Simulate loading time
    progressBar.Message = sprintf('Loading dataset: EEG_epoched_flexion.mat');
    EEG_epoched_flexion = load('EEG_epoched_flexion.mat');
    EEG_epoched_flexion = EEG_epoched_flexion.EEG_epoched_flexion;
    progressBar.Value = 2 / 4;

    pause(1); % Simulate loading time
    progressBar.Message = sprintf('Loading dataset: EEG_epoched_extension.mat');
    EEG_epoched_extension = load('EEG_epoched_extension.mat');
    EEG_epoched_extension = EEG_epoched_extension.EEG_epoched_extension;
    progressBar.Value = 3 / 4;

    pause(1); % Simulate loading time
    progressBar.Message = sprintf('Loading dataset: EEG_epoched_flextoflex.mat');    
    EEG_epoched_flextoflex = load('EEG_epoched_flextoflex.mat');
    EEG_epoched_flextoflex = EEG_epoched_flextoflex.EEG_epoched_flextoflex;
    progressBar.Value = 4 / 4;

    % Close the loading figure
    close(loadingFig);

    cd('C:\Morteza\Analysis\ANSYMB2024\Code\data_processing\Plot functions\Plot epoched EEG')



    % Create a UI figure
    fig = uifigure('Units', 'normalized','Position', [0.2 0.4 0.5 0.32], 'Name', 'Select an Option');
    
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


    % Create a slider
    uilabel(fig, 'Position', [35, 20, 260, 22], 'Text', 'Mean\pm', 'Interpreter', 'latex', 'VerticalAlignment', 'center');
    dd = uidropdown(fig, 'Position', [90, 20, 220, 22], 'Items', {'STD (Standard Deviation)', 'SEM (Standard Error of the Mean)'});


    % Create a button
    uibutton(fig, 'Position', [600, 20, 100, 22], 'Text', 'Plot', ...
        'ButtonPushedFcn', @(btn, event) plotCallback(bg, dd.Value(1:3)));

    % Callback function
    function plotCallback(bg, statMethod)
        selectedOption = bg.SelectedObject.Text;
        switch selectedOption
            
            % Entire Trials
            case 'EEG entire trials - Time Domain - Channels'
                disp(['EEG entire trials - Time Domain - Channels (with ', statMethod, ')']);
                data = EEG_epoched_trials;
                EEG_entire_trials_Time_Domain_channels(data, statMethod);
            case 'EEG entire trials - Time Domain - ICs'
                disp(['EEG entire trials - Time Domain - ICs (with ', statMethod, ')']);
                data = EEG_epoched_trials;
                EEG_entire_trials_Time_Domain_sources(data, statMethod);
            case 'EEG entire trials - Frequency Domain - Channels'
                disp(['EEG entire trials - Frequency Domain - Channels (with ', statMethod, ')']);
                data = EEG_epoched_trials;
                EEG_entire_trials_Frequency_Domain_channels(data, statMethod);
            case 'EEG entire trials - Frequency Domain - ICs'
                disp(['EEG entire trials - Frequency Domain - ICs (with ', statMethod, ')']);
                data = EEG_epoched_trials;
                EEG_entire_trials_Frequency_Domain_sources(data, statMethod);
            
            % Flexion Epochs
            case 'EEG flexion epochs - Time Domain - Channels'
                disp(['EEG flexion epochs - Time Domain - Channels (with ', statMethod, ')']);
                data = EEG_epoched_flexion;
                EEG_flexion_epochs_Time_Domain_channels(data, statMethod);
            case 'EEG flexion epochs - Time Domain - ICs'
                disp(['EEG flexion epochs - Time Domain - ICs (with ', statMethod, ')']);
                data = EEG_epoched_flexion;
                EEG_flexion_epochs_Time_Domain_sources(data, statMethod);
            case 'EEG flexion epochs - Frequency Domain - Channels'
                disp(['EEG flexion epochs - Frequency Domain - Channels (with ', statMethod, ')']);
                data = EEG_epoched_flexion;
                EEG_flexion_epochs_Frequency_Domain_channels(data, statMethod);
            case 'EEG flexion epochs - Frequency Domain - ICs'
                disp(['EEG flexion epochs - Frequency Domain - ICs (with ', statMethod, ')']);
                data = EEG_epoched_flexion;
                EEG_flexion_epochs_Frequency_Domain_sources(data, statMethod);
            
            % Extension Epochs
            case 'EEG extension epochs - Time Domain - Channels'
                disp(['EEG extension epochs - Time Domain - Channels (with ', statMethod, ')']);
                data = EEG_epoched_extension;
                EEG_extension_epochs_Time_Domain_channels(data, statMethod);
            case 'EEG extension epochs - Time Domain - ICs'
                disp(['EEG extension epochs - Time Domain - ICs (with ', statMethod, ')']);
                data = EEG_epoched_extension;
                EEG_extension_epochs_Time_Domain_sources(data, statMethod);
            case 'EEG extension epochs - Frequency Domain - Channels'
                disp(['EEG extension epochs - Frequency Domain - Channels (with ', statMethod, ')']);
                data = EEG_epoched_extension;
                EEG_extension_epochs_Frequency_Domain_channels(data, statMethod);
            case 'EEG extension epochs - Frequency Domain - ICs'
                disp(['EEG extension epochs - Frequency Domain - ICs (with ', statMethod, ')']);
                data = EEG_epoched_extension;
                EEG_extension_epochs_Frequency_Domain_sources(data, statMethod);
            
            % Flexion-to-Flexion Epochs
            case 'EEG flex2flex epochs - Time Domain - Channels'
                disp(['EEG flex2flex epochs - Time Domain - Channels (with ', statMethod, ')']);
                data = EEG_epoched_flextoflex;
                EEG_flextoflex_epochs_Time_Domain_channels(data, statMethod);
            case 'EEG flex2flex epochs - Time Domain - ICs'
                disp(['EEG flex2flex epochs - Time Domain - ICs (with ', statMethod, ')']);
                data = EEG_epoched_flextoflex;
                EEG_flextoflex_epochs_Time_Domain_sources(data, statMethod);
            case 'EEG flex2flex epochs - Frequency Domain - Channels'
                disp(['EEG flex2flex epochs - Frequency Domain - Channels (with ', statMethod, ')']);
                data = EEG_epoched_flextoflex;
                EEG_flextoflex_epochs_Frequency_Domain_channels(data, statMethod);
            case 'EEG flex2flex epochs - Frequency Domain - ICs'
                disp(['EEG flex2flex epochs - Frequency Domain - ICs (with ', statMethod, ')']);
                data = EEG_epoched_flextoflex;
                EEG_flextoflex_epochs_Frequency_Domain_sources(data, statMethod);

            otherwise
                disp('Unknown option');
        end
        
    end
end