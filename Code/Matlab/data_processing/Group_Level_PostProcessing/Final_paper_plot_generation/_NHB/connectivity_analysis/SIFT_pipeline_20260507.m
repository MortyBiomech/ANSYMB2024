%% =====================================================================
%  SIFT Pipeline — Single Subject (subject 18 for now)
%  ---------------------------------------------------------------------
%  Time-Varying Directed Connectivity Between Cortical Sources
%  Re-epoched around ExtS and ExtE; uses sliding-window AMVAR + dDTF08
%  =====================================================================
%
%  Prerequisites:
%    - EEGLAB started (run `eeglab` once before this script)
%    - SIFT plugin installed
%    - SUBJECTS_ICS variable in workspace
%    - build_subject_ic_map.m on path
%    - Continuous dataset for this subject with AMICA already run
%
%  =====================================================================

clc; clear;

% ---- Add the necessary paths --------------------------------------------
addpath(genpath('D:\Morteza\MyProjects\ANSYMB2024\Code'))
EEGLAB_path = 'D:\Morteza\Toolboxes\EEGLAB\eeglab2026.0.0';
data_path = 'D:\Morteza\MyProjects\ANSYMB2024\data\';
processing_path = 'D:\Morteza\MyProjects\ANSYMB2024\Code\Matlab\data_processing\';
continuous_EEG_path = [data_path, '5_single-subject-EEG-analysis\'];
connectivity_path = [processing_path, 'Group_Level_PostProcessing\', ...
    'Final_paper_plot_generation\_NHB\connectivity_analysis'];



% ---- run EEGLAB ---------------------------------------------------------
this_path = pwd;
cd(EEGLAB_path)
if ~exist('ALLEEG', 'var')
    eeglab
end
cd(this_path)



% ---- Load continuous EEG dataset ----------------------------------------
% thisSubjectID = 18;     
% cd([continuous_EEG_path, 'sub-', num2str(thisSubjectID)])
% fileName = ['sub-', num2str(thisSubjectID), '_cleaned_with_ICA.set'];
% EEG = pop_loadset('filename', fileName);
cd(this_path)



% ---- Load SUBJECTS_ICS: Contains Subjects and ICs in Brain regions ------
load("Subjects_ICs_in_clusters.mat") 


% ---- USER PARAMETERS ----------------------------------------------------
thisSubjectID    = 17;
continuousFile   = ['sub-', num2str(thisSubjectID), '_cleaned_with_ICA.set'];
continuousPath   = [continuous_EEG_path, 'sub-', num2str(thisSubjectID)];                  
outputPath       = [connectivity_path, '\sift_results'];          

epochWindow      = [-0.75 0.75];                     % seconds around each event
freqRange        = 3:0.5:30;                           % Hz
windowStepSec    = 0.02;                             % 10 samples @ 500 Hz
modelOrderRange  = [1 50];
finalModelOrder  = [];                               % set after STEP 2 inspection
verbosity        = 2;

windowSafetyFactor = 1.2;
minCyclesLowFreq   = 1.0;

runStats         = false;     % set true to compute phase-randomization stats
nPermutations    = 200;

% Create output folder if needed
if ~exist(outputPath, 'dir'), mkdir(outputPath); end


%% =====================================================================
%  BUILD IC MAP AND SELECT THIS SUBJECT
%  =====================================================================
subjectMap = build_subject_ic_map(SUBJECTS_ICS);

entry = subjectMap([subjectMap.subjectID] == thisSubjectID);
if isempty(entry)
    error('Subject %d not found in SUBJECTS_ICS.', thisSubjectID);
end

srcIdx_orig = entry.icIdx;          % original IC numbers (in continuous EEG)
srcLabels   = entry.regionNames;

fprintf('\n=== Subject %d: %d ICs ===\n', thisSubjectID, numel(srcIdx_orig));
for i = 1:numel(srcIdx_orig)
    fprintf('   IC%2d -> %s\n', srcIdx_orig(i), srcLabels{i});
end


%% =====================================================================
%  LOAD CONTINUOUS DATASET AND RE-EPOCH AROUND FlxS AND ExtS
%  =====================================================================
fprintf('\nLoading continuous dataset...\n');
EEG_cont = pop_loadset('filename', continuousFile, 'filepath', continuousPath);

if isempty(EEG_cont.icaweights)
    error('Continuous dataset has no ICA decomposition.');
end
if max(srcIdx_orig) > size(EEG_cont.icaweights, 1)
    error('IC index %d exceeds available components (%d).', ...
          max(srcIdx_orig), size(EEG_cont.icaweights, 1));
end

% Confirm events exist
eventTypes = unique({EEG_cont.event.type});
fprintf('Event types in dataset: %s\n', strjoin(eventTypes, ', '));
for ev = {'ExtS', 'ExtE'}
    if ~ismember(ev{1}, eventTypes)
        error('Event "%s" not found in EEG.event.', ev{1});
    end
    nEv = sum(strcmp({EEG_cont.event.type}, ev{1}));
    fprintf('  %s: %d occurrences\n', ev{1}, nEv);
end

% Downsample to 250 Hz (anti-aliased via pop_resample)
newFs = 250;
if EEG_cont.srate > newFs
    fprintf('Downsampling from %d to %d Hz...\n', EEG_cont.srate, newFs);
    EEG_cont = pop_resample(EEG_cont, newFs);
    EEG_cont = eeg_checkset(EEG_cont);
    Fs  = EEG_cont.srate;    % update Fs for window-length computation
    fprintf('New sampling rate: %d Hz, samples/trial: %d\n', EEG_cont.srate, EEG_cont.pnts);
end

% Re-epoch
fprintf('\nRe-epoching...\n');
EEG_FlxS = pop_epoch(EEG_cont, {'FlxS'}, epochWindow, ...
                     'newname',   sprintf('Subj%d around FlxS', thisSubjectID), ...
                     'epochinfo', 'yes');
EEG_FlxS = eeg_checkset(EEG_FlxS);

EEG_ExtS = pop_epoch(EEG_cont, {'ExtS'}, epochWindow, ...
                     'newname',   sprintf('Subj%d around ExtS', thisSubjectID), ...
                     'epochinfo', 'yes');
EEG_ExtS = eeg_checkset(EEG_ExtS);

EEG_ExtE = pop_epoch(EEG_cont, {'ExtE'}, epochWindow, ...
                     'newname',   sprintf('Subj%d around ExtE', thisSubjectID), ...
                     'epochinfo', 'yes');
EEG_ExtE = eeg_checkset(EEG_ExtE);

fprintf('  FlxS-locked: %d trials, %d samples/trial\n', EEG_FlxS.trials, EEG_FlxS.pnts);
fprintf('  ExtS-locked: %d trials, %d samples/trial\n', EEG_ExtS.trials, EEG_ExtS.pnts);
fprintf('  ExtE-locked: %d trials, %d samples/trial\n', EEG_ExtE.trials, EEG_ExtE.pnts);


clear EEG_cont;   % free memory


%% =====================================================================
%  RUN SIFT PIPELINE ON EACH RE-EPOCHED DATASET
%  =====================================================================
datasets = {EEG_FlxS, EEG_ExtS};
labels   = {'FlxS',   'ExtS'};
results  = cell(1, numel(datasets));

% Map full source labels -> short labels for the connectivity plot.
% Both naming variants included (PreMot_SuppMot vs PreMotor_SuppMotor).
shortMap = containers.Map('KeyType','char','ValueType','char');
shortMap('Left_PreMot_SuppMot')     = 'LPS';
shortMap('Left_PreMotor_SuppMotor') = 'LPS';
shortMap('Right_PreMot_SuppMot')    = 'RPS';
shortMap('Right_PreMotor_SuppMotor')= 'RPS';
shortMap('Left_Parieto_Occipital')  = 'LPO';
shortMap('Right_Parieto_Occipital') = 'RPO';
shortMap('Left_Prim_Motor')         = 'LM1';
shortMap('Right_Prim_Motor')        = 'RM1';
shortMap('Left_Dorsal_ACC')         = 'dACC';
shortMap('Prime_Visual')             = 'Visual';

for d = 1:numel(datasets)
    EEG     = datasets{d};
    evLabel = labels{d};

    fprintf('\n=========================================================\n');
    fprintf('  SUBJECT %d  |  EVENT-LOCKED: %s\n', thisSubjectID, evLabel);
    fprintf('=========================================================\n');


    %% --- STEP A: Subset components to sources of interest ---------------
    EEG = pop_subcomp(EEG, srcIdx_orig, 0, 1);    % keep only chosen ICs

    [~, new_srcIdx] = sort(srcIdx_orig);
    srcIdx_orig = srcIdx_orig(new_srcIdx);
    srcLabels = srcLabels(new_srcIdx);
    fprintf('\nIC ordering after subset:\n');
    for i = 1:numel(new_srcIdx)
        fprintf('  New IC %d  <-  original IC %2d  (%s)\n', ...
                i, srcIdx_orig(i), srcLabels{i});
    end

    % After subset: srcIdx is just 1:M
    srcIdx = 1:numel(srcIdx_orig);

    % shorten the labels
    srcLabels_short = srcLabels;            % keep the original intact
    for k = 1:numel(srcLabels)
        lbl = srcLabels{k};
        if isKey(shortMap, lbl)
            srcLabels_short{k} = shortMap(lbl);
        else
            warning('No short label for "%s" (index %d) — left unchanged.', lbl, k);
        end
    end


    %% --- STEP B: Compute per-subject window length ---------------------
    M  = numel(srcIdx);
    N  = EEG.trials;
    Fs = EEG.srate;
    pUpper = modelOrderRange(2);

    W_rule    = ceil(windowSafetyFactor * 10 * M^2 * pUpper / N);
    W_freq    = ceil(minCyclesLowFreq * Fs / min(freqRange));
    % W_samples = max(W_rule, W_freq);
    W_samples = 160; % manually increase the window length
    windowLenSec = W_samples / Fs;

    fprintf('\nWindow length: %d samples (%.3f s)  [rule=%d, freq-floor=%d]\n', ...
            W_samples, windowLenSec, W_rule, W_freq);

    if windowLenSec > 0.5 * (EEG.pnts / Fs)
        warning('Window (%.2fs) > half the epoch. Consider fewer ICs.', windowLenSec);
    end

    
    %% --- STEP 1: PREPROCESSING ----------------------------------------
    fprintf('\n--- Preprocessing ---\n');
    EEG = pop_pre_prepData(EEG, 'nogui', ...
        'VerbosityLevel',   verbosity, ...
        'SignalType',       {'Components'}, ...
        'VariableNames',    srcLabels_short, ...
        'Detrend',          {'verb', 0, 'method', {'linear'}}, ...
        'NormalizeData',    {'verb', 0, 'method', {'time','ensemble'}}, ...
        'BadDataSegments',  [], ...
        'TrialsToKeep',     [], ...
        'BalanceTrials',    false);

    % % % % Quick diagnostic: check if normalization was applied
    % % % % Check if the ERP was actually removed (ensemble normalization)
    % % % % After proper ensemble normalization, the trial-average should be ~zero everywhere
    % % % trialAvg = mean(EEG.CAT.srcdata, 3);   % average across trials (dim 3)
    % % % 
    % % % figure;
    % % % plot(trialAvg');
    % % % title('Trial-averaged source activity (should be ~zero if ensemble-normalized)');
    % % % ylabel('Amplitude');
    % % % xlabel('Samples');
    % % % legend(srcLabels, 'Interpreter', 'none');
    % % % 
    % % % % Also check: is the data in EEG.CAT.srcdata or EEG.icaact?
    % % % % SIFT stores its preprocessed data in EEG.CAT.srcdata. 
    % % % % Verify it exists and has the right dimensions:
    % % % fprintf('EEG.CAT.srcdata size: [%s]\n', num2str(size(EEG.CAT.srcdata)));
    % % % fprintf('Expected: [%d x %d x %d]  (sources x samples x trials)\n', ...
    % % %         numel(srcLabels), EEG.pnts, EEG.trials);


    %% --- STEP 2: MODEL ORDER SELECTION --------------------------------
    fprintf('\n--- Model order selection (inspect figure) ---\n');
    modelOrderRange = [1 50];
    EEG = pop_est_selModelOrder(EEG, 'nogui', ...
        'modelingApproach', {'Segmentation VAR', ...
                             'algorithm',       {'Vieira-Morf'}, ...
                             'winStartIdx',     [], ...
                             'winlen',          windowLenSec, ...
                             'winstep',         windowStepSec, ...
                             'taperfcn',        'rectwin', ...
                             'epochTimeLims',   [], ...
                             'prctWinToSample', 100, ...
                             'normalize',       [], ...
                             'detrend',         {'method', 'linear'}, ...
                             'verb',            verbosity}, ...
        'morderRange',     modelOrderRange, ...
        'downdate',        true, ...
        'runPll',          [], ...
        'icselector',      {'sbc','aic','fpe','hq'}, ...
        'winStartIdx',     [], ...
        'epochTimeLims',   [], ...
        'prctWinToSample', 80, ...
        'plot',            {'conditions',[],'icselector',{'sbc','aic','fpe','hq'}, ...
                            'minimizer',{'min'},'prclim',90}, ...
        'verb',            verbosity);


    % % % % resolving the issue with hlp_getFcnPreambleText:
    % % % % line 50 in hlp_getFcnPreambleText was changed from:
    % % % % hlpText = deblank(hlpText(1:endpoint-1)); to
    % % % % hlpText = deblank(hlpText(1:endpoint(1)-1));
    % % % mvarDir = fullfile(fileparts(which('est_fitMVAR')), 'mvar');
    % % % mvarFcns = dir(fullfile(mvarDir, '*.m'));
    % % % 
    % % % fprintf('Checking %d MVAR algorithm files:\n', numel(mvarFcns));
    % % % for k = 1:numel(mvarFcns)
    % % %     [~, fname, ~] = fileparts(mvarFcns(k).name);   % strip the .m
    % % %     try
    % % %         txt = hlp_getFcnPreambleText(fname);
    % % %         fprintf('  OK : %s  (%d chars of help)\n', fname, length(txt));
    % % %     catch ME
    % % %         fprintf('  BAD: %s  -- %s\n', fname, ME.message);
    % % %     end
    % % % end


    % --- PAUSE for user to inspect and choose order --------------------
    if isempty(finalModelOrder)
        fprintf('\n>>> Inspect the model-order figure for %s.\n', evLabel);
        userOrder = input('>>> Enter chosen model order (or press Enter for 15): ');
        if isempty(userOrder), userOrder = 15; end
    else
        userOrder = finalModelOrder;
    end

    fprintf('Using model order p = %d for %s\n', userOrder, evLabel);


    %% --- STEP 3: FIT THE FINAL VAR MODEL ------------------------------
    fprintf('\n--- Fitting VAR[%d] model ---\n', userOrder);
    EEG = pop_est_fitMVAR(EEG, 'nogui', ...
        'algorithm',       'Vieira-Morf', ...
        'morder',          userOrder, ...
        'winStartIdx',     [], ...
        'winlen',          windowLenSec, ...
        'winstep',         windowStepSec, ...
        'taperfcn',        'rectwin', ...
        'epochTimeLims',   [], ...
        'prctWinToSample', 100, ...
        'verb',            verbosity, ...
        'timer',           false);


    % %% --- STEP 4: VALIDATE THE MODEL -----------------------------------
    % fprintf('\n--- Validating model ---\n');
    % EEG = pop_est_validateMVAR(EEG, 'nogui', ...
    %     'checkWhiteness',    {'alpha', 0.05, ...
    %                           'statcorrection', 'none', ...
    %                           'numAcfLags', 50, ...
    %                           'whitenessCriteria', ...
    %                               {'Ljung-Box','ACF','Box-Pierce','Li-McLeod'}, ...
    %                           'winStartIdx', [], ...
    %                           'prctWinToSample', 100, ...
    %                           'verb', 0}, ...
    %     'checkConsistency',  {'winStartIdx', [], ...
    %                           'prctWinToSample', 100, ...
    %                           'Nr', [], ...
    %                           'donorm', 0, ...
    %                           'nlags', [], ...
    %                           'verb', 0}, ...
    %     'checkStability',    {'winStartIdx', [], ...
    %                           'prctWinToSample', 100, ...
    %                           'verb', 0}, ...
    %     'checkResidualVariance', [], ...
    %     'prctWinToSample',   100, ...
    %     'winStartIdx',       [], ...
    %     'verb',              verbosity, ...
    %     'plot',              true);


    %% --- STEP 5: ESTIMATE TIME-VARYING CONNECTIVITY -------------------
    fprintf('\n--- Estimating connectivity ---\n');

    % fmin   = 3;     % keep >0 (log undefined at 0) and within model resolution
    % fmax   = 60;
    % nFreqs = 80;
    % freqRange = logspace(log10(fmin), log10(fmax), nFreqs);

    EEG = pop_est_mvarConnectivity(EEG, 'nogui', ...
        'connmethods',       {'dDTF08','nPDC','nDTF','S'}, ...
        'absvalsq',          true, ...
        'spectraldecibels',  true, ...
        'freqs',             freqRange, ...
        'verb',              verbosity);


    % %% --- STEP 6: STATISTICS (optional) --------------------------------
    % if 1 %runStats
    %     fprintf('\n--- Computing surrogate statistics (this may take a while) ---\n');
    %     EEG.CAT.PConn = stat_surrogateGen(EEG, 'nogui', ...
    %         'modelingApproach',     EEG.CAT.configs.est_fitMVAR, ...
    %         'connectivityModeling', EEG.CAT.configs.est_mvarConnectivity, ...
    %         'mode',                 {'PhaseRand','nperms',nPermutations}, ...
    %         'autosave',             [], ...
    %         'verb',                 verbosity);
    % 
    %     EEG.CAT.Stats = stat_surrogateStats(EEG, 'nogui', ...
    %         'statTest',    {'Hnull'}, ...
    %         'tail',        'right', ...
    %         'computeci',   false, ...
    %         'alpha',       0.05, ...
    %         'mcorrection', 'fdr', ...
    %         'verb',        verbosity);
    % end


    %% --- STEP 7: VISUALIZATION — Time-Frequency Grid ------------------
    fprintf('\n--- Generating Time-Frequency Grid ---\n');
    pop_vis_TimeFreqGrid(EEG, 'nogui', ...
        'MatrixLayout',  {'Partial', ...
                          'UpperTriangle','dDTF08', ...
                          'LowerTriangle','dDTF08', ...
                          'Diagonal','S'}, ...
        'clim',          99.7, ...
        'timeRange',     [], ...
        'freqValues',    freqRange, ...
        'pcontour',      [], ...
        'plotCondDiff',  false, ... 'thresholding',  {'Simple', 'prcthresh', [97.5 3]}, ...
        'baseline',      [], ...
        'fighandles',    [], ...
        'smooth',        false, ...
        'xord',          [], ...
        'yord',          [], ...
        'plotorder',     [], ... 
        'topoplot',      'dipole', ...'topoplot', ... 'dipplot',       {'mri'}, ...
        'titleString',   sprintf('Subj %d: dDTF08 around %s', thisSubjectID, evLabel), ...
        'titleFontSize', 12, ...
        'axesFontSize',  10, ...
        'textColor',     [1 1 1], ...
        'backgroundColor',[0 0 0]);


    %% --- STEP 8: SAVE -------------------------------------------------
    saveName = sprintf('subj%02d_%s_SIFT.set', thisSubjectID, evLabel);
    EEG = pop_saveset(EEG, 'filename', saveName, 'filepath', outputPath);
    fprintf('Saved: %s\n', fullfile(outputPath, saveName));

    results{d} = EEG;
end

fprintf('\n=========================================================\n');
fprintf('  Pipeline complete for subject %d\n', thisSubjectID);
fprintf('  Results saved in: %s\n', outputPath);
fprintf('=========================================================\n');

% Keep results in workspace for follow-up analysis / comparison
EEG_ExtS_result = results{1};
EEG_ExtE_result = results{2};