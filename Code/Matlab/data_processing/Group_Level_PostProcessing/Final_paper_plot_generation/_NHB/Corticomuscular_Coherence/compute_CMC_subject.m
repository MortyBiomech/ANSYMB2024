% COMPUTE_CMC_SUBJECT - Compute and plot time-resolved corticomuscular
%   coherence (CMC) for one subject, one cluster-muscle pair at a time.
%
%   No large structs are saved to disk. Each cluster-muscle pair is
%   computed, plotted, figures saved, and memory cleared before moving
%   to the next pair. If you want to re-investigate a specific pair,
%   re-run for just that cluster and muscle using the 'cluster_filter'
%   and 'muscle_filter' params fields.
%
% Usage:
%   >> compute_CMC_subject(subject_id, EEG_epoched,                   ...
%                EMG_epochs_CMC, epochs_to_keep, condition_vector_final,...
%                SUBJECTS_ICS, warpingvalues, params, savePath);
%
% Inputs:
%   subject_id             - [integer] subject number e.g. 18
%   EEG_epoched            - EEGLAB EEG structure (epoched, ICA activations
%                            and timewarp structure present)
%   EMG_epochs_CMC         - [nMuscles x nSamples x nEpochs] rectified,
%                            bandpass-filtered EMG aligned to EEG epochs
%                            (epochs_to_keep already applied)
%   epochs_to_keep         - [logical or index vector] good epoch indices
%   condition_vector_final - [nEpochs x 1] condition label per epoch
%   SUBJECTS_ICS           - loaded Subjects_ICs_in_clusters.mat variable
%   warpingvalues          - [1 x nEvents] group median event latencies ms
%   params                 - struct with fields:
%       .srate             - EEG sampling rate Hz
%       .cycles            - wavelet cycles e.g. [3 0.8]
%       .freqs             - [min max] Hz e.g. [3 60]
%       .nfreqs            - number of frequencies e.g. 100
%       .ntimesout         - number of time points e.g. 200
%       .subject_id_offset - offset: subject_id - offset = SUBJECTS_ICS index
%       .motor_clusters    - cell array of cluster names to process
%       .emg_labels        - cell array of EMG channel labels
%       .event_labels      - cell array of event names for plot lines
%                            e.g. {'FlxS','FlxE/ExtS','ExtE'}
%       .cluster_filter    - [optional] cell array, process only these clusters
%       .muscle_filter     - [optional] cell array, process only these muscles
%   savePath               - root folder to save figures
%                            Figures saved as:
%                            [savePath]/CMC_figures/sub-[id]/[cluster]_[muscle]/
%
% Author: Morteza Khosrotabar, TU Darmstadt, 2026
% Based on pipeline by Noelle Jacobsen, University of Florida

function compute_CMC_subject(subject_id, EEG_epoched,                 ...
    EMG_epochs_CMC, epochs_to_keep, condition_vector_final,            ...
    SUBJECTS_ICS, warpingvalues, params, savePath)

fprintf('\n========================================================\n');
fprintf('CMC computation: subject %d\n', subject_id);
fprintf('========================================================\n');

% =========================================================================
% Derive tlimits from EEG_epoched.xmin / xmax (seconds -> ms)
% =========================================================================
tlimits  = [EEG_epoched.xmin, EEG_epoched.xmax] * 1000;
nSamples = EEG_epoched.pnts;
fprintf('  tlimits   : [%.1f, %.1f] ms\n', tlimits(1), tlimits(2));
fprintf('  nSamples  : %d\n', nSamples);

% =========================================================================
% Basic validation
% =========================================================================
nEpochs   = size(EMG_epochs_CMC, 3);
nMuscles  = length(params.emg_labels);
nClusters = length(params.motor_clusters);

if size(EMG_epochs_CMC, 2) ~= nSamples
    error('EMG sample count (%d) does not match EEG pnts (%d).', ...
        size(EMG_epochs_CMC,2), nSamples);
end
if ~isfield(EEG_epoched,'timewarp') || isempty(EEG_epoched.timewarp)
    error('EEG_epoched.timewarp not found.');
end
if ~isfield(params,'subject_id_offset'), params.subject_id_offset = 0; end
if ~isfield(params,'event_labels')
    params.event_labels = arrayfun(@(x) sprintf('E%d',x), ...
        1:length(warpingvalues), 'UniformOutput', false);
end

% optional cluster/muscle filters
if isfield(params,'cluster_filter') && ~isempty(params.cluster_filter)
    cluster_list = params.cluster_filter;
else
    cluster_list = params.motor_clusters;
end
if isfield(params,'muscle_filter') && ~isempty(params.muscle_filter)
    muscle_list = params.muscle_filter;
else
    muscle_list = params.emg_labels;
end

% timewarp latencies for good epochs only
timewarp_latencies = EEG_epoched.timewarp.latencies(epochs_to_keep, :);
fprintf('  Timewarp  : [%d epochs x %d events]\n', ...
    size(timewarp_latencies,1), size(timewarp_latencies,2));
fprintf('  Warpvals  :'); fprintf(' %.1f', warpingvalues); fprintf(' ms\n');
fprintf('  Conditions: %s\n', num2str(unique(condition_vector_final)'));
fprintf('  Clusters  : %s\n', strjoin(cluster_list, ', '));
fprintf('  Muscles   : %s\n', strjoin(muscle_list, ', '));

% =========================================================================
% Main loop — one pair at a time
% =========================================================================
pair_count  = 0;
pair_ok     = 0;
pair_skip   = 0;
pair_fail   = 0;

for ci = 1:nClusters

    clusterName = params.motor_clusters{ci};

    % apply cluster filter
    if ~any(strcmpi(cluster_list, clusterName)), continue, end

    fprintf('\n--- Cluster: %s ---\n', clusterName);

    % find cluster in SUBJECTS_ICS
    cluster_indx = find(strcmpi(SUBJECTS_ICS(:,1), clusterName));
    if isempty(cluster_indx)
        fprintf('  Not found in SUBJECTS_ICS. Skipping.\n');
        pair_skip = pair_skip + nMuscles;
        continue
    end

    clusterData = SUBJECTS_ICS{cluster_indx, 2};
    subject_indx_in_cluster = find( ...
        clusterData.Subjects == subject_id - params.subject_id_offset);

    if isempty(subject_indx_in_cluster)
        fprintf('  Subject %d not in this cluster. Skipping.\n', subject_id);
        pair_skip = pair_skip + nMuscles;
        continue
    end

    ic_idx = clusterData.ICs(subject_indx_in_cluster);
    fprintf('  IC index: %d\n', ic_idx);

    % extract EEG IC activation for good epochs -> [nSamples x nEpochs]
    EEG_IC = squeeze(EEG_epoched.icaact(ic_idx, :, epochs_to_keep));
    x      = reshape(EEG_IC, 1, []);   % [1 x nSamples*nEpochs]

    % ---- muscle loop ----------------------------------------------------
    for mi = 1:nMuscles

        muscleName = params.emg_labels{mi};

        % apply muscle filter
        if ~any(strcmpi(muscle_list, muscleName)), continue, end

        pair_count = pair_count + 1;
        fprintf('  [Pair %d] %s — %s ... ', pair_count, clusterName, muscleName);

        % extract EMG -> [nSamples x nEpochs]
        EMG_muscle = squeeze(EMG_epochs_CMC(mi, :, :));
        y          = reshape(EMG_muscle, 1, []);

        % ---- compute CMC ------------------------------------------------
        try
            [coherres, ~, timesout, freqsout, ~, ~,                    ...
             crossspec_warped, alltfX, alltfX_pow_warped,               ...
             alltfY, alltfY_pow_warped]                                 ...
                = my_newcrossf(                                          ...
                    x, y,                                                ...
                    nSamples,                                            ...
                    tlimits,                                             ...
                    params.srate,                                        ...
                    params.cycles,                                       ...
                    'freqs',      params.freqs,                          ...
                    'nfreqs',     params.nfreqs,                         ...
                    'freqscale',  'log',                                 ...
                    'ntimesout',  params.ntimesout,                      ...
                    'timewarp',   timewarp_latencies,                    ...
                    'timewarpms', warpingvalues,                         ...
                    'type',       'coher',                               ...
                    'alpha',      NaN,                                   ...
                    'baseline',   NaN,                                   ...
                    'plotamp',    'off',                                  ...
                    'plotphase',  'off');

            fprintf('done.\n');
            pair_ok = pair_ok + 1;

            % ---- plot and save figures -----------------------------------
            figSavePath = fullfile(savePath, 'CMC_figures',             ...
                sprintf('sub-%d', subject_id),                           ...
                sprintf('%s_%s', clusterName, muscleName));

            plot_CMC_pair(coherres, crossspec_warped,                   ...
                alltfX, alltfX_pow_warped,                              ...
                alltfY, alltfY_pow_warped,                              ...
                timesout, freqsout,                                     ...
                warpingvalues, condition_vector_final,                  ...
                clusterName, muscleName, subject_id,                    ...
                params, figSavePath);

        catch ME
            fprintf('FAILED: %s\n', ME.message);
            pair_fail = pair_fail + 1;
        end

        % clear large arrays before next pair
        clear coherres crossspec_warped alltfX alltfX_pow_warped ...
              alltfY alltfY_pow_warped timesout freqsout
        close all

    end % muscle loop
end % cluster loop

% =========================================================================
% Final summary
% =========================================================================
fprintf('\n--- Subject %d complete ---\n', subject_id);
fprintf('  Pairs computed : %d\n', pair_ok);
fprintf('  Pairs skipped  : %d (no IC for subject)\n', pair_skip);
fprintf('  Pairs failed   : %d (runtime error)\n', pair_fail);

end