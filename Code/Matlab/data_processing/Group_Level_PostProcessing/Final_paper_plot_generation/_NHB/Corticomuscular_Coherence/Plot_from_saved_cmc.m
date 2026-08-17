%% PLOT_FROM_SAVED_CMC
%  Use this script to generate plots from the already-saved sub-18_CMC.mat
%  without recomputing anything. Loads one pair at a time to avoid
%  loading the full 38GB into memory simultaneously.
%
%  Author: Morteza Khosrotabar, TU Darmstadt, 2025

% clear; clc;

%% ---- settings ----------------------------------------------------------
% mat_file    = ['D:\path\to\sub-18_CMC.mat'];   % ← update path
figSavePath = [savePath, '\sub-18']; % ← update path

params.event_labels = {'FlxS', 'FlxE/ExtS', 'ExtE'};
params.freqs        = [3 60];

% which pairs to plot (leave empty {} to plot all)
cluster_filter = {};   % e.g. {'Left_Prim_Motor'} to plot only one
muscle_filter  = {};   % e.g. {'VL','RF'} to plot only some muscles

% %% ---- load struct (metadata only first) ---------------------------------
% fprintf('Loading CMC struct from disk...\n');
% S = load(mat_file, 'cmc');   % loads the full struct — unavoidable
% cmc = S.cmc;
% clear S;
% fprintf('Loaded. Size: [%d x %d]\n', size(cmc,1), size(cmc,2));

%% ---- loop and plot one pair at a time ----------------------------------
[nClusters, nMuscles] = size(cmc);

for ci = 1:nClusters
    for mi = 1:nMuscles

        s = cmc(ci, mi);

        % apply filters
        if ~isempty(cluster_filter) && ...
                ~any(strcmpi(cluster_filter, s.cluster_name)), continue, end
        if ~isempty(muscle_filter) && ...
                ~any(strcmpi(muscle_filter, s.muscle_name)),   continue, end

        % skip pairs that were not computed
        if ~s.computed
            fprintf('Skipping %s — %s (not computed)\n', ...
                s.cluster_name, s.muscle_name);
            continue
        end

        fprintf('\nPlotting: %s — %s\n', s.cluster_name, s.muscle_name);

        pairSavePath = fullfile(figSavePath, ...
            sprintf('%s_%s', s.cluster_name, s.muscle_name));

        plot_CMC_pair(                                                  ...
            s.coherres,           s.crossspec_warped,                  ...
            s.alltfX,             s.alltfX_pow_warped,                 ...
            s.alltfY,             s.alltfY_pow_warped,                 ...
            s.timesout,           s.freqsout,                          ...
            s.warpingvalues,      s.condition_vector,                  ...
            s.cluster_name,       s.muscle_name,                       ...
            s.subject_id,         params,                              ...
            pairSavePath);

        close all;  % clear figures between pairs

    end
end

fprintf('\nAll done. Figures saved to: %s\n', figSavePath);