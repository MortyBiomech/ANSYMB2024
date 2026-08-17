function [TF_struct, freqs, freqs_baseline] = TF_calculation_per_IC(subject, ic, trial_info, subject_data, EEGLAB_source, EEGLAB_events)


    %% baseline calculations
    if subject >= 10

        % finding the start_move indexes
        ids = cellfun(@(x) strcmp(x, 'SM_Start_Move'), {EEGLAB_events.type});
        latency = {EEGLAB_events(ids).latency};
    
        baselines = cellfun(@(x) EEGLAB_source(x-1000:x-1), latency, 'UniformOutput', false);
    
        [cwt_coeffs, freqs] =  cellfun(@(x) cwt(x, 'amor', 500, 'FrequencyLimits', [1 70], 'VoicesPerOctave', 32), baselines, 'UniformOutput', false);
        power_baseline = cellfun(@(x) abs(x).^2, cwt_coeffs, 'UniformOutput', false);
        power_baseline_mean = cellfun(@(x) mean(x, 2), power_baseline, 'UniformOutput', false);

        freqs_baseline = freqs{1,1};

    else

        % finding the start_move indexes
        ids = cellfun(@(x) strcmp(x, 'SB_Start_Beep'), {EEGLAB_events.type});
        latency = {EEGLAB_events(ids).latency};

        baselines = cellfun(@(x) EEGLAB_source(x-1000:x-1), latency, 'UniformOutput', false);
    
        [cwt_coeffs, freqs] =  cellfun(@(x) cwt(x, 'amor', 500, 'FrequencyLimits', [1 70], 'VoicesPerOctave', 32), baselines, 'UniformOutput', false);
        power_baseline = cellfun(@(x) abs(x).^2, cwt_coeffs, 'UniformOutput', false);
        power_baseline_mean = cellfun(@(x) mean(x, 2), power_baseline, 'UniformOutput', false);

        freqs_baseline = freqs{1,1};

    end
    

    
    % %% check baseline power
    % figure(); tiledlayout(4, 2); 
    % nexttile()
    % for p = 1:length(power_baseline_mean)
    %     c = [0.5, 0.5, 0.5];
    %     plot(freqs{1,1}, power_baseline_mean{1, p}, 'Color', c)
    %     hold on
    % end
    % plot(freqs{1,1}, mean(cat(2, power_baseline_mean{:}), 2), 'Color', 'k', 'LineWidth', 2)
    % title('Baseline: 4 second before trial')
    % 
    % nexttile()
    % fill([freqs{1,1}; flipud(freqs{1,1})], ...
    %     [mean(cat(2, power_baseline_mean{:}), 2) + std(cat(2, power_baseline_mean{:}), 0, 2); ...
    %     flipud(mean(cat(2, power_baseline_mean{:}), 2))], ...
    %     [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none')
    % hold on
    % plot(freqs{1,1}, mean(cat(2, power_baseline_mean{:}), 2), 'Color', 'k', 'LineWidth', 2)
    % 
    % 
    % figure(); tiledlayout(5, 1);
    % for p = 1:length(power_baseline_mean)
    %     if p < 30
    %         nexttile(1); hold on
    %         plot(freqs{1,1}, power_baseline_mean{i, p}, 'Color', 'r')
    %     elseif 31 < p && p < 60
    %         nexttile(2); hold on
    %         plot(freqs{1,1}, power_baseline_mean{i, p}, 'Color', 'g')
    %     elseif 60 < p && p < 90
    %         nexttile(3); hold on
    %         plot(freqs{1,1}, power_baseline_mean{i, p}, 'Color', 'b')
    %     elseif 90 < p && p < 120
    %         nexttile(4); hold on
    %         plot(freqs{1,1}, power_baseline_mean{i, p}, 'Color', 'c')
    %     elseif 120 < p && p <= 150
    %         nexttile(5); hold on
    %         plot(freqs{1,1}, power_baseline_mean{i, p}, 'Color', 'k')
    %     end
    % 
    % end
    
  



    %%
    % sub_id = sprintf('sub%d', subject);
    % trial_info = subjects_trialsInfo.(sub_id);
    % subject_data = subjects_data.(sub_id);

    % Get experimental trials
    if subject < 10
        isExperiment = 1:length(trial_info);
    else
        isExperiment = find(cellfun(@(x) strcmp(x.General.Description, 'Experiment'), trial_info));
    end

    
    % Extract pressure levels
    pressures = cellfun(@(x) x.General.Pressure, trial_info, 'UniformOutput', false);
    pressures = cell2mat(pressures);
    P1 = intersect(find(pressures == 1), isExperiment);
    P3 = intersect(find(pressures == 3), isExperiment);
    P6 = intersect(find(pressures == 6), isExperiment);
    
    
    
    % Initialize TF storage
    TF_P1 = cell(3, length(P1));
    TF_P3 = cell(3, length(P3));
    TF_P6 = cell(3, length(P6));
    P1_indx = 1; P3_indx = 1; P6_indx = 1;



    %% finding the outliers and ignore the trials which contain them
    % max_power = [];
    max_power_P1 = [];
    max_power_P3 = [];
    max_power_P6 = [];
    for trial = 1:length(subject_data)

        if ~ismember(trial, isExperiment), continue; end
        % disp(trial)

        % Access source signals and timestamps
        if subject < 10
            sources = subject_data{1, trial}.EEG_stream.Preprocessed.Time_Domain.Sources.Not_Length_Normalized;
            times = subject_data{1, trial}.EEG_stream.Preprocessed.Time_Domain.Times;
        else
            sources = subject_data{1, trial}.EEG_stream.Preprocessed.Sources;
            times = subject_data{1, trial}.EEG_stream.Preprocessed.Times;
        end

        signal = [];
        total_time = [];
        time_sizes = zeros(1, length(sources));

        for i = 1:length(sources)
            signal = cat(2, signal, sources{1, i}(ic, :));  % component ic
            total_time = cat(2, total_time, times{1, i}(1, :));
            time_sizes(i) = length(times{1, i});
        end

        [cwt_coeffs, freqs] = cwt(signal, 'amor', 500, 'FrequencyLimits', [1 70], 'VoicesPerOctave', 32);
        power = abs(cwt_coeffs).^2;

        % max_power = cat(2, max_power, max(max(power)));
        % Sort into pressure groups
        if ismember(trial, P1)
            max_power_P1 = cat(2, max_power_P1, [max(max(power)); trial]);
        elseif ismember(trial, P3)
            max_power_P3 = cat(2, max_power_P3, [max(max(power)); trial]);
        elseif ismember(trial, P6)
            max_power_P6 = cat(2, max_power_P6, [max(max(power)); trial]);
        end

    end

    factor = 3;
    % outliers = isoutlier(max_power, "median", "ThresholdFactor", factor);
    outliers_P1 = isoutlier(max_power_P1(1, :), "median", "ThresholdFactor", factor);
    outliers_P3 = isoutlier(max_power_P3(1, :), "median", "ThresholdFactor", factor);
    outliers_P6 = isoutlier(max_power_P6(1, :), "median", "ThresholdFactor", factor);


    P1 = P1(~outliers_P1);
    P3 = P3(~outliers_P3);
    P6 = P6(~outliers_P6);



    % % optional to check the process. plot TF max power trials
    % helper_function_plot_MaxPowerOutliers(subject, ic, max_power, ...
    %     max_power_P1, max_power_P3, max_power_P6)

    



    %% finding the global average and std over all epochs from all conditoins
    total_power = [];
    total_power_P1 = [];
    total_power_P3 = [];
    total_power_P6 = [];
    for trial = 1:length(subject_data)
    
        if ~ismember(trial, isExperiment), continue; end
        % disp(trial)
    
        % Access source signals and timestamps
        if subject < 10
            sources = subject_data{1, trial}.EEG_stream.Preprocessed.Time_Domain.Sources.Not_Length_Normalized;
            times = subject_data{1, trial}.EEG_stream.Preprocessed.Time_Domain.Times;
        else
            sources = subject_data{1, trial}.EEG_stream.Preprocessed.Sources;
            times = subject_data{1, trial}.EEG_stream.Preprocessed.Times;
        end
    
        signal = [];
        total_time = [];
        time_sizes = zeros(1, length(sources));
        
        for i = 1:length(sources)
            signal = cat(2, signal, sources{1, i}(ic, :));  % component ic
            total_time = cat(2, total_time, times{1, i}(1, :));
            time_sizes(i) = length(times{1, i});
        end
    
        [cwt_coeffs, freqs] = cwt(signal, 'amor', 500, 'FrequencyLimits', [1 70], 'VoicesPerOctave', 32);
        power = abs(cwt_coeffs).^2;

        
        % Sort into pressure groups
        if ismember(trial, P1)
            total_power_P1 = cat(2, total_power_P1, power);
        elseif ismember(trial, P3)
            total_power_P3 = cat(2, total_power_P3, power);
        elseif ismember(trial, P6)
            total_power_P6 = cat(2, total_power_P6, power);
        end

        total_power = cat(2, total_power, power);
        
    end

    grand_mean = mean(total_power, 2);
    % grand_std  = std(total_power, 0, 2);

    grand_mean_P1 = mean(total_power_P1, 2);
    % grand_std_P1  = std(total_power_P1, 0, 2);
    grand_mean_P3 = mean(total_power_P3, 2);
    % grand_std_P3  = std(total_power_P3, 0, 2);
    grand_mean_P6 = mean(total_power_P6, 2);
    % grand_std_P6  = std(total_power_P6, 0, 2);


    %% plot grand means
    figure()
    tiledlayout(2, 1, "TileSpacing", "compact", "Padding", "compact")
    nexttile
    plot(freqs, 10*log10(grand_mean_P1))
    hold on
    plot(freqs, 10*log10(grand_mean_P3))
    plot(freqs, 10*log10(grand_mean_P6))
    plot(freqs, 10*log10(grand_mean), 'k', 'LineWidth', 2, 'LineStyle', '--')
    ylabel('Power (\muV^2)')
    xlabel('Frequency (Hz)')
    legend({'P1', 'P3', 'P6', 'Total'})
    title(['grand mean ', num2str(subject), ' ', num2str(ic)])
    
    % nexttile
    % plot(freqs, 10*log10(grand_std_P1))
    % hold on
    % plot(freqs, 10*log10(grand_std_P3))
    % plot(freqs, 10*log10(grand_std_P6))
    % title(['grand std ', num2str(subject), ' ', num2str(ic)])


    %% main loop over trials
    if exist('power_baseline_mean', 'var')
        final_baseline_average_total = mean(cat(2, power_baseline_mean{:}), 2);
        % final_baseline_std_total     = std(cat(2, power_baseline_mean{:}), 0, 2);
        % final_baseline_average_P1 = mean(cat(2, power_baseline_mean{1, P1}), 2);
        % final_baseline_average_P3 = mean(cat(2, power_baseline_mean{1, P3}), 2);
        % final_baseline_average_P6 = mean(cat(2, power_baseline_mean{1, P6}), 2);
    else
        freqs_baseline = freqs;
    end
    for trial = 1:length(subject_data)

        if ~ismember(trial, isExperiment), continue; end
        % disp(trial)

        % Access source signals and timestamps
        if subject < 10
            sources = subject_data{1, trial}.EEG_stream.Preprocessed.Time_Domain.Sources.Not_Length_Normalized;
            times = subject_data{1, trial}.EEG_stream.Preprocessed.Time_Domain.Times;
        else
            sources = subject_data{1, trial}.EEG_stream.Preprocessed.Sources;
            times = subject_data{1, trial}.EEG_stream.Preprocessed.Times;
        end

        signal = [];
        total_time = [];
        time_sizes = zeros(1, length(sources));

        for i = 1:length(sources)
            signal = cat(2, signal, sources{1, i}(ic, :));  % component ic
            total_time = cat(2, total_time, times{1, i}(1, :));
            time_sizes(i) = length(times{1, i});
        end

        [cwt_coeffs, freqs] = cwt(signal, 'amor', 500, 'FrequencyLimits', [1 70], 'VoicesPerOctave', 32);
        power = abs(cwt_coeffs).^2;

        % Linear time warping
        ev = trial_info{1, trial}.Events.EEG_stream.Preprocessed;
        flxS = ev.flextoflex_start_indx;
        ext  = ev.flexion_end_indx;
        flxE = ev.flextoflex_end_indx;
        events = [reshape(flxS, [], 1), ...
                  reshape(ext(1:length(flxS)), [], 1), ...
                  reshape(flxE, [], 1)];

        flex_lens = events(:,2) - events(:,1) + 1;
        ext_lens = events(:,3) - (events(:,2)+1) + 1;
        median_flx = round(median(flex_lens));
        median_ext = round(median(ext_lens));

        c = [0, cumsum(time_sizes)];
        power_epochs = cell(1, length(sources));
        for i = 1:length(sources)
            power_epochs{i} = power(:, c(i)+1:c(i+1));
        end

        % Interpolation and normalization
        power_warped = cell(size(power_epochs));
        for i = 1:length(power_epochs)
            % Extract original time vectors
            original_time = times{1, i};

            flx_len = flex_lens(i);

            % Flexion
            flx_seg = power_epochs{i}(:, 1:flx_len);
            flx_t = original_time(1, 1:flx_len);
            flx_target_t = linspace(flx_t(1), flx_t(end), median_flx);
            flx_interp = interp1(flx_t', flx_seg', flx_target_t, 'linear', 'extrap')';

            % Extension
            ext_seg = power_epochs{i}(:, flx_len+1:end);
            ext_t = original_time(1, flx_len+1:end);
            ext_target_t = linspace(ext_t(1), ext_t(end), median_ext);
            ext_interp = interp1(ext_t', ext_seg', ext_target_t, 'linear', 'extrap')';

            % Combine
            power_warped{i} = [flx_interp, ext_interp];
        end

        
        % Normalization (baseline correction)
        for i = 1:length(power_warped)
            % mu = mean(power_warped{i}, 2);
            % sigma = std(power_warped{i}, 0, 2);
            % sigma(sigma == 0) = 1;  % prevent divide by zero
            % power_warped{i} = (power_warped{i} - mu) ./ sigma;                    % each epoch as its own baseline 
                                                                                    % (meaningless for comparing ERS/ERD across conditions)

            % power_warped{i} = (power_warped{i} - grand_mean) ./ grand_std;        % Global average and std as a baseline
            % power_warped{i} = (power_warped{i} - grand_mean_P1) ./ grand_std_P1;  % condition P1 as a baseline

            if exist('power_baseline_mean', 'var')
                
                % Same Baseline for all pressure conditions: 
                % average of all trials' baseline powers (per frequency) 
                % dB Change
                power_warped{i} = 10*log10(power_warped{i}(1:length(freqs_baseline), :)...
                        ./ final_baseline_average_total);  
                % z-score normalization
                % power_warped{i} = (power_warped{i}(1:length(freqs_baseline), :) - final_baseline_average_total); %...
                %         % ./ final_baseline_std_total;  

                % % each pressure condition has its own baseline (average of
                % % all trials' baselines which have the same pressure
                % if ismember(trial, P1)
                %     power_warped{i} = 10*log10(power_warped{i}(1:length(freqs_baseline), :)...
                %         ./ final_baseline_average_P1);  
                % elseif ismember(trial, P3)
                %     power_warped{i} = 10*log10(power_warped{i}(1:length(freqs_baseline), :)...
                %         ./ final_baseline_average_P3);  
                % elseif ismember(trial, P6)
                %     power_warped{i} = 10*log10(power_warped{i}(1:length(freqs_baseline), :)...
                %         ./ final_baseline_average_P6);  
                % end

            end

        end

        trial_TF = mean(cat(3, power_warped{:}), 3); % average of all epochs in the trial as the trial's TF representative
        % trial_TF = cat(3, power_warped{:}); % exporting the TF data of all epochs as the representative of this trial

        % Sort into pressure groups
        if ismember(trial, P1)
            TF_P1{1, P1_indx} = trial_TF;
            TF_P1{2, P1_indx} = median_flx;
            TF_P1{3, P1_indx} = median_ext;
            P1_indx = P1_indx + 1;
        elseif ismember(trial, P3)
            TF_P3{1, P3_indx} = trial_TF;
            TF_P3{2, P3_indx} = median_flx;
            TF_P3{3, P3_indx} = median_ext;
            P3_indx = P3_indx + 1;
        elseif ismember(trial, P6)
            TF_P6{1, P6_indx} = trial_TF;
            TF_P6{2, P6_indx} = median_flx;
            TF_P6{3, P6_indx} = median_ext;
            P6_indx = P6_indx + 1;
        end
    end

    % Find non-empty columns
    non_empty_cols_P1 = ~all(cellfun(@isempty, TF_P1), 1);
    non_empty_cols_P3 = ~all(cellfun(@isempty, TF_P3), 1);
    non_empty_cols_P6 = ~all(cellfun(@isempty, TF_P6), 1);

    % Keep only those columns
    TF_P1 = TF_P1(:, non_empty_cols_P1);
    TF_P3 = TF_P3(:, non_empty_cols_P3);
    TF_P6 = TF_P6(:, non_empty_cols_P6);

    % Final struct
    TF_struct = struct('P1', {TF_P1}, 'P3', {TF_P3}, 'P6', {TF_P6});
    % TF_struct = [];

end
