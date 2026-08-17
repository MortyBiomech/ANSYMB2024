% ========================================================================
% FUNCTION: calculate_segmentation_params
% ------------------------------------------------------------------------
function [seg_nFreq, seg_nTime, freq_groups, time_segments] = ...
    calculate_segmentation_params(orig_nFreq, orig_nTime, freq_group_size, time_segment_percent)
    
    % Calculate number of segments
    seg_nTime = round(1 / time_segment_percent);  % 5% = 20 segments
    seg_nFreq = ceil(orig_nFreq / freq_group_size);  % Group every 3 frequencies
    
    % Create time segments (5% each)
    time_edges = round(linspace(1, orig_nTime+1, seg_nTime+1));
    time_segments = cell(seg_nTime, 1);
    for t = 1:seg_nTime
        time_segments{t} = time_edges(t):(time_edges(t+1)-1);
    end
    
    % Create frequency groups (every 3 frequencies)
    freq_groups = cell(seg_nFreq, 1);
    for f = 1:seg_nFreq
        start_freq = (f-1) * freq_group_size + 1;
        end_freq = min(f * freq_group_size, orig_nFreq);
        freq_groups{f} = start_freq:end_freq;
    end
end