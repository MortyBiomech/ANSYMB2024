function plot_EEG_features_heatmap_per_pressure(RMS_Freq_Region, ...
    heatmap_diff_1_2, heatmap_diff_1_3, heatmap_diff_2_3, ...
    significance_1_2, significance_1_3, significance_2_3, subject)

   
    %% Plot heatmaps
    % Define extreme colors
    synch_color = [214, 40, 40]/255; % Replace with your desired RGB value for the maximum
    desynch_color = [58, 134, 255]/255;  % Replace with your desired RGB value for the minimum
    
    % Define the number of colors for the colormap
    num_colors = 256;
    
    % Create a gradient for the negative side (red to white)
    neg_colors = [linspace(desynch_color(1), 1, num_colors/2)', ...
                  linspace(desynch_color(2), 1, num_colors/2)', ...
                  linspace(desynch_color(3), 1, num_colors/2)'];
    
    % Create a gradient for the positive side (white to blue)
    pos_colors = [linspace(1, synch_color(1), num_colors/2)', ...
                  linspace(1, synch_color(2), num_colors/2)', ...
                  linspace(1, synch_color(3), num_colors/2)'];
    
    % Combine the gradients to create the full colormap
    custom_cmap = [neg_colors; pos_colors];
    
    % Calculate the global range for all heatmaps
    global_min = min([heatmap_diff_1_2(:); heatmap_diff_1_3(:); heatmap_diff_2_3(:)]);
    global_max = max([heatmap_diff_1_2(:); heatmap_diff_1_3(:); heatmap_diff_2_3(:)]);
    
    % Set symmetric limits around zero for consistent colormap
    global_limit = max(abs([global_min, global_max]));
    
    
    
    % Labels for brain regions and frequency bands
    % brain_regions_complete = {'Left PreMot SuppMot', 'Left Paracentral Lobule', ...
    %     'Left Dorsal ACC', 'Left VisMotor', 'Left PrimVisual', ...
    %     'Right PreMot SuppMot', 'Right VisMotor', 'Right PrimVisual'}; 
    brain_regions = {'LPS', 'LPL', 'F', 'LP', 'LO', 'RPS', 'RO', 'RP'}; 
    freq_bands = {'\delta', '\theta', '\alpha', '\beta', '\gamma'};
    
    
    % Plot for Condition 1 vs. 2
    figure;
    h = imagesc(flipud(heatmap_diff_1_2)); % Use imagesc for axes-based plotting
    colormap(custom_cmap); % Set colormap
    clim([-global_limit, global_limit]); % Set symmetric color limits based on data

    hold on

    % Set AlphaData to hide NaN regions
    h.AlphaData = ~isnan(flipud(heatmap_diff_1_2));

    % Overlay crosses where NaNs are located
    [row, col] = find(isnan(flipud(heatmap_diff_1_2)));
    for i = 1:length(row)
        % Coordinates for the cross
        x = col(i);
        y = row(i);
        plot(x, y, 'kx', 'MarkerSize', 45, 'LineWidth', 2); % Big black cross
    end


    % Add colorbar and label
    cb = colorbar;
    cb.Label.String = 'Effect Size'; % Label for color bar
    cb.Label.FontSize = 14; % Set font size for the label
    cb.Ticks = linspace(-global_limit, global_limit, 5); % Generate ticks
    cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), cb.Ticks, 'UniformOutput', false); % Format ticks to 2 decimals
    
    title(['Subject ', num2str(subject),': Pressure 3 bar - Pressure 1 bar'], 'FontSize', 14, 'FontWeight', 'normal');
    xlabel('Brain Regions', 'FontSize', 14);
    ylabel('Frequency Bands', 'FontSize', 14);
    
    % Reverse the y-axis labels and align with flipped data
    yticks(1:length(freq_bands));
    yticklabels(flip(freq_bands)); % Flip the labels to match the flipped data
    
    % Customize axis ticks and labels
    xticks(1:length(brain_regions));
    xticklabels(brain_regions);
    
    % Increase font size of axis labels
    set(gca, 'FontSize', 14);
    
    % Overlay significance markers
    [row, col] = find(flipud(significance_1_2)); % Find significant cells
    for i = 1:length(row)
        text(col(i), row(i), '*', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'Color', 'black', ...
            'FontWeight', 'bold', 'FontSize', 20);
    end
    
    % Make the figure more horizontal
    set(gcf, 'Position', [100, 100, 700, 400]); % Adjust figure size
    
    
    % Plot for Condition 1 vs. 3
    figure;
    h = imagesc(flipud(heatmap_diff_1_3)); % Use imagesc for axes-based plotting
    colormap(custom_cmap); % Set colormap
    clim([-global_limit, global_limit]); % Set symmetric color limits based on data
    
    hold on

    % Set AlphaData to hide NaN regions
    h.AlphaData = ~isnan(flipud(heatmap_diff_1_3));

    % Overlay crosses where NaNs are located
    [row, col] = find(isnan(flipud(heatmap_diff_1_3)));
    for i = 1:length(row)
        % Coordinates for the cross
        x = col(i);
        y = row(i);
        plot(x, y, 'kx', 'MarkerSize', 45, 'LineWidth', 2); % Big black cross
    end

    % Add colorbar and label
    cb = colorbar;
    cb.Label.String = 'Effect Size'; % Label for color bar
    cb.Label.FontSize = 14; % Set font size for the label
    cb.Ticks = linspace(-global_limit, global_limit, 5); % Generate ticks
    cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), cb.Ticks, 'UniformOutput', false); % Format ticks to 2 decimals
    
    title(['Subject ', num2str(subject),': Pressure 6 bar - Pressure 1 bar'], 'FontSize', 14, 'FontWeight', 'normal');
    xlabel('Brain Regions');
    ylabel('Frequency Bands');
    
    % Reverse the y-axis labels and align with flipped data
    yticks(1:length(freq_bands));
    yticklabels(flip(freq_bands)); % Flip the labels to match the flipped data
    
    % Customize axis ticks and labels
    xticks(1:length(brain_regions));
    xticklabels(brain_regions);
    
    % Increase font size of axis labels
    set(gca, 'FontSize', 14);
    
    % Overlay significance markers
    [row, col] = find(flipud(significance_1_3)); % Find significant cells
    for i = 1:length(row)
        text(col(i), row(i), '*', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'Color', 'black', ...
            'FontWeight', 'bold', 'FontSize', 20);
    end
    
    % Make the figure more horizontal
    set(gcf, 'Position', [100, 100, 700, 400]); % Adjust figure size
    
    
    
    % Plot for Condition 2 vs. 3
    figure;
    h = imagesc(flipud(heatmap_diff_2_3)); % Use imagesc for axes-based plotting
    colormap(custom_cmap); % Set colormap
    clim([-global_limit, global_limit]); % Set symmetric color limits based on data
    
    hold on

    % Set AlphaData to hide NaN regions
    h.AlphaData = ~isnan(flipud(heatmap_diff_2_3));

    % Overlay crosses where NaNs are located
    [row, col] = find(isnan(flipud(heatmap_diff_2_3)));
    for i = 1:length(row)
        % Coordinates for the cross
        x = col(i);
        y = row(i);
        plot(x, y, 'kx', 'MarkerSize', 45, 'LineWidth', 2); % Big black cross
    end

    % Add colorbar and label
    cb = colorbar;
    cb.Label.String = 'Effect Size'; % Label for color bar
    cb.Label.FontSize = 14; % Set font size for the label
    cb.Ticks = linspace(-global_limit, global_limit, 5); % Generate ticks
    cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), cb.Ticks, 'UniformOutput', false); % Format ticks to 2 decimals
    
    title(['Subject ', num2str(subject),': Pressure 6 bar - Pressure 3 bar'], 'FontSize', 14, 'FontWeight', 'normal');
    xlabel('Brain Regions')
    ylabel('Frequency Bands');
    
    % Reverse the y-axis labels and align with flipped data
    yticks(1:length(freq_bands));
    yticklabels(flip(freq_bands)); % Flip the labels to match the flipped data
    
    % Customize axis ticks and labels
    xticks(1:length(brain_regions));
    xticklabels(brain_regions);
    
    % Increase font size of axis labels
    set(gca, 'FontSize', 14);
    
    % Overlay significance markers
    [row, col] = find(flipud(significance_2_3)); % Find significant cells
    for i = 1:length(row)
        text(col(i), row(i), '*', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'Color', 'black', ...
            'FontWeight', 'bold', 'FontSize', 20);
    end
    
    % Make the figure more horizontal
    set(gcf, 'Position', [100, 100, 700, 400]); % Adjust figure size
    
    
    
    %% Plot bar graphs representing mean+-std of RMS features
    P1_mean_color = [0 0.4470 0.7410];
    % P1_shade_color = [96, 168, 214]/255;
    P3_mean_color = [0.4660 0.6740 0.1880];
    % P3_shade_color = [170, 204, 126]/255;
    P6_mean_color = [0.6350 0.0780 0.1840];
    % P6_shade_color = [221, 168, 177]/255;

    figure();
    t = tiledlayout(5,1);

    for freq = 1:size(RMS_Freq_Region, 1)
        
        means = zeros(size(RMS_Freq_Region, 2), 3);
        stds = zeros(size(RMS_Freq_Region, 2), 3);
        for region = 1:size(RMS_Freq_Region, 2)

            if isempty(RMS_Freq_Region{freq, region})
                continue
            end
            
            % P1
            temp = RMS_Freq_Region{freq, region}.Condition_ID == categorical(1);
            means(region, 1) = mean(RMS_Freq_Region{freq, region}.RMS_value(temp), 1);
            stds(region, 1) = std(RMS_Freq_Region{freq, region}.RMS_value(temp), 0, 1);
            % P3
            temp = RMS_Freq_Region{freq, region}.Condition_ID == categorical(3);
            means(region, 2) = mean(RMS_Freq_Region{freq, region}.RMS_value(temp), 1);
            stds(region, 2) = std(RMS_Freq_Region{freq, region}.RMS_value(temp), 0, 1);
            % P6
            temp = RMS_Freq_Region{freq, region}.Condition_ID == categorical(6);
            means(region, 3) = mean(RMS_Freq_Region{freq, region}.RMS_value(temp), 1);
            stds(region, 3) = std(RMS_Freq_Region{freq, region}.RMS_value(temp), 0, 1);
        
        end

        ax = nexttile(5 - freq + 1);
        
        b = bar(means, 'grouped', 'EdgeColor', 'none');
        % Set colors (optional for clarity)
        colors = [P1_mean_color; P3_mean_color; P6_mean_color];
        for k = 1:3
            b(k).FaceColor = colors(k, :);
        end

        hold on

        % Add error bars
        % Determine the x positions of the bars
        numGroups = size(means, 1);
        numBars = size(means, 2);
        groupWidth = min(0.8, numBars/(numBars + 1.5));
        
        for i = 1:numBars
            % X positions for each bar in group
            x = (1:numGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*numBars);
            non_zero_values = means(:, i) ~= 0;
            errorbar(x(non_zero_values), means(non_zero_values, i), ...
                zeros(size(stds(non_zero_values, i))), stds(non_zero_values, i), ...
                'k', 'linestyle', 'none', 'LineWidth', 1, 'CapSize', 8);        
                       
        end


        if freq == 1
            brain_regions_new = cell(1, 8);
            for i = 1:length(brain_regions)
                num_ICs = length(unique(RMS_Freq_Region{1, i}.IC_ID));
                brain_regions_new{1, i} = [brain_regions{1, i}, ...
                    '(', num2str(num_ICs), ')'];
            end
            ax.XTickLabel = brain_regions_new;
        else
            xticks([])
        end

        set(gca, 'FontSize', 14)
        ylabel(freq_bands{1, freq}, 'Rotation', 0, 'FontSize', 16,...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
        
        % Adjust y-axis position for alignment
        ax.YLabel.Units = 'normalized';
        ax.YLabel.Position(1) = -0.06;  % Adjust this value for horizontal alignment
        ax.Box = "off";

        if freq == 5
            legend({'P1', 'P3', 'P6'}, ...
                'Orientation', 'horizontal', 'Location', 'northwest', 'Box', 'off')
        end
    
    end
    
    title(t, ['Subject ', num2str(subject), ': Frequency bands RMS features (Mean + STD)'], 'FontSize', 16)
    set(gcf, 'Position', [100, 100, 1200, 700])
    
    

end