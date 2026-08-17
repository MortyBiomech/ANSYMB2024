%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Regardless of Brain Region Clustering  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
clear

%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); 

data_processing_path = [main_project_folder, '\Code\Matlab\data_processing'];
data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';

EMG_features_path = [data_processing_path, ...
    '\EMG_processing\structured_EMG_data'];
epoched_data_path = [data_path, '6_Trials_Info_and_Epoched_data\'];
grouplevel_postprocessing_path = [data_processing_path, ...
    '\Group_Level_PostProcessing'];
detailed_per_subject_results = [data_processing_path, ...
    '\Detailed_per_subject_processing'];
cleaned_EEG_data_path = [data_path, '5_single-subject-EEG-analysis\'];


%% Regardless of brain region clustering
% ROIs features containing PSD integrals (power) was created before. 
% Look at this path: 
% C:\Morteza\MyProjects\ANSYMB2024\Code\Matlab\data_processing\
%    Group_Level_PostProcessing\Classification_analysis\
%    regardless_of_GroupLevel_clustering
% This file -> main_classification_regardless_of_clustering_ROIs_creation.m


% Load ROIs features
classisfication_path = [data_path, '8_classification\'];
ROIs_features_path = [classisfication_path, ...
    'ROIs_features\regardless_of_clustering'];
fileName = ['ROIs_0_FlextoFlex_per_trial_pressure_Based_Flx_Ext', ...
    '_separate_NoBrainClustering_PSD_integral.mat'];
ROIs = load(fullfile(ROIs_features_path, fileName));
name = fieldnames(ROIs);
ROIs = ROIs.(name{1});



%% Main Loop

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


k = 1;
neg_colors2 = neg_colors + k*(neg_colors - neg_colors(1, :)).*(neg_colors - 1); 
x = flipud(pos_colors);
pos_colors2 = flipud(x + k*(x - x(1, :)).*(x - 1));
custom_cmap2 = [neg_colors2; pos_colors2];



P1_mean_color = [0 0.4470 0.7410];
P3_mean_color = [0.4660 0.6740 0.1880];
P6_mean_color = [0.6350 0.0780 0.1840];
colors = [P1_mean_color; P3_mean_color; P6_mean_color];
freqBands = {'$\gamma$ (30-50 Hz)', '$\beta$ (13-30 Hz)', ...
    '$\alpha$ (8-13 Hz)', '$\theta$ (4-8 Hz)', ...
    '$\delta$ (2-4 Hz)'};       


eeglabpath = which('eeglab.m');
pathtmp = fileparts(eeglabpath);
dipfits = dir(fullfile(pathtmp, 'plugins', 'dipfit*'));
[~, dipfit_order] = sort(cellfun(@(c) str2double(c(7:end)), {dipfits.name}), 'descend');

dipfit_folder = fullfile(pathtmp, 'plugins', dipfits(1).name);
meshdatapath = fullfile(dipfit_folder, 'standard_BEM', 'standard_vol.mat');
mripath = fullfile(dipfit_folder, 'standard_BEM', 'standard_mri.mat');
 



%% Main Loop over subjects
subject_list = 5:18; %9:18; %5:18;

for i = 1:length(subject_list)

    disp(['Subject ', num2str(subject_list(i))])

    data = load([epoched_data_path, 'sub-', num2str(subject_list(i)), ...
        '\Epochs_FlextoFlex_based.mat']);
    name = fieldnames(data);
    data = data.(name{1});

    Trials_Info = load([epoched_data_path, 'sub-', num2str(subject_list(i)), ...
        '\Trials_Info.mat']);
    name = fieldnames(Trials_Info);
    Trials_Info = Trials_Info.(name{1});
    

    %% Calculate the PSD of Flexion and Extensoin per pressure
    disp('Calculating Power Spectral Density ...')
    [freq_data, frequencies] = calculate_PSD_Flx_Ext_separate(data, Trials_Info, subject_list(i));
    

    %% Store Flexion and Extension PSD per pressure
    condition_indices = condition_indices_identifier(Trials_Info, subject_list(i));
        
    P1_Flx_PSD = [];
    P3_Flx_PSD = [];
    P6_Flx_PSD = [];

    P1_Ext_PSD = [];
    P3_Ext_PSD = [];
    P6_Ext_PSD = [];

    for j = 1:length(freq_data)
        flx_psd = freq_data{1, j}.EEG_stream.Preprocessed.Freq_Domain.Sources.Flexion; 
        ext_psd = freq_data{1, j}.EEG_stream.Preprocessed.Freq_Domain.Sources.Extension;

        if ismember(j, condition_indices.P1)
            P1_Flx_PSD = cat(3, P1_Flx_PSD, flx_psd);
            P1_Ext_PSD = cat(3, P1_Ext_PSD, ext_psd);
        elseif ismember(j, condition_indices.P3)
            P3_Flx_PSD = cat(3, P3_Flx_PSD, flx_psd);
            P3_Ext_PSD = cat(3, P3_Ext_PSD, ext_psd);
        elseif ismember(j, condition_indices.P6)
            P6_Flx_PSD = cat(3, P6_Flx_PSD, flx_psd);
            P6_Ext_PSD = cat(3, P6_Ext_PSD, ext_psd);
        end

    end

    

    %% create flexion and extension features
    features = ROIs.(['TotalBrain_sub_', ...
        num2str(subject_list(i))]);
    ICs = cell2mat(features(:, 2));
    
    flexion_features_P1 = cellfun(@(x) x.P1(:, 1:5), features(:, 3), ...
        'UniformOutput', false);
    flexion_features_P3 = cellfun(@(x) x.P3(:, 1:5), features(:, 3), ...
        'UniformOutput', false);
    flexion_features_P6 = cellfun(@(x) x.P6(:, 1:5), features(:, 3), ...
        'UniformOutput', false);
    extension_features_P1 = cellfun(@(x) x.P1(:, 6:10), features(:, 3), ...
        'UniformOutput', false);
    extension_features_P3 = cellfun(@(x) x.P3(:, 6:10), features(:, 3), ...
        'UniformOutput', false);
    extension_features_P6 = cellfun(@(x) x.P6(:, 6:10), features(:, 3), ...
        'UniformOutput', false);


    %% Load EEG data in EEGLAB
    if ~exist('ALLEEG', 'var')
        eeglab
    end
    % Load the eeg file
    fileName = ['sub-', num2str(subject_list(i)), '_cleaned_with_ICA.set']; 
    filePath = [cleaned_EEG_data_path, 'sub-', num2str(subject_list(i))];
    EEG = pop_loadset('filename', fileName, 'filepath', filePath);


    %% loop over ICs
    % tic
    for j = 1:numel(ICs) 
        
        disp(['Plotting IC ', num2str(ICs(j)), ' detailed frequency/time content ...'])

        %% make the main figure
        P1_flx = fliplr(flexion_features_P1{j, 1});
        P3_flx = fliplr(flexion_features_P3{j, 1});
        P6_flx = fliplr(flexion_features_P6{j, 1});
        P1_ext = fliplr(extension_features_P1{j, 1});
        P3_ext = fliplr(extension_features_P3{j, 1});
        P6_ext = fliplr(extension_features_P6{j, 1});


        
        Main_full_figure = figure('WindowState', 'maximized');
        main_tile = tiledlayout(13, 22, "Padding", "compact");
        title(main_tile, ['Subject ', num2str(subject_list(i)), ...
            ' - IC ', num2str(ICs(j))], 'FontSize', 14)

        % monitors = get(0, 'MonitorPositions');  % [left bottom width height] for each monitor
        % nMonitors = size(monitors, 1);
        % monitorID = 2;
        % pos = monitors(monitorID, :);
        % set(Main_full_figure, 'Units', 'pixels', 'Position', pos);

        % movegui(Main_full_figure, 'onscreen');  % ensure it's fully visible
        % jFig = get(handle(Main_full_figure), 'JavaFrame');
        % jFig.setMaximized(true);
        


        %% Left column: Flexion
        tl_flexion_extension = tiledlayout(main_tile, 10, 4, "TileSpacing", "compact", "Padding", "loose");
        tl_flexion_extension.Layout.Tile = 1; % Explicitly set the position in main layout
        tl_flexion_extension.Layout.TileSpan = [10, 4];
        for k = 1:5
            ax = nexttile(tl_flexion_extension, k*8-7, [2 2]); hold on; 
            means = [mean(P1_flx(:,k), 1), mean(P3_flx(:,k), 1), mean(P6_flx(:,k), 1)];
            stds = [std(P1_flx(:,k), 0, 1), std(P3_flx(:,k), 0, 1), std(P6_flx(:,k), 0, 1)];
            b = barh(means, 'grouped', 'EdgeColor', 'none');
            b.FaceColor = 'flat';
            b.CData = colors;
            
            for e = 1:length(means)
                x = means(e);
                y = e;
                dx = stds(e); % horizontal error
                plot([x-dx, x+dx], [y, y], 'k', 'LineWidth', 1); % horizontal line
                plot([x-dx, x-dx], [y-0.2, y+0.2], 'k', 'LineWidth', 1); % left cap
                plot([x+dx, x+dx], [y-0.2, y+0.2], 'k', 'LineWidth', 1); % right cap
            end
            
            ax.YTick = [];
            ax.FontSize = 10;
            ax.TickLabelInterpreter = "latex";
            % ax.XAxisLocation = "top";
            ylabel(ax, freqBands{1, k});
            ax.YLabel.Interpreter = 'latex';
            ax.YLabel.FontSize = 10;

            if k == 5
                xlabel('Power [$\mu V^2$]', 'Interpreter', ...
                    'latex', 'FontWeight', 'normal','FontSize', 10)
            end
            
            if k == 1
                title('Flexion', 'Interpreter', 'latex', ...
                    'FontSize', 12, 'FontWeight', 'normal');
            end

        end


        % right column: Extension
        for k = 1:5
            ax = nexttile(tl_flexion_extension, k*8-5, [2 2]); hold on; 
            means = [mean(P1_ext(:,k), 1), mean(P3_ext(:,k), 1), mean(P6_ext(:,k), 1)];
            stds = [std(P1_ext(:,k), 0, 1), std(P3_ext(:,k), 0, 1), std(P6_ext(:,k), 0, 1)];
            b = barh(means, 'grouped', 'EdgeColor', 'none');
            b.FaceColor = 'flat';
            b.CData = colors;

            for e = 1:length(means)
                x = means(e);
                y = e;
                dx = stds(e); % horizontal error
                plot([x-dx, x+dx], [y, y], 'k', 'LineWidth', 1); % horizontal line
                plot([x-dx, x-dx], [y-0.2, y+0.2], 'k', 'LineWidth', 1); % left cap
                plot([x+dx, x+dx], [y-0.2, y+0.2], 'k', 'LineWidth', 1); % right cap
            end

            ax.YTick = [];
            ax.FontSize = 10;
            ax.TickLabelInterpreter = "latex";
            % ax.XAxisLocation = "top";

            if k == 5
                xlabel('Power [$\mu V^2$]', 'Interpreter', ...
                    'latex', 'FontWeight', 'normal','FontSize', 10)
            end

            if k == 1
                title('Extension', 'Interpreter', 'latex', ...
                    'FontSize', 12, 'FontWeight', 'normal');
            end

        end

        drawnow

        % %% Create dummy bars (invisible) for the legend
        % % Create axes for dummy bars positioned outside visible range
        % ax_leg = axes('Position',[0 0 0.01 0.01],'Visible','off');
        % hold(ax_leg, 'on');
        % h(1) = bar(nan, 'FaceColor',colors(1, :),'EdgeColor','none');
        % h(2) = bar(nan, 'FaceColor',colors(2, :),'EdgeColor','none');
        % h(3) = bar(nan, 'FaceColor',colors(3, :),'EdgeColor','none');
        % hold(ax_leg, 'off');
        % % Legend at bottom (south) of the figure
        % lg = legend(h, {'P1','P3','P6'},...
        %             'Orientation','horizontal','Interpreter','latex', ...
        %             'FontSize',14, 'Box', 'off');
        % % Precisely position legend at bottom of the entire figure
        % lg.Position = [0.5 0.01 0 0] + [-(lg.Position(3)/2), 0, lg.Position(3), lg.Position(4)];
        % 
        % % title
        % title(main_tl, ['Subject ', num2str(subject_list(i)), ' - IC ', num2str(ICs(j))], 'FontSize', 22)



        %% Plot the Time-Frequency figures for each subject (pressure based)
        Subject_ICs_in_ROIs_format = struct(['Subject', num2str(subject_list(i))], ...
            []);
        Subject_ICs_in_ROIs_format.(['Subject', num2str(subject_list(i))]) = ...
            ROIs.(['TotalBrain_sub_', num2str(subject_list(i))]);

        plot_Time_Frequency_per_pressure_IC_based(Subject_ICs_in_ROIs_format, ...
            subject_list(i), epoched_data_path, ICs(j), main_tile)
        

        drawnow



        %% Total IC power spectral density
        
        % main_figure = figure('WindowState', 'maximized');
        % outerTiled = tiledlayout(3,20);
        
        icaacttmp = EEG.icaact(ICs(j), :, :);

        spec_opt = {'freqrange', [2 50]};

        % specplot PSD
        ax_total_psd = nexttile(main_tile, 231, [3 4]);

        % cla(ax_total_psd)

        spectopo( icaacttmp(1, :), EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,ICs(j)), spec_opt{:} );
        title(ax_total_psd, 'Total Signal PSD', 'units','normalized', 'fontsize', 12, 'FontWeight', 'Normal');
        set(get(ax_total_psd, 'ylabel'), 'string', '10*log_{10}(\muV^2/Hz)', 'fontsize', 10); 
	    set(get(ax_total_psd, 'xlabel'), 'string', 'Frequency (Hz)', 'fontsize', 10, 'fontweight', 'normal'); 
	    set(ax_total_psd, 'fontsize', 10, 'fontweight', 'normal');
        xlims = xlim;
        hfreqline = findobj(ax_total_psd, 'type', 'line');
        xdata = get(hfreqline, 'xdata');
        ydata = get(hfreqline, 'ydata');
        ind = xdata >= xlims(1) & xdata <= xlims(2);
        axis on;
        axis([xlims min(ydata(ind)) max(ydata(ind))])
        box on;
        
        xline([4, 8, 13, 30], 'Color', 0.5*[1 1 1])
        grid off

        drawnow



        %% dipole plots
        ax_dipole = nexttile(main_tile, 225, [3, 6]);

        
        main_figure_dipoles = figure();
        set(gcf, 'Color', 'none')
        set(gcf, 'Position', [300, 400, 915, 290])
        main_tile_dipoles = tiledlayout(1,15, 'TileSpacing', 'compact', 'Padding', 'tight');
        % main_tile_dipoles = tiledlayout(main_tile, 1, 15, 'TileSpacing', 'compact', 'Padding', 'tight');
        % main_tile_dipoles.Layout.Tile = 225;
        % main_tile_dipoles.Layout.TileSpan = [3, 6];

        ax1_dipole = nexttile(main_tile_dipoles, 1, [1, 4]);
        axis(ax1_dipole, 'off')
        
        ax2_dipole = nexttile(main_tile_dipoles, 5, [1, 5]);
        axis(ax2_dipole, 'off')
        
        ax3_dipole = nexttile(main_tile_dipoles, 10, [1, 6]);
        axis(ax3_dipole, 'off')
        

        temp_dipole_figure = figure();
        dipplot(EEG.dipfit.model(ICs(j)), ...
            'meshdata', meshdatapath, ...
            'mri', mripath, ...
            'normlen', 'on', 'coordformat', 'MNI', ...
            'axistight', 'off', 'gui', 'off', ...
            'view', [0 0 1], 'pointout', 'off');
        
        copyobj(allchild(gca), ax1_dipole);
        % Set destination axes properties to maintain original aspect ratio
        axis(ax1_dipole, 'image');                 % locks data aspect ratio and removes stretching
        ax1_dipole.DataAspectRatioMode = 'manual'; % fixes the data aspect ratio
        ax1_dipole.PlotBoxAspectRatioMode = 'manual'; % fixes the plot-box aspect ratio
        
        
        copyobj(allchild(gca), ax2_dipole);
        % Set destination axes properties to maintain original aspect ratio
        axis(ax2_dipole, 'image');                 % locks data aspect ratio and removes stretching
        ax2_dipole.DataAspectRatioMode = 'manual'; % fixes the data aspect ratio
        ax2_dipole.PlotBoxAspectRatioMode = 'manual'; % fixes the plot-box aspect ratio
        view(ax2_dipole, [0 -1 0])


        copyobj(allchild(gca), ax3_dipole);
        % Set destination axes properties to maintain original aspect ratio
        axis(ax3_dipole, 'image');                 % locks data aspect ratio and removes stretching
        ax3_dipole.DataAspectRatioMode = 'manual'; % fixes the data aspect ratio
        ax3_dipole.PlotBoxAspectRatioMode = 'manual'; % fixes the plot-box aspect ratio
        view(ax3_dipole, [1 0 0]) 

        close(temp_dipole_figure)

        pause(0.1)

        
        % Create a panel to hold the copied figure neatly inside the tile
        panel = uipanel('Parent', Main_full_figure, 'Units', 'normalized', ...
            'Position', ax_dipole.Position, 'BorderType', 'none', 'BackgroundColor', 'w');
        axis(ax_dipole, 'off')

        % Move axes from source figure to target panel
        axes_source = findobj(main_figure_dipoles, 'Type', 'Axes');
        
        for k = 1:length(axes_source)
            new_ax = copyobj(axes_source(k), panel);
            set(new_ax, 'Units', 'normalized');
        end
        axis(ax_dipole, 'off');

        % Assuming 'panel' is the uipanel that contains copied axes:
        axes_handles = findobj(panel, 'Type', 'Axes');
        
        % Loop through each axes to find and modify line objects
        for ax = axes_handles'
            lines = findobj(ax, 'Type', 'Line');
            for ln = lines'
                % Adjust marker size and line width here
                set(ln, 'MarkerSize', 14);   % Adjust this value as needed
                set(ln, 'LineWidth', 1.5);  % Adjust this value as needed
            end
        end

        close(main_figure_dipoles)

        rv = num2str(EEG.dipfit.model(ICs(j)).rv*100, '%.1f');

        % Add a static text label inside the panel
        % uicontrol(panel, 'Style', 'text', ...
        %              'String', {['RV: ', rv, '%']}, ...
        %              'Units', 'normalized', ...
        %              'Position', [0.7 0 0.3 0.1], ...
        %              'FontSize', 10, 'BackgroundColor', 'w', ...
        %              'HorizontalAlignment', 'right');

        pos = ax_dipole.OuterPosition;
        annotation(Main_full_figure, 'textbox', pos + [0, -2*pos(2)/3, 0, 0], ...
           'String', ['RV: ', rv, '%'], ...
           'VerticalAlignment', 'bottom', ...  % or 'middle', 'bottom'
           'HorizontalAlignment', 'right', ...
           'FitBoxToText', 'on', ...
           'EdgeColor', 'none', 'FontSize', 10);

        % delete(panel)
        % cla(ax_dipole)

        % text(panel, 0.5, 0, {['RV: ', rv, '%']}, ...
        %     'fontsize', 13, 'Units', 'Normalized', 'HorizontalAlignment', 'center');



        % nexttile(outerTiled, 21, [1, 4])
        % nexttile(outerTiled, 3, [1, 2])
        % nexttile(outerTiled, 17, [1, 4])


        drawnow


        %% topoplot
        % topoplot_tl = tiledlayout(outerTiled, 1, 1);
        % topoplot_tl.Layout.Tile = 41;
        % topoplot_tl.Layout.TileSpan = [1, 3];
        ax_topoplot = nexttile(main_tile, 222, [3 3]);
        
        % cla(ax_topoplot)

        % topoplot_fig = figure();
        % set(topoplot_fig, 'Position', [600, 300, 320, 320])
        topoplot(EEG.icawinv(:,ICs(j)), EEG.chanlocs, ...
            'chaninfo', EEG.chaninfo, 'electrodes','on'); 
        axis(ax_topoplot, 'padded')

        % axis square;
        colormap(custom_cmap)
        
        icaacttmp = EEG.icaact(ICs(j), :, :);
        maxsamp = 1e5;
        n_samp = min(maxsamp, EEG.pnts*EEG.trials);
        try
            samp_ind = randperm(EEG.pnts*EEG.trials, n_samp);
        catch
            samp_ind = randperm(EEG.pnts*EEG.trials);
            samp_ind = samp_ind(1:n_samp);
        end
        if ~isempty(EEG.icachansind)
            icachansind = EEG.icachansind;
        else
            icachansind = 1:EEG.nbchan;
        end
        datavar = mean(var(EEG.data(icachansind, samp_ind), [], 2));
        projvar = mean(var(EEG.data(icachansind, samp_ind) - ...
            EEG.icawinv(:, ICs(j)) * icaacttmp(1, samp_ind), [], 2));
        pvafval = 100 *(1 - projvar/ datavar);
        pvaf = num2str(pvafval, '%3.1f');
    
        text(0.5, -0.05, {['% scalp data var.'], ['accounted for: ' pvaf '%']}, ...
            'fontsize', 10, 'Units', 'Normalized', 'HorizontalAlignment', 'center');
        

        drawnow


        %% Plot the PSD of Flexion and Extension for each Pressure
        P1_Flx_PSD_ic = squeeze(P1_Flx_PSD(ICs(j), :, :));
        P1_Ext_PSD_ic = squeeze(P1_Ext_PSD(ICs(j), :, :));

        P3_Flx_PSD_ic = squeeze(P3_Flx_PSD(ICs(j), :, :));
        P3_Ext_PSD_ic = squeeze(P3_Ext_PSD(ICs(j), :, :));

        P6_Flx_PSD_ic = squeeze(P6_Flx_PSD(ICs(j), :, :));
        P6_Ext_PSD_ic = squeeze(P6_Ext_PSD(ICs(j), :, :));

        [~, min_freq_id] = min(abs(frequencies - 2));
        [~, max_freq_id] = min(abs(frequencies - 50));


        % Flexion part
        ax_flexion_psd = nexttile(main_tile, 235, [3, 4]); hold on;
        % cla(ax_flexion_psd)
        plot(frequencies(min_freq_id:max_freq_id), ...
            mean(10*log10(P1_Flx_PSD_ic(min_freq_id:max_freq_id, :)), 2), ...
            'Color', P1_mean_color, 'LineWidth', 2)
        box on
        
        plot(frequencies(min_freq_id:max_freq_id), ...
            mean(10*log10(P3_Flx_PSD_ic(min_freq_id:max_freq_id, :)), 2), ...
            'Color', P3_mean_color, 'LineWidth', 2) 
        box on

        plot(frequencies(min_freq_id:max_freq_id), ...
            mean(10*log10(P6_Flx_PSD_ic(min_freq_id:max_freq_id, :)), 2), ...
            'Color', P6_mean_color, 'LineWidth', 2) 
        box on

        xline([4, 8, 13, 30], 'Color', 0.5*[1 1 1])

        xlim([2, 50])
        hfreqline = findobj(ax_flexion_psd, 'type', 'line');
        ydata = get(hfreqline, 'ydata');
        ydata_min_f = min(cell2mat(cellfun(@(x) min(x), ydata, 'UniformOutput', false)));
        ydata_max_f = max(cell2mat(cellfun(@(x) max(x), ydata, 'UniformOutput', false)));
        ylim([ydata_min_f ydata_max_f])
        
        ax_flexion_psd.FontSize = 10;
        ax_flexion_psd.TickLabelInterpreter = "tex";
        title('Flexion PSD', 'FontWeight', 'normal')
        xlabel('Frequency (Hz)')


        % Extension part
        ax_extension_psd = nexttile(main_tile, 239, [3, 4]); hold on;
        % cla(ax_extension_psd)
        plot(frequencies(min_freq_id:max_freq_id), ...
            mean(10*log10(P1_Ext_PSD_ic(min_freq_id:max_freq_id, :)), 2), ...
            'Color', P1_mean_color, 'LineWidth', 2) 
        box on
        
        plot(frequencies(min_freq_id:max_freq_id), ...
            mean(10*log10(P3_Ext_PSD_ic(min_freq_id:max_freq_id, :)), 2), ...
            'Color', P3_mean_color, 'LineWidth', 2) 
        box on

        plot(frequencies(min_freq_id:max_freq_id), ...
            mean(10*log10(P6_Ext_PSD_ic(min_freq_id:max_freq_id, :)), 2), ...
            'Color', P6_mean_color, 'LineWidth', 2) 
        box on

        xline([4, 8, 13, 30], 'Color', 0.5*[1 1 1])

        xlim([2, 50])
        hfreqline = findobj(ax_extension_psd, 'type', 'line');
        ydata = get(hfreqline, 'ydata');
        ydata_min_e = min(cell2mat(cellfun(@(x) min(x), ydata, 'UniformOutput', false)));
        ydata_max_e = max(cell2mat(cellfun(@(x) max(x), ydata, 'UniformOutput', false)));
        % ylim([ydata_min_e ydata_max_e])
        
        ax_extension_psd.FontSize = 10;
        ax_extension_psd.TickLabelInterpreter = "tex";
        title('Extension PSD', 'FontWeight', 'normal')
        xlabel('Frequency (Hz)')

        ydata_min_global = min(ydata_min_f, ydata_min_e);
        ydata_max_global = max(ydata_max_f, ydata_max_e);
        ylim(ax_flexion_psd, [ydata_min_global, ydata_max_global])
        ylim(ax_extension_psd, [ydata_min_global, ydata_max_global])

        legend({'P1', 'P3', 'P6', '' '' '' ''}, 'Location', 'northeast', 'Box', 'off')

        drawnow


        %% Save the figure
        figName_fig = ['sub-', num2str(subject_list(i)), ...
            '_IC', num2str(ICs(j)), '_PSD_integral.fig'];
        figName_png = ['sub-', num2str(subject_list(i)), ...
            '_IC', num2str(ICs(j)), '_PSD_integral.png'];
        folderName = [detailed_per_subject_results, ...
                '\sub-', num2str(subject_list(i))];
        if ~isfolder(folderName)
            mkdir(folderName)
        end
        saveas(Main_full_figure, fullfile(folderName, figName_fig))
        saveas(Main_full_figure, fullfile(folderName, figName_png))
        close(Main_full_figure)

    end
    % time = toc/60;
end




