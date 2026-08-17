clc
clear

%% Add and Define Necessary Paths
% main folder containing all codes and data:
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); 
% main processing and data path:
data_processing_path = [main_project_folder, '\Code\Matlab\data_processing'];
data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
% Required path for the code
EMG_features_path = [data_processing_path, ...
    '\EMG_processing\structured_EMG_data'];
cortico_muscular_analysis_path = [data_processing_path, ...
    '\Cortico_Muscular_Coherence'];
ROIs_feature_path = [data_path, '8_Classification\ROIs_features'];
epoched_data_path = [data_path, '6_Trials_Info_and_Epoched_data\'];
grouplevel_postprocessing_path = [data_processing_path, ...
    '\Group_Level_PostProcessing'];
detailed_per_subject_results = [data_processing_path, ...
    '\Detailed_per_subject_processing'];
cleaned_EEG_data_path = [data_path, '5_single-subject-EEG-analysis\'];


%% create cortical and muscular feature
% integral of PSD of each IC for each frequency band gives us the power 
% integral of EMG signals: iEMG as an indicator of total effort.

muscles_names = {'Vastus Medialis', 'Rectus Femoris', 'Gastrocnemius', 'Biceps Femoris'};
freqbands_names = {'\delta', '\theta', '\alpha', '\beta', '\gamma'};

P1_mean_color = [0 0.4470 0.7410];
P3_mean_color = [0.4660 0.6740 0.1880];
P6_mean_color = [0.6350 0.0780 0.1840];

freq_bands(1, :) = [2, 4];   % Delta
freq_bands(2, :) = [4, 8];   % Theta
freq_bands(3, :) = [8, 13];  % Alpha
freq_bands(4, :) = [13, 30]; % Beta
freq_bands(5, :) = [30, 50]; % Gamma


% Load the ROIs regardless of brain region clustering (subjects with potential brain ICs)
classisfication_path = [data_path, '8_classification\'];
ROIs_features_path = [classisfication_path, ...
    'ROIs_features\regardless_of_clustering'];
fileName = 'ROIs_0_FlextoFlex_per_trial_pressure_Based_Flx_Ext_separate_NoBrainClustering.mat';
ROIs = load(fullfile(ROIs_features_path, fileName));
name = fieldnames(ROIs);
ROIs = ROIs.(name{1});


subject_list = [5,6,7,9,11,12,13,14,15,16,17,18];


cortico_muscular_pearson_r = struct();


%%
for i = subject_list

    %% load epoched data
    disp(['Subject ', num2str(i)])

    data = load([epoched_data_path, 'sub-', num2str(i), ...
        '\Epochs_FlextoFlex_based.mat']);
    name = fieldnames(data);
    data = data.(name{1});

    Trials_Info = load([epoched_data_path, 'sub-', num2str(i), ...
        '\Trials_Info.mat']);
    name = fieldnames(Trials_Info);
    Trials_Info = Trials_Info.(name{1});


    %% Calculate the Power Spectral Density (PSD)
    [PSD, frequencies] = calculate_PSD_Flx_Ext_separate(data, Trials_Info, i);


    %% Compute the spectral power (Integral of Power Spectral Density)
    inner_struct = struct('Flexion', [], 'Extension', [], 'FlextoFlex', []);
    spectral_powers = repmat({inner_struct}, size(PSD));

    for j = 1:length(PSD)

        % Flexion part
        psd_epochs = PSD{1, j}.EEG_stream.Preprocessed.Freq_Domain.Sources.Flexion;
        if isempty(psd_epochs)
            continue
        end

        ics_bands_epochs = zeros(size(psd_epochs, 1), 5, size(psd_epochs, 3));
        for k = 1:size(freq_bands, 1)
            freq_indx = find(frequencies >= freq_bands(k,1) & frequencies <= freq_bands(k,2));

            for e = 1:size(psd_epochs, 3)
                power = trapz(frequencies(freq_indx), squeeze(psd_epochs(:, freq_indx, e)), 2);
                ics_bands_epochs(:, k, e) = power;
            end

        end
        spectral_powers{1, j}.Flexion = mean(ics_bands_epochs, 3);


        % Extension part
        psd_epochs = PSD{1, j}.EEG_stream.Preprocessed.Freq_Domain.Sources.Extension;
        
        ics_bands_epochs = zeros(size(psd_epochs, 1), 5, size(psd_epochs, 3));
        for k = 1:size(freq_bands, 1)
            freq_indx = find(frequencies >= freq_bands(k,1) & frequencies <= freq_bands(k,2));

            for e = 1:size(psd_epochs, 3)
                power = trapz(frequencies(freq_indx), squeeze(psd_epochs(:, freq_indx, e)), 2);
                ics_bands_epochs(:, k, e) = power;
            end

        end

        spectral_powers{1, j}.Extension = mean(ics_bands_epochs, 3);


        % FlextoFlex part
        psd_epochs = PSD{1, j}.EEG_stream.Preprocessed.Freq_Domain.Sources.FlextoFlex;
        
        ics_bands_epochs = zeros(size(psd_epochs, 1), 5, size(psd_epochs, 3));
        for k = 1:size(freq_bands, 1)
            freq_indx = find(frequencies >= freq_bands(k,1) & frequencies <= freq_bands(k,2));

            for e = 1:size(psd_epochs, 3)
                power = trapz(frequencies(freq_indx), squeeze(psd_epochs(:, freq_indx, e)), 2);
                ics_bands_epochs(:, k, e) = power;
            end

        end

        spectral_powers{1, j}.FlextoFlex = mean(ics_bands_epochs, 3);
    
    end


    %% Compute the iEMG of right knee extensors and flexors
    % load Structured EMG data
    disp(['subject ', num2str(i), ', muscles iEMG calculation ...'])

    filename = ['sub-', num2str(i), '\sub-', num2str(i),'_structured_EMG_data.mat'];
    structured_EMG_data = load(fullfile(EMG_features_path, filename));
    name = fieldnames(structured_EMG_data);
    structured_EMG_data = structured_EMG_data.(name{1});


    inner_struct = struct('Flexion', [], 'Extension', [], 'FlextoFlex', []);
    iEMG_muscles = repmat({inner_struct}, size(structured_EMG_data));
    for j = 1:length(structured_EMG_data)
        
        flexion_muscles_epochs = zeros(4, length(structured_EMG_data{1, j}.Signal));
        extension_muscles_epochs = zeros(4, length(structured_EMG_data{1, j}.Signal));
        flextoflex_muscles_epochs = zeros(4, length(structured_EMG_data{1, j}.Signal));

        for k = 1:length(structured_EMG_data{1, j}.Signal)
            EMG_signal = structured_EMG_data{1, j}.Signal{1, k};
            if isempty(EMG_signal)
                continue
            end

            t = data{1, j}.EMG_stream.Times{1, k};

            total_length = Trials_Info{1, j}.Events.EMG_stream.flextoflex_end_indx(k) - ...
            Trials_Info{1, j}.Events.EMG_stream.flextoflex_start_indx(k) + 1;
            
            flexion_length = Trials_Info{1, j}.Events.EMG_stream.extension_start_indx(k) - ...
            Trials_Info{1, j}.Events.EMG_stream.flextoflex_start_indx(k) + 1;

            % Flexion
            iEMG = trapz(t(1:flexion_length), EMG_signal(:, 1:flexion_length), 2);
            flexion_muscles_epochs(:, k) = iEMG;

            % Extension
            iEMG = trapz(t(flexion_length+1:end), EMG_signal(:, flexion_length+1:end), 2);
            extension_muscles_epochs(:, k) = iEMG;

            % FlextoFlex
            iEMG = trapz(t, EMG_signal, 2);
            flextoflex_muscles_epochs(:, k) = iEMG;

        end

        temp_iEMG = sum(~structured_EMG_data{1, j}.Outlier .* flexion_muscles_epochs, 2) ./ sum(~structured_EMG_data{1, j}.Outlier, 2);
        iEMG_muscles{1, j}.Flexion = temp_iEMG;
        temp_iEMG = sum(~structured_EMG_data{1, j}.Outlier .* extension_muscles_epochs, 2) ./ sum(~structured_EMG_data{1, j}.Outlier, 2);
        iEMG_muscles{1, j}.Extension = temp_iEMG;
        temp_iEMG = sum(~structured_EMG_data{1, j}.Outlier .* flextoflex_muscles_epochs, 2) ./ sum(~structured_EMG_data{1, j}.Outlier, 2);
        iEMG_muscles{1, j}.FlextoFlex = temp_iEMG;

    end



    %% Plot the main cortico-muscular features correlation based on potential Brain ICs
    % extract the ICs 
    % group the cortical powers and muscles iEMG

    ICs = cell2mat(ROIs.(['TotalBrain_sub_', num2str(i)])(:, 2));

    condition_indices = condition_indices_identifier(Trials_Info, i);

    P1_flx_EEGPower = [];
    P3_flx_EEGPower = [];
    P6_flx_EEGPower = [];

    P1_ext_EEGPower = [];
    P3_ext_EEGPower = [];
    P6_ext_EEGPower = [];


    P1_flx_iEMG = [];
    P3_flx_iEMG = [];
    P6_flx_iEMG = [];

    P1_ext_iEMG = [];
    P3_ext_iEMG = [];
    P6_ext_iEMG = [];

    for j = 1:length(spectral_powers)
        if ismember(j, condition_indices.P1)
            P1_flx_EEGPower = cat(3, P1_flx_EEGPower, spectral_powers{1, j}.Flexion);
            P1_ext_EEGPower = cat(3, P1_ext_EEGPower, spectral_powers{1, j}.Extension);
            P1_flx_iEMG = cat(2, P1_flx_iEMG, iEMG_muscles{1, j}.Flexion);
            P1_ext_iEMG = cat(2, P1_ext_iEMG, iEMG_muscles{1, j}.Extension);

        elseif ismember(j, condition_indices.P3)
            P3_flx_EEGPower = cat(3, P3_flx_EEGPower, spectral_powers{1, j}.Flexion);
            P3_ext_EEGPower = cat(3, P3_ext_EEGPower, spectral_powers{1, j}.Extension);
            P3_flx_iEMG = cat(2, P3_flx_iEMG, iEMG_muscles{1, j}.Flexion);
            P3_ext_iEMG = cat(2, P3_ext_iEMG, iEMG_muscles{1, j}.Extension);

        elseif ismember(j, condition_indices.P6)
            P6_flx_EEGPower = cat(3, P6_flx_EEGPower, spectral_powers{1, j}.Flexion);
            P6_ext_EEGPower = cat(3, P6_ext_EEGPower, spectral_powers{1, j}.Extension);
            P6_flx_iEMG = cat(2, P6_flx_iEMG, iEMG_muscles{1, j}.Flexion);
            P6_ext_iEMG = cat(2, P6_ext_iEMG, iEMG_muscles{1, j}.Extension);

        end
    end


    %% plot the final figures per IC
    table_flexion = table( ...
        'Size',      [0, 5], ...
        'VariableTypes', {'double','string','string','double','double'}, ...
        'VariableNames', {'IC','Frequency Band','Muscle','r','p value'});
    table_extension = table( ...
        'Size',      [0, 5], ...
        'VariableTypes', {'double','string','string','double','double'}, ...
        'VariableNames', {'IC','Frequency Band','Muscle','r','p value'});

    for ic = 1:length(ICs)
        
        %% create the frame - Flexion
        % Main_figure_flexion = figure('WindowState', 'maximized');
        % Main_layout = tiledlayout(1, 4);
        % title(Main_layout, ['Subject ', num2str(i), ' - Flexion Part', ' - IC ', num2str(ICs(ic))], ...
        %     'FontSize', 14)

        for muscle = 1:4 % four muscles
            
            % row_layout = tiledlayout(Main_layout, 5, 1);
            % row_layout.Layout.Tile = muscle;   
            % xlabel(row_layout, muscles_names{muscle})

            spec_order = [5, 4, 3, 2, 1];
            for spect = 1:5 % five spectral bands
                

                % mark outliers
                S = 2;
                
                x1 = 1000*[P1_flx_iEMG(muscle, :)];
                y1 = squeeze(P1_flx_EEGPower(ICs(ic), spect, :));
                mean_P1 = mean(y1);
                std_P1  = std(y1);
                outlierIdx_P1 = or(y1<(mean_P1-S*std_P1) , y1>(mean_P1+S*std_P1));

                x3 = 1000*[P3_flx_iEMG(muscle, :)];
                y3 = squeeze(P3_flx_EEGPower(ICs(ic), spect, :));
                mean_P3 = mean(y3);
                std_P3  = std(y3);
                outlierIdx_P3 = or(y3<(mean_P3-S*std_P3) , y3>(mean_P3+S*std_P3));

                x6 = 1000*[P6_flx_iEMG(muscle, :)];
                y6 = squeeze(P6_flx_EEGPower(ICs(ic), spect, :));
                mean_P6 = mean(y6);
                std_P6  = std(y6);
                outlierIdx_P6 = or(y6<(mean_P6-S*std_P6) , y6>(mean_P6+S*std_P6));



                % ax = nexttile(row_layout, spec_order(spect)); hold on;
                % 
                % scatter(x1(~outlierIdx_P1)', y1(~outlierIdx_P1), ...
                %     10, P1_mean_color, "filled", "o")
                % scatter(x3(~outlierIdx_P3)', y3(~outlierIdx_P3), ...
                %     10, P3_mean_color, "filled", "o")
                % scatter(x6(~outlierIdx_P6)', y6(~outlierIdx_P6), ...
                %     10, P6_mean_color, "filled", "o")

                
                % if muscle == 1
                %     ylabel(ax, [freqbands_names{spect},  ' Power (\muV^2)'])
                % end
                % if spect == 1
                %     xlabel(ax, 'iEMG (mV.second)')
                % end



                x_clean = [x1(~outlierIdx_P1), x3(~outlierIdx_P3), x6(~outlierIdx_P6)]';
                y_clean = [y1(~outlierIdx_P1); y3(~outlierIdx_P3); y6(~outlierIdx_P6)];

                % mdl = fitlm(x_clean, y_clean, 'linear');
                % 
                % [xs, ~] = sort(x_clean); 
                % 
                % % get the fit + CI at each xs
                % [yfit, yci] = predict(mdl, xs);     
                % 
                % % solid fit line
                % plot(xs, yfit,  'k-',  'LineWidth',2);  
                % % dashed CIs
                % plot(xs, yci(:,1), 'k--', 'LineWidth',1);
                % plot(xs, yci(:,2), 'k--', 'LineWidth',1);


                [x,y] = deal(x_clean, y_clean); 
                [r, pval] = corr(x, y, 'Type','Pearson');
                % txt = sprintf('r = %.2f, p = %.2g', r, pval);
                % 
                % title(ax, txt, 'FontWeight', 'normal')


                % R2   = mdl.Rsquared.Ordinary;
                % txt = sprintf('R^2 = %.1f%%', 100*R2);

                % % put it at 5% in from the left, 90% up from the bottom
                % text(0.6, 0.95, txt, ...
                %      'Units','normalized', ...
                %      'FontSize', 10, ...
                %      'FontWeight', 'normal', ...
                %      'Interpreter','tex',...
                %      'FontAngle','italic',...
                %      'HorizontalAlignment','right',...
                %      'BackgroundColor', 'w', ...
                %      'EdgeColor','none', ...
                %      'VerticalAlignment','top', ...
                %      'HorizontalAlignment','left');


                if abs(r) > 0.5
                    table_flexion(end+1, :) = {ICs(ic), ...
                        freqbands_names{spect}(2:end), ...
                        muscles_names{muscle}, r, pval};
                end
                
                


            end
        end

        

        % legend({'P1', 'P3', 'P6'}, 'Box', 'off', 'Location', 'best')


        %% create the frame - Extension
        % Main_figure_extension = figure('WindowState', 'maximized');
        % Main_layout = tiledlayout(1, 4);
        % title(Main_layout, ['Subject ', num2str(i), ' - Extension Part', ' - IC ', num2str(ICs(ic))], ...
        %     'FontSize', 14)

        for muscle = 1:4 % four muscles
            
            % row_layout = tiledlayout(Main_layout, 5, 1);
            % row_layout.Layout.Tile = muscle;   
            % xlabel(row_layout, muscles_names{muscle})

            % spec_order = [5, 4, 3, 2, 1];
            for spect = 1:5 % five spectral bands
                

                % mark outliers
                S = 2;
                
                x1 = 1000*[P1_ext_iEMG(muscle, :)];
                y1 = squeeze(P1_ext_EEGPower(ICs(ic), spect, :));
                mean_P1 = mean(y1);
                std_P1  = std(y1);
                outlierIdx_P1 = or(y1<(mean_P1-S*std_P1) , y1>(mean_P1+S*std_P1));

                x3 = 1000*[P3_ext_iEMG(muscle, :)];
                y3 = squeeze(P3_ext_EEGPower(ICs(ic), spect, :));
                mean_P3 = mean(y3);
                std_P3  = std(y3);
                outlierIdx_P3 = or(y3<(mean_P3-S*std_P3) , y3>(mean_P3+S*std_P3));

                x6 = 1000*[P6_ext_iEMG(muscle, :)];
                y6 = squeeze(P6_ext_EEGPower(ICs(ic), spect, :));
                mean_P6 = mean(y6);
                std_P6  = std(y6);
                outlierIdx_P6 = or(y6<(mean_P6-S*std_P6) , y6>(mean_P6+S*std_P6));



                % ax = nexttile(row_layout, spec_order(spect)); hold on;
                % 
                % scatter(x1(~outlierIdx_P1)', y1(~outlierIdx_P1), ...
                %     10, P1_mean_color, "filled", "o")
                % scatter(x3(~outlierIdx_P3)', y3(~outlierIdx_P3), ...
                %     10, P3_mean_color, "filled", "o")
                % scatter(x6(~outlierIdx_P6)', y6(~outlierIdx_P6), ...
                %     10, P6_mean_color, "filled", "o")

                
                % if muscle == 1
                %     ylabel(ax, [freqbands_names{spect},  ' Power (\muV^2)'])
                % end
                % if spect == 1
                %     xlabel(ax, 'iEMG (mV.second)')
                % end



                x_clean = [x1(~outlierIdx_P1), x3(~outlierIdx_P3), x6(~outlierIdx_P6)]';
                y_clean = [y1(~outlierIdx_P1); y3(~outlierIdx_P3); y6(~outlierIdx_P6)];

                % mdl = fitlm(x_clean, y_clean, 'linear');
                % 
                % [xs, ~] = sort(x_clean); 
                % 
                % % get the fit + CI at each xs
                % [yfit, yci] = predict(mdl, xs);     
                % 
                % % solid fit line
                % plot(xs, yfit,  'k-',  'LineWidth',2);  
                % % dashed CIs
                % plot(xs, yci(:,1), 'k--', 'LineWidth',1);
                % plot(xs, yci(:,2), 'k--', 'LineWidth',1);


                [x,y] = deal(x_clean, y_clean); 
                [r, pval] = corr(x, y, 'Type','Pearson');
                % txt = sprintf('r = %.2f\np = %.2g', r, pval);
                % 
                % title(ax, txt, 'FontWeight', 'normal')


                % R2   = mdl.Rsquared.Ordinary;
                % txt = sprintf('R^2 = %.1f%%', 100*R2);

                % % put it at 5% in from the left, 90% up from the bottom
                % text(0.6, 0.95, txt, ...
                %      'Units','normalized', ...
                %      'FontSize', 10, ...
                %      'FontWeight', 'normal', ...
                %      'Interpreter','tex',...
                %      'FontAngle','italic',...
                %      'HorizontalAlignment','right',...
                %      'BackgroundColor', 'w', ...
                %      'EdgeColor','none', ...
                %      'VerticalAlignment','top', ...
                %      'HorizontalAlignment','left');
                
                if abs(r) > 0.5
                    table_extension(end+1, :) = {ICs(ic), ...
                        freqbands_names{spect}(2:end), ...
                        muscles_names{muscle}, r, pval};
                end
                


            end
        end
        


        %% save figures

        % outFolder = [cortico_muscular_analysis_path, '\sub-', num2str(i)];
        % 
        % if ~exist(outFolder, 'dir')
        %     mkdir(outFolder);
        % end
        % 
        % figName = ['Subject ', num2str(i), ' - Flexion Part', ' - IC ', num2str(ICs(ic))];
        % savefig(Main_figure_flexion, fullfile(outFolder, [figName '.fig']));
        % saveas(Main_figure_flexion, fullfile(outFolder, [figName '.png']));
        % 
        % 
        % figName = ['Subject ', num2str(i), ' - Extension Part', ' - IC ', num2str(ICs(ic))];
        % savefig(Main_figure_extension, fullfile(outFolder, [figName '.fig']));
        % saveas(Main_figure_extension, fullfile(outFolder, [figName '.png']));



        % close all



    end


    cortico_muscular_pearson_r.(['sub_', num2str(i)]).Flexion = table_flexion;
    cortico_muscular_pearson_r.(['sub_', num2str(i)]).Extension = table_extension;




end





%%
% figure(); hold on;
% plot(frequencies, PSD{1, 20}.EEG_stream.Preprocessed.Freq_Domain.Sources.Flexion(1, :, 1));
% plot(frequencies, PSD{1, 20}.EEG_stream.Preprocessed.Freq_Domain.Sources.Extension(1, :, 1));
% plot(frequencies, PSD{1, 20}.EEG_stream.Preprocessed.Freq_Domain.Sources.FlextoFlex(1, :, 1));
% xlim([2, 50])

% figure(); hold on;
% plot(spectral_powers{1, 20}.Flexion(1, :), 'Marker', 'o', 'LineStyle', 'none');
% plot(spectral_powers{1, 20}.Extension(1, :), 'Marker', 'o', 'LineStyle', 'none');
% plot(spectral_powers{1, 20}.FlextoFlex(1, :), 'Marker', 'o', 'LineStyle', 'none');