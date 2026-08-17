clc
clear


%% Add necessary paths
main_project_path = 'D:\Morteza\MyProjects\ANSYMB2024\';

addpath(genpath([main_project_path, 'Code']));
addpath(genpath([main_project_path, 'data\7_STUDY\Epoched_data']));

data_path         = [main_project_path, 'data\'];
Code_path         = [main_project_path, 'Code\Matlab\data_processing\'];
all_STUDY_PATH    = [data_path, '7_STUDY\Epoched_data\', ...
                        'multiple_clustering\'];

icatimef_path     = [data_path, '5_single-subject-EEG-analysis\', ...
                        'timewarp_test\Epoched_data'];
epoched_data_path = [data_path, '6_Trials_Info_and_Epoched_data\'];
ersp_data_path    = [data_path, '7_STUDY\Epoched_data\Final_figures', ...
                        '\ERSP\Three Pressure Conditions\', ...
                        'p 0.01 ersp results\'];
Subject_ICs_in_clusters_path = [Code_path, ...
    'Group_Level_PostProcessing\Final_paper_plot_generation\', ...
    'Detailed_Analysis_on_TF_regions\', ...
    'extracting Subjects and ICs in the brain clusters'];

current_path = ['D:\Morteza\MyProjects\ANSYMB2024\Code\Matlab', ...
    '\data_processing\Group_Level_PostProcessing\', ...
    'Final_paper_plot_generation\', ...
    '_NHB\sPCA_denoising_of_TF_data'];

titles = {'Low Pressure', 'Medium Pressure', 'High Pressure'};
studyNames = {'Left Dorsal ACC', 'Left Parieto Occipital', ...
    'Left PreMot SuppMot', 'Left Prim Motor', 'Prime Visual', ...
    'Right Parieto Occipital', 'Right PreMot SuppMot', 'Right Prim Motor'};



%% ========================================================================
%                          Main Loop on Studies
%  ========================================================================

% load the subject-IC pairs for our STUDY
cd(Subject_ICs_in_clusters_path)
load("Subjects_ICs_in_clusters.mat")
cd(current_path)

for study = 1:length(studyNames)

    disp([studyNames{study}, ' ...'])

    idx_cluster = find(cellfun(@(x) strcmp(x, SUBJECTS_ICS{study, 1}), ...
        SUBJECTS_ICS(:, 1)));
    Subjects = SUBJECTS_ICS{idx_cluster, 2}.Subjects + 4;
    [Subjects_sorted, idx_subject_sort] = sort(Subjects, 2, "ascend");

    ICs                    = SUBJECTS_ICS{idx_cluster, 2}.ICs;
    ICs_sorted             = ICs(idx_subject_sort);


    for sub = 1:length(Subjects_sorted)

        cd(icatimef_path)
        fileExt = '.icatimef';
        
        % load icatime (just load the com_X, times, freqs)
        fileBaseName = ['S', num2str(Subjects_sorted(sub))];
        chanList = ['comp', num2str(ICs_sorted(sub))];
        disp(['Loading S', num2str(Subjects_sorted(sub)),'.icatimef ...']);
        icatimef = load('-mat', [ fileBaseName fileExt ], chanList, ...
            'times', 'freqs', 'trialinfo', 'parameters');
        trialinfo = icatimef.trialinfo;
        ic = icatimef.(['comp', num2str(ICs_sorted(sub))]);
        ic = ic.*conj(ic);
        % crop the ic data based on the timewarpms 
        idx = find(strcmp(icatimef.parameters, 'timewarpms'));
        timewarpms = icatimef.parameters{1, idx+1};
        times = icatimef.times;
        idx_to_keep = times < timewarpms(end);
        ic = ic(:, idx_to_keep, :);
        new_times = 100*(times(idx_to_keep)/timewarpms(end));


        cd(current_path)
        
       

        %% ===============================================
        %  ERSP Calculation (ERSP: nFreq x nTime x nEpoch)
        %  ===============================================
        
        tf_cells = squeeze(num2cell(ic, [1 2]));   
        tbl = table(tf_cells, 'VariableNames', {'tf'});  
        nEpoch = height(tbl);
        
        % trial numbers in tbl
        tbl.trial = cellfun(@(x) str2num(x), {trialinfo.trial}');
        % trial condition in tbl
        tbl.cond = cellfun(@(x) str2num(x), {trialinfo.cond}');
        
        % creating baseline
        mean_baselines = cell(1, 3);
        C = [1, 3, 6];
        for i = 1:3
            idx_pX = find(tbl.cond == C(i));
            pX_epochs = cat(3, tbl.tf{idx_pX});
            mean_baselines{i} = mean(pX_epochs, 3);
        end
        main_baseline = (1/3)*(sum(cat(3, mean_baselines{:}), 3));
        main_baseline = mean(main_baseline, 2);
        
        % compute the ERSP per epoch
        % --- preallocate
        % trial ERSP will be [nFreq x nTime x nEpoch]
        tf0 = tbl.tf{1};
        [nFreq, nTime] = size(tf0);
        ERSP_epoch = nan(nFreq, nTime, nEpoch);
        
        B = repmat(main_baseline, 1, nTime);
        for iE = 1:nEpoch
            P = tbl.tf{iE};
            ERSP_epoch(:, :, iE) = P ./ B;    
        end

        

        %% ===============================================
        %  Spectral PCA across observations (trial x time)
        %  ===============================================
        
        % Build X: rows are (time, trial), columns are frequency
        % ERSP_trial: 
        % [F x T x R] -> permute to [T x R x F] -> reshape to [(T*R) x F]
        X = reshape(permute(ERSP_epoch, [2 3 1]), [], nFreq);
        
        % mean-center per frequency
        mu = mean(X, 1);
        Xc = X - mu;
        
        % PCA via SVD (no toolbox dependency)
        [~, S, V] = svd(Xc, 'econ');
        
        pc1   = V(:,1);            % [nFreq x 1] dominant spectral pattern
        score1 = Xc * pc1;         % [(T*R) x 1]
        
        % remove PC1 (rank-1 reconstruction)
        Xc_den = Xc - score1 * pc1';
        
        % optional: explained variance of PC1
        singvals = diag(S);
        expVar = (singvals.^2) / sum(singvals.^2);
        expVar1 = expVar(1);
        
        % reshape back to [F x T x R]
        X_den = Xc_den + mu;
        ERSP_den = permute(reshape(X_den, nTime, nEpoch, nFreq), [3 1 2]);
        
        % store the conditined data
        ERSP_epoch_cond = cell(1, 3);
        ERSP_den_cond   = cell(1, 3);
        for i = 1:3
            idx_pX = find(tbl.cond == C(i));
            ERSP_epoch_cond{1, i} = ERSP_epoch(:, :, idx_pX);
            ERSP_den_cond{1, i} = ERSP_den(:, :, idx_pX);
        end
        
        % find the best climits
        data = [];
        for i = 1:3
            data = [data, ...
                reshape(10*log10(mean(ERSP_epoch_cond{1, i}, 3)).', 1, [])];
            data = [data, ...
                reshape(10*log10(mean(ERSP_den_cond{1, i}, 3)).', 1, [])];
        end
        idx = find(imag(data));
        data(idx) = 0;
        IQR = iqr(data); % interquartile range
        Q1 = quantile(data,0.25);
        myMin = round(Q1-1.5*IQR,1);
        erspdata_clim = [myMin myMin*(-1)];


        %% =========================================================
        %  Plot the Subject IC TF with and without specPCA deniosing
        %  =========================================================

        % figure()
        monitors = get(0, 'MonitorPositions');
        fig = figure('name', ['TF with and without specPCA denoising'], ...
            'InvertHardcopy', 'off', 'PaperType', 'a2', ...
            'PaperOrientation', 'landscape', ...
            'Resize', 'off');
        
        % For second monitor (row 2), add drawnow before setting position
        drawnow;  % Let MATLAB finish drawing on primary monitor first
        pause(0.1);  % Short pause helps
        
        set(fig, 'Position', [monitors(1,1)+100, monitors(1,2)+700, 1100, 800]);
        for i = 1:3
            fh(1, i).h = subplot(2, 3, i);
            tf = 10*log10(mean(ERSP_epoch_cond{1, i}, 3));
            idx = find(imag(tf) ~= 0);
            tf(idx) = 0;
            contourf(new_times, icatimef.freqs, tf, 200, ...
                'linecolor','none')
            set(gca,'clim', erspdata_clim, ...
                'xlim', [new_times(1) new_times(end)], ...
                'ydir', 'norm', ...
                'ylim', [icatimef.freqs(1) icatimef.freqs(200)], ...
                'yscale','log')
            set(gcf, 'Colormap', calldefinedcolormap(), 'Color', [1 1 1]);
            
            fh(2, i).h = subplot(2, 3, i+3);
            tf = 10*log10(mean(ERSP_den_cond{1, i}, 3));
            idx = find(imag(tf) ~= 0);
            tf(idx) = 0;
            contourf(new_times, icatimef.freqs, tf, 200, ...
                'linecolor','none')
            set(gca, 'clim', erspdata_clim, ...
                'xlim', [new_times(1) new_times(end)], ...
                'ydir', 'norm', ...
                'ylim', [icatimef.freqs(1) icatimef.freqs(200)], ...
                'yscale', 'log')
            set(gcf, 'Colormap', calldefinedcolormap(), 'Color', [1 1 1]);
        
            % resize plot to fit title
            pos = fh(1, i).h.Position;
            fh(1, i).h.Position = [pos(1) pos(2)-0.07 pos(3)*.9 pos(4)*.8];
        
            pos = fh(2, i).h.Position;
            fh(2, i).h.Position = [pos(1) pos(2) pos(3)*.9 pos(4)*.8];
        
        
            if i == 3
                for j = 1:2
                    pos = fh(j, i).h.Position;
                    c = colorbar('Position', ...
                        [pos(1)+pos(3)+0.01  pos(2) 0.01 pos(4)]);
                    c.Limits = clim;
                    % make the Ticks symmetric
                    maxAbs = max(abs(c.Ticks));
                    % If maxAbs isn't in v, append it to make symmetric
                    if ~ismember(maxAbs, c.Ticks)
                        c.Ticks = sort([c.Ticks(1) 0 maxAbs]); 
                    end
                    hL = ylabel(c,[{'Power (dB)'}],...
                        'fontweight', 'bold', 'FontName', 'Arial', ...
                        'FontSize', 14, 'Rotation', 90);
                    hL.Position(1) = 4;
                    hL.Position(2) = 0;
                end
            end
        
            for j = 1:2
                axes(fh(1, i).h)
                set(gca, ...
                    'XTick',[new_times(1) new_times(66) new_times(end)], ...
                    'XTickLabel', {'0', '50', '100'}, ...
                    'ytick', [4 8 14 30 60], ...
                    'fontsize',10);
                xtickangle(45)
                set(gca,'Fontsize',16);
                axes(fh(2, i).h)
                set(gca, ...
                    'XTick',[new_times(1) new_times(66) new_times(end)], ...
                    'XTickLabel', {'0', '50', '100'}, ...
                    'ytick', [4 8 14 30 60], ...
                    'fontsize',10);
                xtickangle(45)
                set(gca,'Fontsize',16);
            end
            
            % ylabel
            if i == 1
                axes(fh(1, i).h)
                ylh = ylabel(sprintf('No sPCA denoising\n\nFrequency (Hz)'), ...
                    'fontsize', 16, 'fontweight', 'bold', 'FontName', 'Arial');
                ylh.Position(1) = ylh.Position(1); 
        
                axes(fh(2, i).h)
                ylh = ylabel(sprintf('With sPCA denoising\n\nFrequency (Hz)'), ...
                    'fontsize', 16, 'fontweight', 'bold', 'FontName', 'Arial');
                ylh.Position(1) = ylh.Position(1); 
            else
                set(gca,'YTickLabel',[]);
            end
            
            xlh = xlabel('Cycle (%)','Fontsize',16,'fontweight','bold');
            xlh.Position(2) = 1.5;
        
            
            axes(fh(1, i).h)
            T = title(titles{i},'FontSize',16);
            T.Position(2) = T.Position(2)+80;
        
        
            for j = 1:2
                axes(fh(j, i).h)
                evPlotLines = [new_times(1) new_times(66) new_times(end)];
                eventLabels = ...
                    {sprintf('FlxS'), sprintf('FlxE\nExtS'), sprintf('ExtE')};
                hold on;
                for L = 1:length(evPlotLines)
                    if L == 1 || L == length(evPlotLines)
                        v = vline(evPlotLines(L),'-k', eventLabels{1,L}); 
                        set(v,'LineWidth',1); % solid line
                    else
                        v = vline(evPlotLines(L),':k',eventLabels{1,L}); 
                        set(v,'LineWidth',1.2);
                    end
                end
            
                H = findobj(gca);
                tb = findobj(H,'Type','text');
        
                for textbox = 1:3 % 1:size(tb,1)
                    if     mod(textbox, 3) == 1
                        % Ext End
                        pos = tb(textbox).Position;
                        tb(textbox).Position = [pos(1)+3 70 0];
                        set(tb(textbox),'Rotation',90) % rotate 90 degrees
                        set(tb(textbox),'FontSize',8, 'FontWeight', 'bold') 
                    elseif mod(textbox, 3) == 2
                        pos = tb(textbox).Position;
                        tb(textbox).Position = [pos(1)-1 70 0];
                        set(tb(textbox),'Rotation',90) % rotate 90 degrees
                        set(tb(textbox),'FontSize',8, 'FontWeight', 'bold') 
                    elseif mod(textbox, 3) == 0
                        pos = tb(textbox).Position;
                        tb(textbox).Position = [pos(1)+1 70 0];
                        set(tb(textbox),'Rotation',90) % rotate 90 degrees
                        set(tb(textbox),'FontSize',8, 'FontWeight', 'bold') 
                    end
                end
            end
        
        end

        mainT = sgtitle(['S', num2str(Subjects_sorted(sub)), ...
                 ' IC', num2str(ICs_sorted(sub)), ...
                 ' - ', studyNames{study}]);
        mainT.FontSize = 18;
        mainT.FontWeight = "bold";
        mainT.FontName = 'Arial';


        % save the figure
        cd([current_path, '\figures'])
        folderName = studyNames{study};
        if ~exist(folderName, "dir")
            mkdir(folderName)
        end
        
        figname = mainT.String;
        savethisfig(gcf, strcat(figname,'.png'), ...
            [current_path,'\figures\', folderName,'\png'],'png')
        savethisfig(gcf, strcat(figname,'.fig'), ...
            [current_path,'\figures\', folderName,'\fig'],'fig')
        savethisfig(gcf, strcat(figname,'.svg'), ...
            [current_path,'\figures\', folderName,'\svg'],'svg')
        
        close;
        clear tbl

    end

end