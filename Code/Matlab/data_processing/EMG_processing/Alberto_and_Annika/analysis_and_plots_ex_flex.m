%% Load data
path = '/Users/annika/Documents/MATLAB/EMG/Matlab/sub-18'; % adjust path and subject!
load(fullfile(path, 'Epochs_Extension_based.mat'))
load(fullfile(path, 'Epochs_Flexion_based.mat'))

data_ex = Epochs_Extension_based;
data_flex = Epochs_Flexion_based;
data = [data_ex; data_flex]; 


%% Errors, conditions, search for session indexes
index_mistakes_ex = [];
index_mistakes_flex = [];
index_session = [];
index_cond = cell(1,3);

for b = 1:2 % ex or flex
    for d = 1:length(data) % implementation
        implementation = data{b,d}.EMG_stream.Sensors_Preprocessed;

        % remove if implementation has no data
        if isempty(implementation)
            if b == 1
                index_mistakes_ex = [index_mistakes_ex, d];
            elseif b == 2
                index_mistakes_flex = [index_mistakes_flex, d];
            end
            continue
        end

        % save index of session-start in vector
        % -> same for ex and flex
        session = data{b,d}.General.Session;
        if ~ismember(session,index_session)
           index_session = [index_session; session, d];
        end
       

        % divide for conditions (NoPam = 1)
        % -> same for ex and flex
        description = data{b,d}.General.Description;
        if contains(description, 'Experiment')
            pressure = data{b,d}.General.Pressure;
            if pressure == 1
                index_cond{1,1}{end+1} = d;
            elseif pressure == 3
                index_cond{1,2}{end+1} = d;
            elseif pressure == 6
                index_cond{1,3}{end+1} = d;
            end
        elseif contains(description, 'No_PAM')
            index_cond{1,1}{end+1} = d;
        end

        % search for errors in epochs
        for ep = 1:length(implementation) % epochs
            one_swing = implementation{1,ep};    % 6x... matrix
            if length(one_swing) < 50
                if b == 1
                    index_mistakes_ex = [index_mistakes_ex, d];
                elseif b == 2
                    index_mistakes_flex = [index_mistakes_flex, d];
                end
            end
        end
    end
end

% index of error-trials
index_trial_mistakes = cell(1,2);
index_trial_mistakes{1,1} = unique(index_mistakes_ex);
index_trial_mistakes{1,2} = unique(index_mistakes_flex);


%% Save newly sorted data without errors in trials
for b = 1:2
    for ses = 1:4
        for sen = 1:4
            exflex_session_sensor_cond{1,b}{1,ses}{1,sen} = cell(1,3);
        end
    end
end


for b = 1:2
    for d = 1:length(data) % implementation
        implementation = data{b,d}.EMG_stream.Sensors_Preprocessed;
        
        % skip if implementation has errors
        if ismember(d,index_trial_mistakes{1,b})
            continue
        end
        
        % divide data in sessions with previously created session-index
        for ses = 1:3 % sessions 1-3
            if d >= index_session(ses,2) && d < index_session(ses+1,2)
                % conditions
                for con =1:3
                    if ismember(d,cell2mat(index_cond{1,con}))
                        for ep = 1:length(implementation) % epochs
                            one_swing = implementation{1,ep};
                            % sensors
                            for sen = 1:4
                                exflex_session_sensor_cond{1,b}{1,ses}{1,sen}{1,con}{end+1,1} = one_swing(sen,:);
                            end
                        end
                    end
                end 
            end
        end
        if d >= index_session(4,2) % session 4
            % conditions
            for con =1:3
                if ismember(d,cell2mat(index_cond{1,con}))
                    for ep = 1:length(implementation) % epochs
                        one_swing = implementation{1,ep};
                        % sensors
                        for sen = 1:4
                            exflex_session_sensor_cond{1,b}{1,4}{1,sen}{1,con}{end+1,1} = one_swing(sen,:);
                        end
                    end
                end
            end
        end
    end
end


%% Outlier detection for every ex/flex/session/sensor/cond individually
% segment epochs and calculate mean
for b = 1:2
    for ses = 1:5 % 4 sessions + total
        for sen = 1:4
            exflex_session_sensor_cond_seg{1,b}{1,ses}{1,sen} = cell(1,3);
        end
    end
end
epochs_seg_mat = [];
epochs_equal_mat_total = [];

for b = 1:2 % ex/flex
    for ses = 1:4 % session
        for sen = 1:4 % sensor
            for con = 1:3 % condition
                for ep = 1:size(exflex_session_sensor_cond{1,b}{1,ses}{1,sen}{1,con},1)
                    l = length(exflex_session_sensor_cond{1,b}{1,ses}{1,sen}{1,con}{ep,1}); % length of epoch
                    % length must be dividable by 10
                    while mod(l, 10) ~= 0
                        l = l+1;
                    end
                    laenge_seg = l/10;
                    % divide epochs into 10 segmentes and take mean
                    epoch = exflex_session_sensor_cond{1,b}{1,ses}{1,sen}{1,con}{ep,1};
                    seg_vec = [];
                    for seg = 1:10
                        if seg ~= 10
                            segment = epoch(1,laenge_seg*(seg-1)+1:laenge_seg*seg);
                            seg_mean = mean(segment);
                            seg_vec = [seg_vec, seg_mean];
                        else
                            segment = epoch(1,laenge_seg*(seg-1)+1:length(epoch));
                            seg_mean = mean(segment);
                            seg_vec = [seg_vec, seg_mean];
                        end
                    end
                    epochs_seg_mat = [epochs_seg_mat;seg_vec];
                end
                exflex_session_sensor_cond_seg{1,b}{1,ses}{1,sen}{1,con} = epochs_seg_mat;
                % for total
                epochs_equal_mat_total = exflex_session_sensor_cond_seg{1,b}{1,5}{1,sen}{1,con};
                epochs_equal_mat_total = [epochs_equal_mat_total;epochs_seg_mat];
                exflex_session_sensor_cond_seg{1,b}{1,5}{1,sen}{1,con} = epochs_equal_mat_total;
                epochs_seg_mat = [];
            end
        end
    end
end

% calculate robust covariance and mahalanobis distance
for b = 1:2 
    for ses = 1:5 % 4 sessions + total
        for sen = 1:4
            exflex_session_sensor_cond_outlier{1,b}{1,ses}{1,sen} = cell(1,3);
        end
    end
end

for b = 1:2
    for ses = 1:5 % session
        for sen = 1:4 % sensor
            for con = 1:3 % condition
                [sig,mu,mah,outliers] = robustcov(exflex_session_sensor_cond_seg{1,b}{1,ses}{1,sen}{1,con});
                exflex_session_sensor_cond_outlier{1,b}{1,ses}{1,sen}{1,con} = outliers;
            end
        end
    end
end


%% Bring to the same length
% calculate median
all_len = [];
exflex_med = [];
for b = 1:2
    for ses = 1:4
        for sen = 1:4
            for con = 1:3
                for ep = 1:length(exflex_session_sensor_cond{1,b}{1,ses}{1,sen}{1,con})
                    l = length(exflex_session_sensor_cond{1,b}{1,ses}{1,sen}{1,con}{ep,1});
                    all_len = [all_len, l];
                end
            end
        end
    end
    med = round(median(all_len));
    exflex_med = [exflex_med,med];    
end


% bring all to same length (original signal + outliers out)
for o = 1:2 % original , ohne outlier
    for b = 1:2
        for ses = 1:5 % 4 sessions + total
            for sen = 1:4
                exflex_session_sensor_cond_equal{1,o}{1,b}{1,ses}{1,sen} = cell(1,3); 
            end
        end
    end
end
epochs_equal_mat = [];
epochs_equal_mat_out = [];
epochs_equal_mat_total = [];

for  o = 1:2
    for b = 1:2
        for ses = 1:4
            for sen = 1:4
                for con = 1:3
                    for ep = 1:length(exflex_session_sensor_cond{1,b}{1,ses}{1,sen}{1,con})
                        emg_signal = exflex_session_sensor_cond{1,b}{1,ses}{1,sen}{1,con}{ep,1};
                        int_signal = interp1(1:length(emg_signal), emg_signal, linspace(1, length(emg_signal), exflex_med(1,b)), 'linear');
                        epochs_equal_mat = [epochs_equal_mat;int_signal];
                        if exflex_session_sensor_cond_outlier{1,b}{1,ses}{1,sen}{1,con}(ep,1) == 0
                            epochs_equal_mat_out = [epochs_equal_mat_out;int_signal];
                        end
                    end
                    exflex_session_sensor_cond_equal{1,1}{1,b}{1,ses}{1,sen}{1,con} = epochs_equal_mat;
                    exflex_session_sensor_cond_equal{1,2}{1,b}{1,ses}{1,sen}{1,con} = epochs_equal_mat_out;
        
                    % for total
                    epochs_equal_mat_total = exflex_session_sensor_cond_equal{1,1}{1,b}{1,5}{1,sen}{1,con};
                    epochs_equal_mat_total = [epochs_equal_mat_total;epochs_equal_mat];
                    exflex_session_sensor_cond_equal{1,1}{1,b}{1,5}{1,sen}{1,con} = epochs_equal_mat_total;
        
                    epochs_equal_mat_total = exflex_session_sensor_cond_equal{1,2}{1,b}{1,5}{1,sen}{1,con};
                    epochs_equal_mat_total = [epochs_equal_mat_total;epochs_equal_mat_out];
                    exflex_session_sensor_cond_equal{1,2}{1,b}{1,5}{1,sen}{1,con} = epochs_equal_mat_total;
        
                    % empty matrices for next condition
                    epochs_equal_mat = [];
                    epochs_equal_mat_out = [];
                end
            end
        end
        epochs_equal_mat_total = [];
    end
end




%% Plot
% creates folder in path
status = mkdir(path, 'plots');

% calculate mean and std
for ses = 1:5 % 4 sessions + total
    for o = 1:2 % with outliers and without outliers
    max_signal{1,o}{1,ses} = 0;
    min_signal{1,o}{1,ses} = 0;
        for b = 1:2
            for sen = 1:4
                out_exflex_session_sensor_cond_mean_std{1,o}{1,b}{1,ses}{1,sen} = cell(1,3);
            end
        end
    end
end

for o = 1:2
    for ses = 1:5
        for b = 1:2
            for sen = 1:4
                for con = 1:3
                    
                    data = exflex_session_sensor_cond_equal{1,o}{1,b}{1,ses}{1,sen}{1,con};
                    m = mean(data);     % mean of each column
                    s = std(data);
                    n = numel(data);
                    sem = s / sqrt(n);
    
                    out_exflex_session_sensor_cond_mean_std{1,o}{1,b}{1,ses}{1,sen}{1,con}{1,1} = m;
                    out_exflex_session_sensor_cond_mean_std{1,o}{1,b}{1,ses}{1,sen}{1,con}{1,2} = sem;
    
                    % max std
                    if max(m + sem) > max_signal{1,o}{1,ses}
                        max_signal{1,o}{1,ses} = max(m + sem);
                    end
                    if min(m + sem) < min_signal{1,o}{1,ses}
                        min_signal{1,o}{1,ses} = min(m + sem);
                    end
                end
            end
        end
    end
end


% definitions
muscle = ["Knee Extensor", "Knee Extensor", "Knee Flexor", "Knee Flexor"];
muscle_prasi = ["Vastus med R", "Rectus femoris R", "Gastrocnemius R", "Biceps femoris R"];
motion = {'Extension','Flexion'};
outliers = {' – No Outlier Detection', ' – Outliers Ignored'};

%% Plotting data 
for o = 1:2
    for ses = 1:5
        figure('Name', "S" + int2str(ses) + outliers(o) ,"Position",[100,1,1200,800]);
        for b = 1:2 
            for sen = 1:4
                if b == 1
                    subplot(2,4,2*sen-1)
                else
                    subplot(2,4,2*sen)
                end
                for con = 1:3
                    % mean and std
                    m = out_exflex_session_sensor_cond_mean_std{1,o}{1,b}{1,ses}{1,sen}{1,con}{1,1};
                    s = out_exflex_session_sensor_cond_mean_std{1,o}{1,b}{1,ses}{1,sen}{1,con}{1,2};
                    t = linspace(0, 100, exflex_med(1,b));
                    
                    % plot std as shaded area
                    x = [t, fliplr(t)];
                    y = [m + s, fliplr(m - s)];
                    fill(x, y, 'o-', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off'); % transparent shading
                    hold on;
                    
                    % plot mean
                    plot(t, m);
                    hold on;
                end
                % label plots
                sgtitle("S" + int2str(ses) +  outliers(o), 'FontSize', 20, 'FontWeight', 'bold') 
                title({"M" + int2str(sen) + " – " + muscle(sen), cell2mat(motion(1,b))})
                xlabel("% cycle")
                ylabel("mV")
                grid on;
                % maimum y-Achse
                ylim([min_signal{1,o}{1,ses}-0.002 max_signal{1,o}{1,ses}+0.005]);
            end
        end
        lg = legend("cond1","cond2","cond3", 'Position', [0.925 0.9 0.05 0.05]); % adjust position
        lg.Title.String = 'Conditions'; % optional: titel for legend
        
    
        % save plots in new folder
        figurename = fullfile(path, '/plots', "S" + int2str(ses) +  outliers(o) + "_som.jpg");
        saveas(gcf, figurename)
    end
end

%% Plot for presentation
% o = 2, ses = 5
figure('Name', "Results" ,"Position",[100,1,1200,800]);
for b = 1:2
    for sen = 1:4
        if b == 1
            subplot(2,4,2*sen-1)
        else
            subplot(2,4,2*sen)
        end
        for con = 1:3
            % mean and std
            m = out_exflex_session_sensor_cond_mean_std{1,2}{1,b}{1,5}{1,sen}{1,con}{1,1};
            s = out_exflex_session_sensor_cond_mean_std{1,2}{1,b}{1,5}{1,sen}{1,con}{1,2};
            t = linspace(0, 100, exflex_med(1,b));
            
            % plot std as shaded area
            x = [t, fliplr(t)];
            y = [m + s, fliplr(m - s)];
            fill(x, y, 'o-', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off'); % transparent shading
            hold on;
            
            % plot mean
            plot(t, m);
            hold on;
        end
        % label plots
        sgtitle("Results", 'FontSize', 20, 'FontWeight', 'bold') 
        title({"M" + int2str(sen) + " – " + muscle_prasi(sen), cell2mat(motion(1,b))})
        xlabel("% cycle")
        ylabel("mV")
        grid on;
        % maimum y-Achse
        ylim([min_signal{1,2}{1,5}-0.002 max_signal{1,2}{1,5}+0.005]);
    end
end
lg = legend("pres1","pres2","pres3", 'Position', [0.925 0.9 0.05 0.05]); % adjust position
lg.Title.String = 'Pressure'; % optional: titel for legend

% save plots in new folder
figurename = fullfile(path, '/plots', "Results.jpg");
saveas(gcf, figurename)

