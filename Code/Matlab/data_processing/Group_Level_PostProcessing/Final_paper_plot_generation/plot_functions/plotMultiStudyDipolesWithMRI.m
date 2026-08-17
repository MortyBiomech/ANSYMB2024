

%% --- Locate EEGLAB’s standard MRI (MNI avg152) ---
bemdir  = fullfile(fileparts(which('pop_dipfit_settings')), 'standard_BEM');
mriFile = fullfile(bemdir, 'standard_mri.mat');   

coords = [];
all_colors = [];
colors = [0, 134, 103; ...
          82, 48, 192; ...
          192, 60, 66; ...
          202, 182, 0; ...
          83, 0, 143; ...
          0, 218, 234; ...
          0, 76, 145; ...
          223, 128, 0]./255;
for i = 1:8

    coords = cat(1, coords, all_dipoles{i, 1});
    all_colors = cat(1, all_colors, repmat(colors(i, :), length(all_dipoles{i, 1}), 1));

end


% --- Build a minimal DIPFIT model array for dipplot ---
N = size(coords,1);
model = repmat(struct('posxyz',[],'momxyz',[0 0 1],'rv',0,'active',1), N, 1);
for k = 1:N
    model(k).posxyz = coords(k,:);   % MNI coords (mm)
    model(k).color = all_colors(i, :);
end

all_colors = mat2cell(all_colors, ones(103,1), 3);

%% --- Plot in three views using tiledlayout ---
fig = figure('Color','w'); 
tl = tiledlayout(fig,1,2); %,'Padding','compact','TileSpacing','compact');
views = { [90 90], [90 0]}; 
titles = {'Sagittal (R↔L)','Coronal (A↔P)','Axial (S↔I)'};

axs = gobjects(1,2);
for i = 1:2
    axs(i) = nexttile(tl,i);
    axes(axs(i)); %#ok<LAXES>  % make ax current; dipplot draws into current axes
    dipplot(model, ...
        'coordformat','mni', ...
        'mri', mriFile, ...
        'spheres','on', ...
        'projimg','off', ...     % show MRI slice under projection
        'gui','off', ...
        'color', all_colors, ...
        'cornermri', 'on', ...
        'drawedges', 'off', ...
        'projimg', 'off', ...
        'dipolesize', 30, ...
        'view', views{i});
    

    % ------ make background not black ------
    set(axs(i),'Color','w');               % axes background
    set(ancestor(axs(i),'figure'),'Color','w');           % figure background white

    % Make zero-intensity voxels transparent (standard_mri has zeros outside brain)
    imgs = findobj(axs,'Type','image');
    for im = reshape(imgs,1,[])
        I = get(im,'CData');
        set(im,'AlphaData', double(I>0));   % transparent where MRI==0
    end

    % title(titles{i});
    axis(axs(i),'image'); 
    axis(axs(i),'off');  % neat framing

    camva(axs(i),5);                          % same angle for all axes
    set(axs(i),'CameraViewAngleMode','manual');

end

% camlight(axs(1), 10, 10);
camlight(axs(2), 40, 50);

% keep the layout uniform
tl.TileSpacing = 'loose';
tl.Padding     = 'compact';
set(gcf, 'Position', [140 250 1000 420])


%% plot the cluster mean topoplots


maplimits  = 'absmax';                  
electrodes = 'on';                    

% Use the first dataset’s chanlocs as the (common) montage
T = ALLEEG(1).chanlocs;

new_sets = sets;
new_sets(5) = [];

for i = 1:numel(new_sets)

    subjects = new_sets{1, i}(:, 1) - 4;
    ics = new_sets{1, i}(:, 2);

    % Collect maps (each column is one IC map)
    maps = [];

    for j = 1:numel(subjects)
        m = ALLEEG(subjects(j)).icawinv(:, ics(j));
        maps = cat(2, maps, m);
    end

    % --- Fix sign ambiguity (align polarities to the first map)
    ref = maps(:,1);
    for k = 2:size(maps,2)
        c = corr(ref, maps(:,k), 'rows','pairwise');
        if ~isnan(c) && c < 0 
            maps(:,k) = -maps(:,k); 
        end
    end

    % --- Average across IC maps (ignore NaNs if any)
    meanMap = mean(maps, 2, 'omitmissing');
    
    % --- Plot
    figure('Color','w');
    topoplot(meanMap, T, 'maplimits', maplimits, 'electrodes', electrodes);
    colormap(calldefinedcolormap())

end








