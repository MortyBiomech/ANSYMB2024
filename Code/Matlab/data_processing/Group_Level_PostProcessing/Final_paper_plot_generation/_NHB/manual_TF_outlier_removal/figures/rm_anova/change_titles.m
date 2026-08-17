
savePath = ['D:\Morteza\MyProjects\ANSYMB2024\Code\Matlab', ...
    '\data_processing\Group_Level_PostProcessing\', ...
    'Final_paper_plot_generation\', ...
    '_NHB\manual_TF_outlier_removal\figures\rm_anova\'];

studyNames = {'Left Dorsal ACC', 'Left Parieto Occipital', ...
    'Left PreMot SuppMot', 'Left Prim Motor', 'Prime Visual', ...
    'Right Parieto Occipital', 'Right PreMot SuppMot', 'Right Prim Motor'};

monitors = get(0, 'MonitorPositions');
    
numCond = 3;
fig_width = 1.75*(numCond+1); 
fig_height = 2*fig_width/2.857;

for study = 1:length(studyNames)
    
    figureName = [studyNames{study}, '.fig'];
    openfig(figureName, 'visible');
    set(gcf, 'Position', ...
        [monitors(1,1)-200 monitors(1,2)+600 150*fig_width 150*fig_height]);

    mainT = sgtitle(studyNames{study});
    mainT.FontSize = 20;
    mainT.FontWeight = "bold";
    mainT.FontName = 'Arial';

    savethisfig(gcf, strcat(studyNames{study},'.png'), savePath,'png')
    savethisfig(gcf, strcat(studyNames{study},'.fig'), savePath,'fig')
    savethisfig(gcf, strcat(studyNames{study},'.svg'), savePath,'svg')

    close;

end

