% EEGLAB history file generated on the 25-Aug-2025
% ------------------------------------------------
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

STUDY = std_topoplot(STUDY,ALLEEG,'clusters',[2   3   4   5   6   7   8   9  10  11], 'design', 1);
pop_saveh( ALLCOM, 'eeglabhist3.m', 'D:\Morteza\MyProjects\ANSYMB2024\Code\Matlab\data_processing\Group_Level_PostProcessing\Final_paper_plot_generation\');

STUDY = std_dipplot(STUDY,ALLEEG,'clusters',[2   3   4   5   6   7   8   9  10  11], 'design', 1);
pop_saveh( ALLCOM, 'eeglabhist3.m', 'D:\Morteza\MyProjects\ANSYMB2024\Code\Matlab\data_processing\Group_Level_PostProcessing\Final_paper_plot_generation\');

STUDY = std_dipplot(STUDY,ALLEEG,'clusters',8, 'design', 1);
STUDY = pop_dipparams(STUDY, 'centrline','off');
STUDY = std_dipplot(STUDY,ALLEEG,'clusters',8, 'design', 1);
STUDY = pop_dipparams(STUDY, 'density','on');
STUDY = std_dipplot(STUDY,ALLEEG,'clusters',8, 'design', 1);
STUDY = pop_dipparams(STUDY, 'projlines','on','centrline','on','density','off');
STUDY = pop_dipparams(STUDY, 'projlines','off','centrline','off');
STUDY = std_dipplot(STUDY,ALLEEG,'clusters',8, 'design', 1);
pop_saveh( ALLCOM, 'eeglabhist3.m', 'D:\Morteza\MyProjects\ANSYMB2024\Code\Matlab\data_processing\Group_Level_PostProcessing\Final_paper_plot_generation\');
CURRENTSTUDY = 0;[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, [1:14] ,'retrieve',1,'study',1); 
pop_dipplot( EEG, [1:54] ,'mri','standard_mri.mat','normlen','on');
pop_dipplot( EEG, [1:54] ,'mri','standard_mri.mat');
eeglab redraw;
