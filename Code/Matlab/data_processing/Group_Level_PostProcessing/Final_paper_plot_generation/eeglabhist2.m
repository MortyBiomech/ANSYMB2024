% EEGLAB history file generated on the 26-Jul-2025
% ------------------------------------------------
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

STUDY = std_topoplot(STUDY,ALLEEG,'clusters',6, 'design', 1, 'plotsubjects', 'on' );
STUDY = std_dipplot(STUDY,ALLEEG,'clusters',6, 'design', 1, 'plotsubjects', 'on' );
STUDY = std_specplot(STUDY,ALLEEG,'clusters',6, 'design', 1, 'plotsubjects', 'on' );
STUDY = pop_specparams(STUDY, 'freqrange',[1 50] );
STUDY = std_specplot(STUDY,ALLEEG,'clusters',6, 'design', 1, 'plotsubjects', 'on' );
STUDY = std_specplot(STUDY,ALLEEG,'clusters',6, 'comps', 8, 'design', 1 );
STUDY = std_specplot(STUDY,ALLEEG,'clusters',6, 'design', 1, 'plotsubjects', 'on' );
STUDY = std_specplot(STUDY,ALLEEG,'clusters',6, 'comps', 6, 'design', 1 );
STUDY = std_specplot(STUDY,ALLEEG,'clusters',6, 'design', 1, 'plotsubjects', 'on' );
STUDY = std_specplot(STUDY,ALLEEG,'clusters',6, 'design', 1, 'plotsubjects', 'on' );
pop_saveh( ALLCOM, 'eeglabhist.m', 'D:\Morteza\MyProjects\ANSYMB2024\Code\Matlab\data_processing\Group_Level_PostProcessing\Final_paper_plot_generation\');

STUDY = std_specplot(STUDY,ALLEEG,'clusters',6, 'design', 1, 'plotsubjects', 'on' );
CURRENTSTUDY = 0;[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, [1:14] ,'retrieve',13,'study',1); 
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 13,'retrieve',8,'study',1); 

figure; 

[ersp, ~, ~, times, freqs, ~, ~, ~, ~] = pop_newtimef( EEG, 0, 7, [0 10000], [4 0.8] , ...
    'topovec', EEG.icawinv(:,7), 'elocs', EEG.chanlocs, ...
    'chaninfo', EEG.chaninfo, 'caption', ['IC 7'], ...
    'baseline',[NaN], 'basenorm', 'off', 'trialbase', 'off', 'freqs', [1 50], 'plotitc' , 'off', ...
    'plotphase', 'off', 'plotphasesign', 'off', 'plotersp', 'off', ...
    'padratio', 2, 'scale', 'abs', 'ntimesout', 800);

ersp2 = ersp - mean(ersp, 2);

eeglab redraw;
