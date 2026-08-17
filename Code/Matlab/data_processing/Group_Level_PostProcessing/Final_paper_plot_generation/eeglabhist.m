% EEGLAB history file generated on the 25-Jul-2025
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
eeglab redraw;
