clc

a1 = ALLEEG(10).event(3).latency;
a2 = ALLEEG(10).event(5).latency;
base = ALLEEG(10).icaact(1, a1:a2);

sig = subjects_data.sub14{1, 1}.EEG_stream.Preprocessed.Sources{1, 1}(1,:);

sig2 = ALLEEG(10).icaact(1, a2:a2+length(sig));

figure()
plot([base sig], 'LineWidth', 3) 
hold on
xline(length(base))

plot(length(base)+1:length(base)+length(sig2), sig2, 'LineStyle', '-',  'Color',    'r')

