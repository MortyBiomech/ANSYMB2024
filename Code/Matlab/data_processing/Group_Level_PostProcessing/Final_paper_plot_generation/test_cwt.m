clc

%%
f = 500;
freq_limits = [1 70];
voices_per_octave = 32;

t1 = 0:1/f:1;
x1 = 2*sin(2*pi*10*t1);

t2 = 1.002:1/f:3;
x2 = 3*sin(2*pi*14*t2) + [3*cos(2*pi*40*t2(t2 < 2)) , 2*cos(2*pi*55*t2(t2 >= 2))];

t3 = 3.002:1/f:4;
x3 = zeros(1, length(t3));

t4 = 4.002:1/f:5;
x4 = sin(2*pi*5*t4);


ttt = [t1 t2 t3 t4];
xxx = [x1 x2 x3 x4] + 0.3*randn(1, length(ttt));

[cwt_coeffs_test, freqs_test] = cwt(xxx, 'amor', f, ...
        'FrequencyLimits', freq_limits, 'VoicesPerOctave', voices_per_octave);
power_test = abs(cwt_coeffs_test).^2;

figure()
pcolor(ttt, freqs_test, power_test);
shading flat; 
colormap('turbo')
clim([0 max(power_test(:))])
set(gca, 'YScale', 'log');  
yticks([2 4 8 14 20 30 40 50])

hold on
yline([5, 10, 14, 40, 55], 'Color', 'w')


%% 

ttt2 = [ttt, ttt+5.002];
xxx2 = [xxx, xxx];

[cwt_coeffs_test2, freqs_test2] = cwt(xxx2, 'amor', f, ...
        'FrequencyLimits', freq_limits, 'VoicesPerOctave', voices_per_octave);
power_test2 = abs(cwt_coeffs_test2).^2;


figure()

plot(freqs_test, mean(power_test, 2), 'LineStyle','--', 'LineWidth', 2)
hold on
plot(freqs_test2, mean(power_test2, 2))
