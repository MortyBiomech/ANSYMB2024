function plot_muscle_activity_per_Pressure(X_axis_cycle, end_of_flexion_indx, ...
    Vastus_med_R, Rectus_femoris_R, Gastrocnemius_R, Biceps_femoris_R, subject)
    
    P1_mean_color = [0 0.4470 0.7410];
    P1_shade_color = [96, 168, 214]/255;
    P3_mean_color = [0.4660 0.6740 0.1880];
    P3_shade_color = [170, 204, 126]/255;
    P6_mean_color = [0.6350 0.0780 0.1840];
    P6_shade_color = [221, 168, 177]/255;


    %% with SEM

    t = tiledlayout(2,2, "TileSpacing", "compact");
    
    nexttile(1)
    hold on
    plot(X_axis_cycle, 1000*mean(Vastus_med_R.P1, 1), 'LineWidth', 2, 'Color', P1_mean_color);
    shaded_y = [mean(Vastus_med_R.P1, 1)+std(Vastus_med_R.P1, 0, 1)/sqrt(size(Vastus_med_R.P1, 1)), ...
        fliplr(mean(Vastus_med_R.P1, 1)-std(Vastus_med_R.P1, 0, 1)/sqrt(size(Vastus_med_R.P1, 1)))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P1_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)
    
    plot(X_axis_cycle, 1000*mean(Vastus_med_R.P3, 1), 'LineWidth', 2, 'Color', P3_mean_color);
    shaded_y = [mean(Vastus_med_R.P3, 1)+std(Vastus_med_R.P3, 0, 1)/sqrt(size(Vastus_med_R.P3, 1)), ...
        fliplr(mean(Vastus_med_R.P3, 1)-std(Vastus_med_R.P3, 0, 1)/sqrt(size(Vastus_med_R.P3, 1)))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P3_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)
    
    plot(X_axis_cycle, 1000*mean(Vastus_med_R.P6, 1), 'LineWidth', 2, 'Color', P6_mean_color);
    shaded_y = [mean(Vastus_med_R.P6, 1)+std(Vastus_med_R.P6, 0, 1)/sqrt(size(Vastus_med_R.P6, 1)), ...
        fliplr(mean(Vastus_med_R.P6, 1)-std(Vastus_med_R.P6, 0, 1)/sqrt(size(Vastus_med_R.P6, 1)))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P6_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)
    
    xline(X_axis_cycle(end_of_flexion_indx), 'LineWidth', 1, 'LineStyle', '--')
    title('Vastus Medialis')
    xlabel('Cycle [%]')
    ylabel('mV')
    set(gca, 'FontSize', 12)


    nexttile(2)
    hold on
    plot(X_axis_cycle, 1000*mean(Rectus_femoris_R.P1, 1), 'LineWidth', 2, 'Color', P1_mean_color);
    shaded_y = [mean(Rectus_femoris_R.P1, 1)+std(Rectus_femoris_R.P1, 0, 1)/sqrt(size(Rectus_femoris_R.P1, 1)), ...
        fliplr(mean(Rectus_femoris_R.P1, 1)-std(Rectus_femoris_R.P1, 0, 1)/sqrt(size(Rectus_femoris_R.P1, 1)))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P1_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)
    
    plot(X_axis_cycle, 1000*mean(Rectus_femoris_R.P3, 1), 'LineWidth', 2, 'Color', P3_mean_color);
    shaded_y = [mean(Rectus_femoris_R.P3, 1)+std(Rectus_femoris_R.P3, 0, 1)/sqrt(size(Rectus_femoris_R.P3, 1)), ...
        fliplr(mean(Rectus_femoris_R.P3, 1)-std(Rectus_femoris_R.P3, 0, 1)/sqrt(size(Rectus_femoris_R.P3, 1)))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P3_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)
    
    plot(X_axis_cycle, 1000*mean(Rectus_femoris_R.P6, 1), 'LineWidth', 2, 'Color', P6_mean_color);
    shaded_y = [mean(Rectus_femoris_R.P6, 1)+std(Rectus_femoris_R.P6, 0, 1)/sqrt(size(Rectus_femoris_R.P6, 1)), ...
        fliplr(mean(Rectus_femoris_R.P6, 1)-std(Rectus_femoris_R.P6, 0, 1)/sqrt(size(Rectus_femoris_R.P6, 1)))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P6_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)
    
    xline(X_axis_cycle(end_of_flexion_indx), 'LineWidth', 1, 'LineStyle', '--')
    title('Rectus Femoris')
    xlabel('Cycle [%]')
    ylabel('mV')
    set(gca, 'FontSize', 12)



    nexttile(3);
    hold on
    plot(X_axis_cycle, 1000*mean(Gastrocnemius_R.P1, 1), 'LineWidth', 2, 'Color', P1_mean_color);
    shaded_y = [mean(Gastrocnemius_R.P1, 1)+std(Gastrocnemius_R.P1, 0, 1)/sqrt(size(Gastrocnemius_R.P1, 1)), ...
        fliplr(mean(Gastrocnemius_R.P1, 1)-std(Gastrocnemius_R.P1, 0, 1)/sqrt(size(Gastrocnemius_R.P1, 1)))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P1_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    plot(X_axis_cycle, 1000*mean(Gastrocnemius_R.P3(a:b, :), 1), 'LineWidth', 2, 'Color', P3_mean_color);
    shaded_y = [mean(Gastrocnemius_R.P3(a:b, :), 1)+std(Gastrocnemius_R.P3, 0, 1)/sqrt(size(Gastrocnemius_R.P1, 1)), ...
        fliplr(mean(Gastrocnemius_R.P3(a:b, :), 1)-std(Gastrocnemius_R.P3, 0, 1)/sqrt(size(Gastrocnemius_R.P1, 1)))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P3_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)
    
    plot(X_axis_cycle, 1000*mean(Gastrocnemius_R.P6(a:b, :), 1), 'LineWidth', 2, 'Color', P6_mean_color);
    shaded_y = [mean(Gastrocnemius_R.P6(a:b, :), 1)+std(Gastrocnemius_R.P6, 0, 1)/sqrt(size(Gastrocnemius_R.P1, 1)), ...
        fliplr(mean(Gastrocnemius_R.P6(a:b, :), 1)-std(Gastrocnemius_R.P6, 0, 1)/sqrt(size(Gastrocnemius_R.P1, 1)))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P6_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)
    
    xline(X_axis_cycle(end_of_flexion_indx), 'LineWidth', 1, 'LineStyle', '--')
    title('Gastrocnemius Medial')
    xlabel('Cycle [%]')
    ylabel('mV')
    set(gca, 'FontSize', 12)
    hold off
    


    nexttile(4)
    hold on
    plot(X_axis_cycle, 1000*mean(Biceps_femoris_R.P1, 1), 'LineWidth', 2, 'Color', P1_mean_color);
    shaded_y = [mean(Biceps_femoris_R.P1, 1)+std(Biceps_femoris_R.P1, 0, 1)/sqrt(size(Biceps_femoris_R.P1, 1)), ...
        fliplr(mean(Biceps_femoris_R.P1, 1)-std(Biceps_femoris_R.P1, 0, 1)/sqrt(size(Biceps_femoris_R.P1, 1)))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P1_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)
    
    plot(X_axis_cycle, 1000*mean(Biceps_femoris_R.P3, 1), 'LineWidth', 2, 'Color', P3_mean_color);
    shaded_y = [mean(Biceps_femoris_R.P3, 1)+std(Biceps_femoris_R.P3, 0, 1)/sqrt(size(Biceps_femoris_R.P3, 1)), ...
        fliplr(mean(Biceps_femoris_R.P3, 1)-std(Biceps_femoris_R.P3, 0, 1)/sqrt(size(Biceps_femoris_R.P3, 1)))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P3_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)
    
    plot(X_axis_cycle, 1000*mean(Biceps_femoris_R.P6, 1), 'LineWidth', 2, 'Color', P6_mean_color);
    shaded_y = [mean(Biceps_femoris_R.P6, 1)+std(Biceps_femoris_R.P6, 0, 1)/sqrt(size(Biceps_femoris_R.P6, 1)), ...
        fliplr(mean(Biceps_femoris_R.P6, 1)-std(Biceps_femoris_R.P6, 0, 1)/sqrt(size(Biceps_femoris_R.P6, 1)))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P6_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)
    
    xline(X_axis_cycle(end_of_flexion_indx), 'LineWidth', 1, 'LineStyle', '--')
    title('Biceps Femoris')
    xlabel('Cycle [%]')
    ylabel('mV')
    set(gca, 'FontSize', 12)
    

    title(t, ['Subject ', num2str(subject), ...
        ' Muscle Activity: Flexion + Extension '], 'FontSize', 16, 'FontWeight', 'bold')

    set(gcf, 'Position', [100, 100, 1200, 700])




    %% with STD
    % figure();
    t = tiledlayout(2,2, "TileSpacing", "compact");

    nexttile
    hold on

    shaded_y = [mean(Vastus_med_R.P1, 1)+std(Vastus_med_R.P1, 0, 1), ...
        fliplr(mean(Vastus_med_R.P1, 1)-std(Vastus_med_R.P1, 0, 1))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P1_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    plot(X_axis_cycle, 1000*mean(Vastus_med_R.P1, 1), 'LineWidth', 2, 'Color', P1_mean_color);

    shaded_y = [mean(Vastus_med_R.P3, 1)+std(Vastus_med_R.P3, 0, 1), ...
        fliplr(mean(Vastus_med_R.P3, 1)-std(Vastus_med_R.P3, 0, 1))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P3_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    plot(X_axis_cycle, 1000*mean(Vastus_med_R.P3, 1), 'LineWidth', 2, 'Color', P3_mean_color);

    shaded_y = [mean(Vastus_med_R.P6, 1)+std(Vastus_med_R.P6, 0, 1), ...
        fliplr(mean(Vastus_med_R.P6, 1)-std(Vastus_med_R.P6, 0, 1))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P6_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    plot(X_axis_cycle, 1000*mean(Vastus_med_R.P6, 1), 'LineWidth', 2, 'Color', P6_mean_color);

    xline(X_axis_cycle(end_of_flexion_indx), 'LineWidth', 1, 'LineStyle', '--')
    title('Vastus Medialis')
    xlabel('Cycle [%]')
    ylabel('mV')
    set(gca, 'FontSize', 12)


    nexttile
    hold on

    shaded_y = [mean(Rectus_femoris_R.P1, 1)+std(Rectus_femoris_R.P1, 0, 1), ...
        fliplr(mean(Rectus_femoris_R.P1, 1)-std(Rectus_femoris_R.P1, 0, 1))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P1_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    plot(X_axis_cycle, 1000*mean(Rectus_femoris_R.P1, 1), 'LineWidth', 2, 'Color', P1_mean_color);

    shaded_y = [mean(Rectus_femoris_R.P3, 1)+std(Rectus_femoris_R.P3, 0, 1), ...
        fliplr(mean(Rectus_femoris_R.P3, 1)-std(Rectus_femoris_R.P3, 0, 1))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P3_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    plot(X_axis_cycle, 1000*mean(Rectus_femoris_R.P3, 1), 'LineWidth', 2, 'Color', P3_mean_color);

    shaded_y = [mean(Rectus_femoris_R.P6, 1)+std(Rectus_femoris_R.P6, 0, 1), ...
        fliplr(mean(Rectus_femoris_R.P6, 1)-std(Rectus_femoris_R.P6, 0, 1))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P6_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    plot(X_axis_cycle, 1000*mean(Rectus_femoris_R.P6, 1), 'LineWidth', 2, 'Color', P6_mean_color);

    xline(X_axis_cycle(end_of_flexion_indx), 'LineWidth', 1, 'LineStyle', '--')
    title('Rectus Femoris')
    xlabel('Cycle [%]')
    ylabel('mV')
    set(gca, 'FontSize', 12)


    nexttile(3);
    hold on
    shaded_y = [mean(Gastrocnemius_R.P1, 1)+std(Gastrocnemius_R.P1, 0, 1), ...
        fliplr(mean(Gastrocnemius_R.P1, 1)-std(Gastrocnemius_R.P1, 0, 1))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P1_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    plot(X_axis_cycle, 1000*mean(Gastrocnemius_R.P1, 1), 'LineWidth', 2, 'Color', P1_mean_color);

    shaded_y = [mean(Gastrocnemius_R.P3, 1)+std(Gastrocnemius_R.P3, 0, 1), ...
        fliplr(mean(Gastrocnemius_R.P3, 1)-std(Gastrocnemius_R.P3, 0, 1))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P3_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    plot(X_axis_cycle, 1000*mean(Gastrocnemius_R.P3, 1), 'LineWidth', 2, 'Color', P3_mean_color);

    shaded_y = [mean(Gastrocnemius_R.P6, 1)+std(Gastrocnemius_R.P6, 0, 1), ...
        fliplr(mean(Gastrocnemius_R.P6, 1)-std(Gastrocnemius_R.P6, 0, 1))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P6_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    plot(X_axis_cycle, 1000*mean(Gastrocnemius_R.P6, 1), 'LineWidth', 2, 'Color', P6_mean_color);

    xline(X_axis_cycle(end_of_flexion_indx), 'LineWidth', 1, 'LineStyle', '--')
    title('Gastrocnemius Medial')
    xlabel('Cycle [%]')
    ylabel('mV')
    set(gca, 'FontSize', 12)


    nexttile
    hold on

    shaded_y = [mean(Biceps_femoris_R.P1, 1)+std(Biceps_femoris_R.P1, 0, 1), ...
        fliplr(mean(Biceps_femoris_R.P1, 1)-std(Biceps_femoris_R.P1, 0, 1))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P1_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    plot(X_axis_cycle, 1000*mean(Biceps_femoris_R.P1, 1), 'LineWidth', 2, 'Color', P1_mean_color);

    shaded_y = [mean(Biceps_femoris_R.P3, 1)+std(Biceps_femoris_R.P3, 0, 1), ...
        fliplr(mean(Biceps_femoris_R.P3, 1)-std(Biceps_femoris_R.P3, 0, 1))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P3_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    plot(X_axis_cycle, 1000*mean(Biceps_femoris_R.P3, 1), 'LineWidth', 2, 'Color', P3_mean_color);

    shaded_y = [mean(Biceps_femoris_R.P6, 1)+std(Biceps_femoris_R.P6, 0, 1), ...
        fliplr(mean(Biceps_femoris_R.P6, 1)-std(Biceps_femoris_R.P6, 0, 1))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, P6_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    plot(X_axis_cycle, 1000*mean(Biceps_femoris_R.P6, 1), 'LineWidth', 2, 'Color', P6_mean_color);

    xline(X_axis_cycle(end_of_flexion_indx), 'LineWidth', 1, 'LineStyle', '--')
    title('Biceps Femoris')
    xlabel('Cycle [%]')
    ylabel('mV')
    set(gca, 'FontSize', 12)
    legend({'', 'P1', '', 'P3', '', 'P6'})

    title(t, ['Subject ', num2str(subject), ...
        ' Muscle Activity: Flexion + Extension '], 'FontSize', 16, 'FontWeight', 'bold')

    set(gcf, 'Position', [100, 100, 1200, 700])


end