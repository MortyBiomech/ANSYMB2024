function plot_muscle_activity_NoPAM(X_axis_cycle, ...
            Vastus_med_R_ref, Rectus_femoris_R_ref, ...
            Gastrocnemius_R_ref, Biceps_femoris_R_ref)

    No_PAM_color = 0.5*[1, 1, 1];

    %% with SEM
    nexttile(1)
    hold on
    plot(X_axis_cycle, 1000*mean(Vastus_med_R_ref, 1), 'LineWidth', 2, 'Color', 'k');
    shaded_y = [mean(Vastus_med_R_ref, 1)+std(Vastus_med_R_ref, 0, 1)/sqrt(size(Vastus_med_R_ref, 1)), ...
        fliplr(mean(Vastus_med_R_ref, 1)-std(Vastus_med_R_ref, 0, 1)/sqrt(size(Vastus_med_R_ref, 1)))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, No_PAM_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)
    

    nexttile(2)
    hold on
    plot(X_axis_cycle, 1000*mean(Rectus_femoris_R_ref, 1), 'LineWidth', 2, 'Color', 'k');
    shaded_y = [mean(Rectus_femoris_R_ref, 1)+std(Rectus_femoris_R_ref, 0, 1)/sqrt(size(Rectus_femoris_R_ref, 1)), ...
        fliplr(mean(Rectus_femoris_R_ref, 1)-std(Rectus_femoris_R_ref, 0, 1)/sqrt(size(Rectus_femoris_R_ref, 1)))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, No_PAM_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)
    
    
    nexttile(3)
    hold on
    plot(X_axis_cycle, 1000*mean(Gastrocnemius_R_ref, 1), 'LineWidth', 2, 'Color', 'k');
    shaded_y = [mean(Gastrocnemius_R_ref, 1)+std(Gastrocnemius_R_ref, 0, 1)/sqrt(size(Gastrocnemius_R_ref, 1)), ...
        fliplr(mean(Gastrocnemius_R_ref, 1)-std(Gastrocnemius_R_ref, 0, 1)/sqrt(size(Gastrocnemius_R_ref, 1)))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, No_PAM_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)
    
    
    nexttile(4)
    hold on
    plot(X_axis_cycle, 1000*mean(Biceps_femoris_R_ref, 1), 'LineWidth', 2, 'Color', 'k');
    shaded_y = [mean(Biceps_femoris_R_ref, 1)+std(Biceps_femoris_R_ref, 0, 1)/sqrt(size(Biceps_femoris_R_ref, 1)), ...
        fliplr(mean(Biceps_femoris_R_ref, 1)-std(Biceps_femoris_R_ref, 0, 1)/sqrt(size(Biceps_femoris_R_ref, 1)))];
    shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    fill(shaded_x, 1000*shaded_y, No_PAM_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)
    
    
    % %% with STD
    % figure();
    % t = tiledlayout(2,2, "TileSpacing", "compact");
    % 
    % nexttile
    % hold on
    % 
    % shaded_y = [mean(Vastus_med_R, 1)+std(Vastus_med_R, 0, 1), ...
    %     fliplr(mean(Vastus_med_R, 1)-std(Vastus_med_R, 0, 1))];
    % shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    % fill(shaded_x, 1000*shaded_y, P1_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    % plot(X_axis_cycle, 1000*mean(Vastus_med_R, 1), 'LineWidth', 2, 'Color', P1_mean_color);
    % 
    % shaded_y = [mean(Vastus_med_R, 1)+std(Vastus_med_R, 0, 1), ...
    %     fliplr(mean(Vastus_med_R, 1)-std(Vastus_med_R, 0, 1))];
    % shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    % fill(shaded_x, 1000*shaded_y, P3_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    % plot(X_axis_cycle, 1000*mean(Vastus_med_R, 1), 'LineWidth', 2, 'Color', P3_mean_color);
    % 
    % shaded_y = [mean(Vastus_med_R, 1)+std(Vastus_med_R, 0, 1), ...
    %     fliplr(mean(Vastus_med_R, 1)-std(Vastus_med_R, 0, 1))];
    % shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    % fill(shaded_x, 1000*shaded_y, P6_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    % plot(X_axis_cycle, 1000*mean(Vastus_med_R, 1), 'LineWidth', 2, 'Color', P6_mean_color);
    % 
    % xline(X_axis_cycle(end_of_flexion_indx), 'LineWidth', 1, 'LineStyle', '--')
    % title('Vastus Medialis')
    % xlabel('Cycle [%]')
    % ylabel('mV')
    % set(gca, 'FontSize', 12)
    % 
    % 
    % nexttile
    % hold on
    % 
    % shaded_y = [mean(Rectus_femoris_R, 1)+std(Rectus_femoris_R, 0, 1), ...
    %     fliplr(mean(Rectus_femoris_R, 1)-std(Rectus_femoris_R, 0, 1))];
    % shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    % fill(shaded_x, 1000*shaded_y, P1_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    % plot(X_axis_cycle, 1000*mean(Rectus_femoris_R, 1), 'LineWidth', 2, 'Color', P1_mean_color);
    % 
    % shaded_y = [mean(Rectus_femoris_R, 1)+std(Rectus_femoris_R, 0, 1), ...
    %     fliplr(mean(Rectus_femoris_R, 1)-std(Rectus_femoris_R, 0, 1))];
    % shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    % fill(shaded_x, 1000*shaded_y, P3_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    % plot(X_axis_cycle, 1000*mean(Rectus_femoris_R, 1), 'LineWidth', 2, 'Color', P3_mean_color);
    % 
    % shaded_y = [mean(Rectus_femoris_R, 1)+std(Rectus_femoris_R, 0, 1), ...
    %     fliplr(mean(Rectus_femoris_R, 1)-std(Rectus_femoris_R, 0, 1))];
    % shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    % fill(shaded_x, 1000*shaded_y, P6_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    % plot(X_axis_cycle, 1000*mean(Rectus_femoris_R, 1), 'LineWidth', 2, 'Color', P6_mean_color);
    % 
    % xline(X_axis_cycle(end_of_flexion_indx), 'LineWidth', 1, 'LineStyle', '--')
    % title('Rectus Femoris')
    % xlabel('Cycle [%]')
    % ylabel('mV')
    % set(gca, 'FontSize', 12)
    % 
    % 
    % nexttile
    % hold on
    % shaded_y = [mean(Gastrocnemius_R, 1)+std(Gastrocnemius_R, 0, 1), ...
    %     fliplr(mean(Gastrocnemius_R, 1)-std(Gastrocnemius_R, 0, 1))];
    % shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    % fill(shaded_x, 1000*shaded_y, P1_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    % plot(X_axis_cycle, 1000*mean(Gastrocnemius_R, 1), 'LineWidth', 2, 'Color', P1_mean_color);
    % 
    % shaded_y = [mean(Gastrocnemius_R, 1)+std(Gastrocnemius_R, 0, 1), ...
    %     fliplr(mean(Gastrocnemius_R, 1)-std(Gastrocnemius_R, 0, 1))];
    % shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    % fill(shaded_x, 1000*shaded_y, P3_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    % plot(X_axis_cycle, 1000*mean(Gastrocnemius_R, 1), 'LineWidth', 2, 'Color', P3_mean_color);
    % 
    % shaded_y = [mean(Gastrocnemius_R, 1)+std(Gastrocnemius_R, 0, 1), ...
    %     fliplr(mean(Gastrocnemius_R, 1)-std(Gastrocnemius_R, 0, 1))];
    % shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    % fill(shaded_x, 1000*shaded_y, P6_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    % plot(X_axis_cycle, 1000*mean(Gastrocnemius_R, 1), 'LineWidth', 2, 'Color', P6_mean_color);
    % 
    % xline(X_axis_cycle(end_of_flexion_indx), 'LineWidth', 1, 'LineStyle', '--')
    % title('Gastrocnemius Medial')
    % xlabel('Cycle [%]')
    % ylabel('mV')
    % set(gca, 'FontSize', 12)
    % 
    % 
    % nexttile
    % hold on
    % 
    % shaded_y = [mean(Biceps_femoris_R, 1)+std(Biceps_femoris_R, 0, 1), ...
    %     fliplr(mean(Biceps_femoris_R, 1)-std(Biceps_femoris_R, 0, 1))];
    % shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    % fill(shaded_x, 1000*shaded_y, P1_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    % plot(X_axis_cycle, 1000*mean(Biceps_femoris_R, 1), 'LineWidth', 2, 'Color', P1_mean_color);
    % 
    % shaded_y = [mean(Biceps_femoris_R, 1)+std(Biceps_femoris_R, 0, 1), ...
    %     fliplr(mean(Biceps_femoris_R, 1)-std(Biceps_femoris_R, 0, 1))];
    % shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    % fill(shaded_x, 1000*shaded_y, P3_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    % plot(X_axis_cycle, 1000*mean(Biceps_femoris_R, 1), 'LineWidth', 2, 'Color', P3_mean_color);
    % 
    % shaded_y = [mean(Biceps_femoris_R, 1)+std(Biceps_femoris_R, 0, 1), ...
    %     fliplr(mean(Biceps_femoris_R, 1)-std(Biceps_femoris_R, 0, 1))];
    % shaded_x = [X_axis_cycle, fliplr(X_axis_cycle)];
    % fill(shaded_x, 1000*shaded_y, P6_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    % plot(X_axis_cycle, 1000*mean(Biceps_femoris_R, 1), 'LineWidth', 2, 'Color', P6_mean_color);
    % 
    % xline(X_axis_cycle(end_of_flexion_indx), 'LineWidth', 1, 'LineStyle', '--')
    % title('Biceps Femoris')
    % xlabel('Cycle [%]')
    % ylabel('mV')
    % set(gca, 'FontSize', 12)
    % legend({'', 'P1', '', 'P3', '', 'P6'})
    % 
    % title(t, ['Subject ', num2str(subject), ...
    %     ' Muscle Activity: Flexion + Extension '], 'FontSize', 16, 'FontWeight', 'bold')
    % 
    % set(gcf, 'Position', [100, 100, 1200, 700])



end