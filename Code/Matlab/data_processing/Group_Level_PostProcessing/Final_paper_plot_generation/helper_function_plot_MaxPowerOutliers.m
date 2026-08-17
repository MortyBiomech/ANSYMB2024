function helper_function_plot_MaxPowerOutliers(subject, ic, max_power, ...
    max_power_P1, max_power_P3, max_power_P6)

    factor = 3;
    outliers = isoutlier(max_power, "median", "ThresholdFactor", factor);
    outliers_P1 = isoutlier(max_power_P1(1, :), "median", "ThresholdFactor", factor);
    outliers_P3 = isoutlier(max_power_P3(1, :), "median", "ThresholdFactor", factor);
    outliers_P6 = isoutlier(max_power_P6(1, :), "median", "ThresholdFactor", factor);

    figure()
    tiledlayout(2, 3)
    nexttile(1, [1 3])
    plot(max_power, 'LineWidth', 2); hold on
    x = 1:length(max_power);
    plot(x(outliers), max_power(outliers), 'Marker', 'o', 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', 'r')
    title(['Subject ', num2str(subject), ' IC ', num2str(ic), ', '...
        num2str(length(max_power)), ' trials, ', ...
        num2str(sum(outliers)), ' outliers'], 'FontSize', 14, 'FontWeight', 'normal')
    ylabel('Max Power (\muV^2)', 'FontSize', 14)
    xlabel('Trials', 'FontSize', 14)
    set(gca, 'FontSize', 14)

    nexttile(4)
    plot(max_power_P1(1, :), 'LineWidth', 2); hold on
    x = 1:length(max_power_P1);
    plot(x(outliers_P1), max_power_P1(1, outliers_P1), 'Marker', 'o','Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', 'r')
    title(['P1, ', num2str(length(max_power_P1)), ' trials, ', ...
        num2str(sum(outliers_P1)), ' outliers'], 'FontSize', 14, 'FontWeight', 'normal')
    ylabel('Max Power (\muV^2)', 'FontSize', 14)
    xlabel('Trials', 'FontSize', 14)
    set(gca, 'FontSize', 14)
    
    nexttile(5)
    plot(max_power_P3(1, :), 'LineWidth', 2); hold on
    x = 1:length(max_power_P3);
    plot(x(outliers_P3), max_power_P3(1, outliers_P3), 'Marker', 'o', 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', 'r')
    title(['P3, ', num2str(length(max_power_P3)), ' trials, ', ...
        num2str(sum(outliers_P3)), ' outliers'], 'FontSize', 14, 'FontWeight', 'normal')
    xlabel('Trials', 'FontSize', 14)
    set(gca, 'FontSize', 14)
    
    nexttile(6)
    plot(max_power_P6(1, :), 'LineWidth', 2); hold on
    x = 1:length(max_power_P6);
    plot(x(outliers_P6), max_power_P6(1, outliers_P6), 'Marker', 'o', 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', 'r')
    title(['P6, ', num2str(length(max_power_P6)), ' trials, ', ...
        num2str(sum(outliers_P6)), ' outliers'], 'FontSize', 14, 'FontWeight', 'normal')
    xlabel('Trials', 'FontSize', 14)
    set(gca, 'FontSize', 14)

end