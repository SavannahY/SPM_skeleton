%% ===========================================================================================
%% plot_ocps - PLOT ELECTRODE OCPs AND CELL OCV WITH STOICHIOMETRY LIMITS
%% ===========================================================================================
%
% DESCRIPTION:
%   Plots the open circuit potentials (OCPs) of both electrodes and the cell OCV,
%   marking the stoichiometric limits (theta_0 and theta_100) for both electrodes.
%
% USAGE:
%   plot_ocps()
%
% OUTPUTS:
%   Creates a figure with three subplots showing:
%   1. Cathode OCP vs stoichiometry with theta_p_0 and theta_p_100 marked
%   2. Anode OCP vs stoichiometry with theta_n_0 and theta_n_100 marked
%   3. Cell OCV vs SOC showing how the electrodes combine

function plot_ocps()

    addpath('cell_model');
    param = modelParameters(101);
    
    % Extract theta values
    theta_n_0 = param.theta_n_0;
    theta_n_100 = param.theta_n_100;
    theta_p_0 = param.theta_p_0;
    theta_p_100 = param.theta_p_100;
    
    theta_p_plot = linspace(0, 1, 1000);
    theta_n_plot = linspace(0, 1, 1000);
    
    % Compute OCPs
    U_p_plot = zeros(size(theta_p_plot));
    U_n_plot = zeros(size(theta_n_plot));
    
    for i = 1:length(theta_p_plot)
        U_p_plot(i) = U_p(theta_p_plot(i), 298.15, param);
    end
    
    for i = 1:length(theta_n_plot)
        U_n_plot(i) = U_n(theta_n_plot(i), 298.15, param);
    end
    
    % Compute OCP values at key theta points
    U_p_0 = U_p(theta_p_0, 298.15, param);
    U_p_100 = U_p(theta_p_100, 298.15, param);
    U_n_0 = U_n(theta_n_0, 298.15, param);
    U_n_100 = U_n(theta_n_100, 298.15, param);
    
    % Create SOC array for cell OCV plot
    soc_array = linspace(0, 1, 1000);
    
    % Map SOC to electrode stoichiometries
    % At SOC = 0: theta_n = theta_n_0, theta_p = theta_p_0
    % At SOC = 1: theta_n = theta_n_100, theta_p = theta_p_100
    theta_n_vs_soc = theta_n_0 + soc_array * (theta_n_100 - theta_n_0);
    theta_p_vs_soc = theta_p_0 + soc_array * (theta_p_100 - theta_p_0);
    
    % Compute cell OCV
    OCV_plot = zeros(size(soc_array));
    for i = 1:length(soc_array)
        U_p_i = U_p(theta_p_vs_soc(i), 298.15, param);
        U_n_i = U_n(theta_n_vs_soc(i), 298.15, param);
        OCV_plot(i) = U_p_i - U_n_i;
    end
    
    % OCV at 0% and 100% SOC
    OCV_0 = U_p_0 - U_n_0;
    OCV_100 = U_p_100 - U_n_100;
    
    % Create figure with three subplots
    figure('Position', [361,214,654,685]);
    
    %% Subplot 1: Cathode OCP
    subplot(3,1,1);
    plot(theta_p_plot, U_p_plot, 'b-', 'LineWidth', 2);
    hold on;
    
    plot(theta_p_0, U_p_0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    plot([theta_p_0, theta_p_0], [min(U_p_plot), U_p_0], 'r--', 'LineWidth', 1.5);
    text(theta_p_0, min(U_p_plot)-0.05, sprintf('$\\theta_p^{0} = %.4f$', theta_p_0), ...
        'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', 'r', 'FontWeight', 'bold', 'Interpreter', 'latex');
    
    plot(theta_p_100, U_p_100, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    plot([theta_p_100, theta_p_100], [min(U_p_plot), U_p_100], 'g--', 'LineWidth', 1.5);
    text(theta_p_100, min(U_p_plot)-0.05, sprintf('$\\theta_p^{100} = %.4f$', theta_p_100), ...
        'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', 'g', 'FontWeight', 'bold', 'Interpreter', 'latex');
    
    grid on;
    xlabel('Cathode Stoichiometry \theta_p [-]', 'FontSize', 12);
    ylabel('U_{p} [V]', 'FontSize', 12);
    title('Cathode OCP', 'FontSize', 14, 'FontWeight', 'bold');
    legend('U_{p}(\theta_p)', '0% SOC', '', '100% SOC', '', 'Location', 'best');
    set(gca, 'FontSize', 11, 'XDir', 'reverse');
    %reverse axis for cathode
    xlim([theta_p_100, theta_p_0]);
    
    %% Subplot 2: Anode OCP
    subplot(3,1,2);
    plot(theta_n_plot, U_n_plot, 'b-', 'LineWidth', 2);
    hold on;
    
    plot(theta_n_0, U_n_0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    plot([theta_n_0, theta_n_0], [min(U_n_plot), U_n_0], 'r--', 'LineWidth', 1.5);
    text(theta_n_0, min(U_n_plot)-0.02, sprintf('$\\theta_n^{0} = %.4f$', theta_n_0), ...
        'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', 'r', 'FontWeight', 'bold', 'Interpreter', 'latex');
    
    plot(theta_n_100, U_n_100, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    plot([theta_n_100, theta_n_100], [min(U_n_plot), U_n_100], 'g--', 'LineWidth', 1.5);
    text(theta_n_100, min(U_n_plot)-0.02, sprintf('$\\theta_n^{100} = %.4f$', theta_n_100), ...
        'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', 'g', 'FontWeight', 'bold', 'Interpreter', 'latex');
    
    grid on;
    xlabel('Anode Stoichiometry \theta_n [-]', 'FontSize', 12);
    ylabel('U_{n} [V]', 'FontSize', 12);
    title('Anode OCP', 'FontSize', 14, 'FontWeight', 'bold');
    legend('U_{n}(\theta_n)', '0% SOC', '', '100% SOC', '', 'Location', 'best');
    set(gca, 'FontSize', 11);
    xlim([theta_n_0, theta_n_100]);
    
    %% Subplot 3: Cell OCV
    subplot(3,1,3);
    plot(soc_array, OCV_plot, 'k-', 'LineWidth', 2);
    hold on;
    
    plot(0, OCV_0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    plot([0, 0], [min(OCV_plot)-0.1, OCV_0], 'r--', 'LineWidth', 1.5);
    text(0, min(OCV_plot)-0.15, sprintf('SOC = 0'), ...
        'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', 'r', 'FontWeight', 'bold');
    
    plot(1, OCV_100, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    plot([1, 1], [min(OCV_plot)-0.1, OCV_100], 'g--', 'LineWidth', 1.5);
    text(1, min(OCV_plot)-0.15, sprintf('SOC = 1'), ...
        'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', 'g', 'FontWeight', 'bold');
    
    grid on;
    xlabel('State of Charge [-]', 'FontSize', 12);
    ylabel('Cell OCV [V]', 'FontSize', 12);
    title('OCV(SOC) = U_{p}(\theta_p(SOC)) - U_{n}(\theta_n(SOC))', ...
        'FontSize', 14, 'FontWeight', 'bold');
    legend('OCV', '0% SOC', '', '100% SOC', '', 'Location', 'best');
    set(gca, 'FontSize', 11);
    xlim([0, 1]);
    
end
