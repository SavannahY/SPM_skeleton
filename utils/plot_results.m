function plot_results(all_data)
    %PLOT_RESULTS Visualize SPM simulation results
    %
    % INPUTS:
    %   all_data - Structure containing simulation results with units:
    %              .t           - Time [s]
    %              .V           - Voltage [V]
    %              .soc_bulk_n  - Anode SOC [-]
    %              .soc_bulk_p  - Cathode SOC [-]
    %              .I           - Current [A]
    
    figure();
    
    ax(1) = subplot(3, 1, 1);
    plot(all_data.t, all_data.V, 'b', 'LineWidth', 1.5);
    hold on
    plot(all_data.t, all_data.V_ref, 'r--', 'LineWidth', 1.5);
    legend('Simulated Voltage', 'Experimental Voltage');
    ylabel('Voltage [V]');
    grid minor;
    title('Simulation Results');

    ax(2) = subplot(3, 1, 2); hold on;
    plot(all_data.t, all_data.soc_bulk_n, 'b-', 'LineWidth', 1.5, 'DisplayName', 'SOC_n');
    plot(all_data.t, all_data.soc_bulk_p, 'r--', 'LineWidth', 1.5, 'DisplayName', 'SOC_p');
    ylabel('SOC, [-]');
    grid minor;
    legend()

    ax(3) = subplot(3, 1, 3);
    plot(all_data.t, all_data.I, 'm', 'LineWidth', 1.5);
    ylabel('Current [A]');
    xlabel('Time [s]');
    grid minor;

    linkaxes(ax, 'x');
    xlim([all_data.t(1), all_data.t(end)]);

end