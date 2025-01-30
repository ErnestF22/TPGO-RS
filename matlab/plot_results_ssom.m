function plot_results_ssom(sigmas, results, mode)
%PLOT_RESULTS Plot rotation, translation errors and execution times for
%Manopt and Procrustes in order to compare results graphically
%Results are passed to this function as members the results struct


if (mode == "ssom")
    %manopt and procrustes together on the same graph (easier to compare)

    figure("Name", "rot errors"); %figure 1
    plot (sigmas, results.ssom_rot_errs, 'rD', ...
        "DisplayName", "ssom mean rot error", 'markersize', 15, ...
        'LineWidth',10);
    hold on
    xlabel('sigma')
    % ylabel('[°]')
    % plot (sigmas, results.procrustes_rot_errs, 'bs', ...
    %     "DisplayName", "procrustes mean rot error", 'markersize', 15, ...
    %     'LineWidth',10)
    % plot (sigmas, results.manopt_rs_rot_errs, 'g+', ...
    %     "DisplayName", "manopt\_rs mean rot error", 'markersize', 15, ...
    %     'LineWidth',15);
    set(legend,'FontSize',30);
    % legend;
    hold off
    
    figure("Name", "transl errors"); %figure 2
    plot (sigmas, results.ssom_transl_errs, 'rD', ...
        "DisplayName", "ssom mean transl error", 'markersize', 10,'LineWidth',10)
    hold on
    xlabel('sigma')
    % ylabel('[m]')
    % plot (sigmas, results.procrustes_transl_errs, 'bs', ...
    %     "DisplayName", "procrustes mean transl error", 'markersize', 10,'LineWidth',10)
    % plot (sigmas, results.manopt_rs_transl_errs, 'g+', ...
    %     "DisplayName", "manopt\_rs mean transl error", 'markersize', 10,'LineWidth',15)
    set(legend,'FontSize',30);
    legend;
    hold off
    
    figure("Name", "execution times"); %figure 3
    plot (sigmas, results.ssom_exec_times, 'rD', ...
        "DisplayName", "ssom mean exec time", 'markersize', 15,'LineWidth',10)
    hold on
    xlabel('sigma')
    % ylabel('[s]')
    % plot (sigmas, results.procrustes_exec_times, 'bs', ...
    %     "DisplayName", "procrustes mean exec time", 'markersize', 10,'LineWidth',10)
    % plot (sigmas, results.manopt_rs_exec_times, 'g+', ...
    %     "DisplayName", "manopt\_rs mean exec time", 'markersize', 15,'LineWidth',15)
    % set(legend,'FontSize',30);
    legend
    hold off

    figure("Name", "scale ratios"); %figure 4
    plot (sigmas, results.ssom_scale_ratios, 'rD', ...
        "DisplayName", "ssom mean scale ratios", 'markersize', 15,'LineWidth',10)
    hold on
    xlabel('sigma')
    % ylabel('[s]')
    legend
    hold off

    figure("Name", "translation errors normalized"); %figure 5
    plot (sigmas, results.ssom_transl_errs_norm, 'rD', ...
        "DisplayName", "ssom mean scale ratios", 'markersize', 15,'LineWidth',10)
    hold on
    xlabel('sigma')
    % ylabel('[s]')
    legend
    hold off
else
    disp("Plot mode unknown!");
end


end %function


