function plot_results(sigmas, results, mode)
%PLOT_RESULTS Plot rotation, translation errors and execution times for
%Manopt and Procrustes in order to compare results graphically
%Results are passed to this function as members the results struct


if (mode == "manopt_sep_vs_procrustes")
    %manopt and procrustes together on the same graph (easier to compare)

    figure("Name", "rot errors"); %figure 1
    plot (sigmas, results.manopt_sep_rot_errs, 'r.', ...
        "DisplayName", "manopt\_sep mean rot error", 'markersize', 15);
    hold on
    xlabel('sigma')
    % ylabel('[°]')
    plot (sigmas, results.procrustes_rot_errs, 'bs', ...
        "DisplayName", "procrustes mean rot error", 'markersize', 15)
    legend;
    hold off
    
    figure("Name", "transl errors"); %figure 2
    plot (sigmas, results.manopt_sep_transl_errs, 'r.', ...
        "DisplayName", "manopt\_sep mean transl error", 'markersize', 10)
    hold on
    xlabel('sigma')
    % ylabel('[m]')
    plot (sigmas, results.procrustes_transl_errs, 'bs', ...
        "DisplayName", "procrustes mean transl error", 'markersize', 10)
    legend;
    hold off
    
    figure("Name", "execution times"); %figure 3
    plot (sigmas, results.manopt_sep_exec_times, 'r.', ...
        "DisplayName", "manopt\_sep mean exec time", 'markersize', 15)
    hold on
    xlabel('sigma')
    ylabel('[s]')
    plot (sigmas, results.procrustes_exec_times, 'bs', ...
        "DisplayName", "procrustes mean exec time", 'markersize', 10)
    legend
    hold off
elseif (mode == "manopt_sep_vs_gen")
    %manopt and procrustes together on the same graph (easier to compare)

    figure("Name", "rot errors"); %figure 1
    plot (sigmas, results.manopt_sep_rot_errs, 'r.', ...
        "DisplayName", "manopt\_sep mean rot error", 'markersize', 15);
    hold on
    xlabel('sigma')
    % ylabel('[°]')
    plot (sigmas, results.manopt_gen_rot_errs, 'bs', ...
        "DisplayName", "manopt\_gen mean rot error", 'markersize', 15)
    legend;
    hold off
    
    figure("Name", "transl errors"); %figure 2
    plot (sigmas, results.manopt_sep_transl_errs, 'r.', ...
        "DisplayName", "manopt\_sep mean transl error", 'markersize', 10)
    hold on
    xlabel('sigma')
    % ylabel('[m]')
    plot (sigmas, results.manopt_gen_transl_errs, 'bs', ...
        "DisplayName", "manopt\_gen mean transl error", 'markersize', 10)
    legend;
    hold off
    
    figure("Name", "execution times"); %figure 3
    plot (sigmas, results.manopt_sep_exec_times, 'r.', ...
        "DisplayName", "manopt\_sep mean exec time", 'markersize', 15)
    hold on
    xlabel('sigma')
    ylabel('[s]')
    plot (sigmas, results.manopt_gen_exec_times, 'bs', ...
        "DisplayName", "manopt\_gen mean exec time", 'markersize', 10)
    legend
    hold off
elseif (mode == "rgrad")
    %manopt and manopt_ad together on the same graph (easier to compare)

    figure("Name", "rot errors"); %figure 1
    plot (sigmas, results.manopt_rot_errs, 'r.', ...
        "DisplayName", "manopt mean rot error", 'markersize', 15);
    hold on
    xlabel('sigma')
    % ylabel('[°]')
    plot (sigmas, results.manopt_ad_rot_errs, 'bs', ...
        "DisplayName", "manopt\_ad mean rot error", 'markersize', 15)
    legend;
    hold off
    
    figure("Name", "transl errors"); %figure 2
    plot (sigmas, results.manopt_transl_errs, 'r.', ...
        "DisplayName", "manopt mean transl error", 'markersize', 10)
    hold on
    xlabel('sigma')
    % ylabel('[m]')
    plot (sigmas, results.manopt_ad_transl_errs, 'bs', ...
        "DisplayName", "manopt\_ad mean transl error", 'markersize', 10)
    legend;
    hold off
    
    figure("Name", "execution times"); %figure 3
    plot (sigmas, results.manopt_exec_times, 'r.', ...
        "DisplayName", "manopt mean exec time", 'markersize', 15)
    hold on
    xlabel('sigma')
    ylabel('[s]')
    plot (sigmas, results.manopt_ad_exec_times, 'bs', ...
        "DisplayName", "manopt\_ad mean exec time", 'markersize', 10)
    legend
    hold off
elseif (mode == "manopt_icp_maxiter")
    %compare (supposed) improvements on Manopt pipeline results when
    %increasing max number of ICP iterations, which should come at the 
    %cost of an increased execution time

    colors = ["r."; "g."; "b."];
    fprintf("maxiter 5 -> red\nmaxiter 10 -> green\nmaxiter 15 -> blue\n");

    disp(size(results));

    figure("Name", "rot errors"); %figure 1
    for ii = 1:size(results, 2)
        result_i = results(ii);
        
        plot (sigmas, result_i.manopt_rot_errs, colors(ii), ...
            "DisplayName", "manopt mean rot error", 'markersize', 15);
        hold on
        xlabel('sigma')
        % ylabel('[°]')
    
        legend;
    end
        hold off
    
    figure("Name", "transl errors"); %figure 2
    for ii = 1:size(results, 2)
        result_i = results(ii);

        plot (sigmas, result_i.manopt_transl_errs, colors(ii), ...
            "DisplayName", "manopt mean transl error", 'markersize', 10)
        hold on
        xlabel('sigma')
        % ylabel('[m]')
        
        legend;
    end
        hold off
        
    figure("Name", "execution times"); %figure 3
    for ii = 1:size(results, 2)
        result_i = results(ii);
        
        plot (sigmas, result_i.manopt_exec_times, colors(ii), ...
            "DisplayName", "manopt mean exec time", 'markersize', 15)
        hold on
        xlabel('sigma')
        ylabel('[s]')
        
        legend
    end
        hold off

elseif (mode == "cp2m_noiseinit")
    %manopt and procrustes together on the same graph (easier to compare)

    figure("Name", "rot errors"); %figure 1
    plot (sigmas, results.manopt_sep_rot_errs, 'r.', ...
        "DisplayName", "manopt\_sep mean rot error", 'markersize', 15);
    hold on
    xlabel('sigma\_init')
    % ylabel('[°]')
    plot (sigmas, results.manopt_gen_rot_errs, 'bs', ...
        "DisplayName", "manopt\_gen mean rot error", 'markersize', 15)
    plot (sigmas, results.procrustes_rot_errs, 'g+', ...
        "DisplayName", "procrustes mean rot error", 'markersize', 15)
    plot (sigmas, results.som_riemstair_rot_errs, 'g+', ...
        "DisplayName", "som_riemstair mean rot error", 'markersize', 15)
    legend;
    hold off
    
    figure("Name", "transl errors"); %figure 2
    plot (sigmas, results.manopt_sep_transl_errs, 'r.', ...
        "DisplayName", "manopt\_sep mean transl error", 'markersize', 10)
    hold on
    xlabel('sigma\_init')
    % ylabel('[m]')
    plot (sigmas, results.manopt_gen_transl_errs, 'bs', ...
        "DisplayName", "manopt\_gen mean transl error", 'markersize', 10)
    plot (sigmas, results.procrustes_transl_errs, 'g+', ...
        "DisplayName", "procrustes mean transl error", 'markersize', 15)
    plot (sigmas, results.som_riemstair_transl_errs, 'g+', ...
        "DisplayName", "som_riemstair mean transl error", 'markersize', 15)
    legend;
    hold off
    
    figure("Name", "execution times"); %figure 3
    plot (sigmas, results.manopt_sep_exec_times, 'r.', ...
        "DisplayName", "manopt\_sep mean exec time", 'markersize', 15)
    hold on
    xlabel('sigma\_init')
    ylabel('[s]')
    plot (sigmas, results.manopt_gen_exec_times, 'bs', ...
        "DisplayName", "manopt\_gen mean exec time", 'markersize', 10)
    plot (sigmas, results.procrustes_exec_times, 'g+', ...
        "DisplayName", "procrustes mean exec time", 'markersize', 15)
    plot (sigmas, results.som_riemstair_exec_times, 'g+', ...
        "DisplayName", "som_riemstair mean exec time", 'markersize', 15)
    legend
    hold off

elseif (mode == "procrustes_manopt_riemstaircase")
    %manopt and procrustes together on the same graph (easier to compare)

    figure("Name", "rot errors"); %figure 1
    plot (sigmas, results.manopt_sep_rot_errs, 'r.', ...
        "DisplayName", "manopt\_sep mean rot error", 'markersize', 15);
    hold on
    xlabel('sigma')
    % ylabel('[°]')
    plot (sigmas, results.procrustes_rot_errs, 'bs', ...
        "DisplayName", "procrustes mean rot error", 'markersize', 15)
    plot (sigmas, results.manopt_rc_rot_errs, 'g+', ...
        "DisplayName", "manopt\_rc mean rot error", 'markersize', 15);
    legend;
    hold off
    
    figure("Name", "transl errors"); %figure 2
    plot (sigmas, results.manopt_sep_transl_errs, 'r.', ...
        "DisplayName", "manopt\_sep mean transl error", 'markersize', 10)
    hold on
    xlabel('sigma')
    % ylabel('[m]')
    plot (sigmas, results.procrustes_transl_errs, 'bs', ...
        "DisplayName", "procrustes mean transl error", 'markersize', 10)
    plot (sigmas, results.manopt_rc_transl_errs, 'g+', ...
        "DisplayName", "manopt\_rc mean transl error", 'markersize', 10)
    legend;
    hold off
    
    figure("Name", "execution times"); %figure 3
    plot (sigmas, results.manopt_sep_exec_times, 'r.', ...
        "DisplayName", "manopt\_sep mean exec time", 'markersize', 15)
    hold on
    xlabel('sigma')
    ylabel('[s]')
    plot (sigmas, results.procrustes_exec_times, 'bs', ...
        "DisplayName", "procrustes mean exec time", 'markersize', 10)
    plot (sigmas, results.manopt_rc_exec_times, 'g+', ...
        "DisplayName", "manopt\_rc mean exec time", 'markersize', 15)
    legend
    hold off
elseif(mode == "procrustes_manopt_sesyncriemstair")
    %manopt and procrustes together on the same graph (easier to compare)

    figure("Name", "rot errors"); %figure 1
    plot (sigmas, results.manopt_sep_rot_errs, 'r.', ...
        "DisplayName", "manopt\_sep mean rot error", 'markersize', 15);
    hold on
    xlabel('sigma')
    % ylabel('[°]')
    plot (sigmas, results.procrustes_rot_errs, 'bs', ...
        "DisplayName", "procrustes mean rot error", 'markersize', 15)
    plot (sigmas, results.sesync_riemstair_rot_errs, 'g+', ...
        "DisplayName", "sesync\_riemstair mean rot error", 'markersize', 15);
    legend;
    hold off
    
    figure("Name", "transl errors"); %figure 2
    plot (sigmas, results.manopt_sep_transl_errs, 'r.', ...
        "DisplayName", "manopt\_sep mean transl error", 'markersize', 10)
    hold on
    xlabel('sigma')
    % ylabel('[m]')
    plot (sigmas, results.procrustes_transl_errs, 'bs', ...
        "DisplayName", "procrustes mean transl error", 'markersize', 10)
    plot (sigmas, results.sesync_riemstair_transl_errs, 'g+', ...
        "DisplayName", "sesync\_riemstair mean transl error", 'markersize', 10)
    legend;
    hold off
    
    figure("Name", "execution times"); %figure 3
    plot (sigmas, results.manopt_sep_exec_times, 'r.', ...
        "DisplayName", "manopt\_sep mean exec time", 'markersize', 15)
    hold on
    xlabel('sigma')
    ylabel('[s]')
    plot (sigmas, results.procrustes_exec_times, 'bs', ...
        "DisplayName", "procrustes mean exec time", 'markersize', 10)
    plot (sigmas, results.sesync_riemstair_exec_times, 'g+', ...
        "DisplayName", "sesync\_riemstair mean exec time", 'markersize', 15)
    legend
    hold off
else
    disp("Plot mode unknown!");
end


end %function


