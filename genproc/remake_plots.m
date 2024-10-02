function remake_plots

close all;

ws_filename = "_n9_mindeg2.mat";

load(ws_filename, "results");
load(ws_filename, "sigmas");
load(ws_filename, "test_str");
sz = [4 3];

for ii = [1]
    sigmas(end) = [];
    results.manopt_sep_rot_errs(end) = [];
    results.procrustes_rot_errs(end) = [];
    results.manopt_rs_rot_errs(end) = [];
    
    results.manopt_sep_transl_errs(end) = [];
    results.procrustes_transl_errs(end) = [];
    results.manopt_rs_transl_errs(end) = [];
    
    results.manopt_sep_exec_times(end) = [];
    results.procrustes_exec_times(end) = [];
    results.manopt_rs_exec_times(end) = [];
end

close all;

figure("Name", "rot errors"); %figure 1
plot (sigmas, results.manopt_sep_rot_errs, 'Red.', ...
    "DisplayName", "SOM-ICP", 'markersize', 15);
hold on
xlabel('sigma')
ylabel('mean rotation error')
plot (sigmas, results.procrustes_rot_errs, 'Blues', ...
    "DisplayName", "SOM-PROCR", 'markersize', 10)
plot (sigmas, results.manopt_rs_rot_errs, 'Green+', ...
    "DisplayName", "SOM-RS", 'markersize', 15);
% legend('Location','eastoutside');
lgd = legend;
fontsize(lgd,6,'points')
rot_fig_name = convertStringsToChars(strcat("rot_errors", test_str));
savefigure(rot_fig_name,'epsc',sz)
hold off

figure("Name", "transl errors"); %figure 2
plot (sigmas, results.manopt_sep_transl_errs, 'Red.', ...
    "DisplayName", "SOM-ICP", 'markersize', 15)
hold on
xlabel('sigma')
ylabel('mean translation error')
plot (sigmas, results.procrustes_transl_errs, 'Blues', ...
    "DisplayName", "SOM-PROCR", 'markersize', 10)
plot (sigmas, results.manopt_rs_transl_errs, 'Green+', ...
    "DisplayName", "SOM-RS", 'markersize', 15)
lgd = legend;
fontsize(lgd,6,'points')
transl_fig_name = convertStringsToChars(strcat('transl_errors', test_str));
savefigure(transl_fig_name,'epsc',sz)
hold off

figure("Name", "execution times"); %figure 3
plot (sigmas, results.manopt_sep_exec_times, 'Red.', ...
    "DisplayName", "SOM-ICP mean exec time", 'markersize', 15)
hold on
xlabel('sigma')
ylabel('[s]')
plot (sigmas, results.procrustes_exec_times, 'Blues', ...
    "DisplayName", "SOM-PROCR mean exec time", 'markersize', 10)
plot (sigmas, results.manopt_rs_exec_times, 'Green+', ...
    "DisplayName", "SOM-RS mean exec time", 'markersize', 15)
lgd = legend;
fontsize(lgd,6,'points')
exectimes_fig_name = convertStringsToChars(strcat('exec_times', test_str));
savefigure(exectimes_fig_name,'epsc',sz)
hold off

%%

y_transl_icp_mindeg3 = zeros(5,4);
y_transl_procrustes_mindeg3 = zeros(5,4);
y_transl_rs_mindeg3 = zeros(5,4);

y_rot_icp_mindeg3 = zeros(5,4);
y_rot_procrustes_mindeg3 = zeros(5,4);
y_rot_rs_mindeg3 = zeros(5,4);

y_rot_icp_mindeg2 = zeros(5,4);
y_rot_procrustes_mindeg2 = zeros(5,4);
y_rot_rs_mindeg2 = zeros(5,4);

y_transl_icp_mindeg2 = zeros(5,4);
y_transl_procrustes_mindeg2 = zeros(5,4);
y_transl_rs_mindeg2 = zeros(5,4);


ns = [5 6 7 8 9];
for n = ns
    idx = n-4;
    % mindeg 3
    mindeg = 3;
    ws_filename = strcat("_n", string(n), "_mindeg", string(mindeg), ".mat");
    load(ws_filename, "results");
%     load(ws_filename, "sigmas");
%     load(ws_filename, "test_str");    
    
    y_rot_icp_mindeg3(:,idx) =  results.manopt_sep_rot_errs(1:end-1);
    y_rot_procrustes_mindeg3(:,idx) = results.procrustes_rot_errs(1:end-1);
    y_rot_rs_mindeg3(:,idx) = results.manopt_rs_rot_errs(1:end-1);
    
    y_transl_icp_mindeg3(:,idx) = results.manopt_sep_transl_errs(1:end-1);
    y_transl_procrustes_mindeg3(:,idx) = results.procrustes_transl_errs(1:end-1);
    y_transl_rs_mindeg3(:,idx) = results.manopt_rs_transl_errs(1:end-1);

    % mindeg 2
    mindeg = 2;
    ws_filename = strcat("_n", string(n), "_mindeg", string(mindeg), ".mat");
    load(ws_filename, "results");
%     load(ws_filename, "sigmas");
%     load(ws_filename, "test_str");

    
    y_rot_icp_mindeg2(:,idx) = results.manopt_sep_rot_errs(1:end-1);
    y_rot_procrustes_mindeg2(:,idx) = results.procrustes_rot_errs(1:end-1);
    y_rot_rs_mindeg2(:,idx) = results.manopt_rs_rot_errs(1:end-1);

    y_transl_icp_mindeg2(:,idx) = results.manopt_sep_transl_errs(1:end-1);
    y_transl_procrustes_mindeg2(:,idx) = results.procrustes_transl_errs(1:end-1);
    y_transl_rs_mindeg2(:,idx) = results.manopt_rs_transl_errs(1:end-1);

end

sigmas(end) = [];

for s = 1:length(sigmas)
    sigma = sigmas(s);

    figure("Name", strcat("Rotation error progression for sigma = ", string(sigma), ""));
    plot (ns, y_rot_icp_mindeg3(s,:), 'o', 'Color',[1 0 0], ...
        "DisplayName", "SOM-ICP mindeg 3", 'markersize', 11);
    hold on
    xlabel('n')
    ylabel('mean rotation error')
    title_text = strcat("sigma ", string(sigma));
    title(title_text)
    plot (ns, y_rot_icp_mindeg2(s,:), 'o', 'Color',[.6 0 0], ...
        "DisplayName", "SOM-ICP mindeg 2", 'markersize', 11);
    plot (ns, y_rot_procrustes_mindeg3(s,:), 's', 'Color',[0 1 0], ...
        "DisplayName", "SOM-PROCR mindeg 3", 'markersize', 11)
    plot (ns, y_rot_procrustes_mindeg2(s,:), 's', 'Color',[0 0.6 0], ...
        "DisplayName", "SOM-PROCR mindeg 2", 'markersize', 11)
    plot (ns, y_rot_rs_mindeg3(s,:), '^', 'Color',[0 0 1], ...
        "DisplayName", "SOM-RS mindeg 3", 'markersize', 11);
    plot (ns, y_rot_rs_mindeg2(s,:), '^', 'Color',[0 0 0.6], ...
        "DisplayName", "SOM-RS mindeg 2", 'markersize', 11);
    lgd = legend;
    fontsize(lgd,6,'points')
    rot_fig_name = convertStringsToChars(strcat("rot_err_progression_sigma", string(sigma)));
    rot_fig_name(rot_fig_name=='.')=[];
    savefigure(rot_fig_name,'epsc',sz)
    hold off
    
    figure("Name", strcat("Translation error progression for sigma = ", string(sigma), ""));
    plot (ns, y_transl_icp_mindeg3(s,:), 'o', 'Color',[1 0 0], ...
        "DisplayName", "SOM-ICP mindeg 3", 'markersize', 11)
    hold on
    xlabel('n')
    ylabel('mean translation error')
    title_text = strcat("sigma ", string(sigma));
    title(title_text)
    plot (ns, y_transl_icp_mindeg2(s,:), 'o', 'Color',[.6 0 0], ...
        "DisplayName", "SOM-ICP mindeg 2", 'markersize', 11);
    plot (ns, y_transl_procrustes_mindeg3(s,:), 's', 'Color',[0 1 0], ...
        "DisplayName", "SOM-PROCR mindeg 3", 'markersize', 11)
    plot (ns, y_transl_procrustes_mindeg2(s,:), 's', 'Color',[0 0.6 0], ...
        "DisplayName", "SOM-PROCR mindeg 2", 'markersize', 11)
    plot (ns, y_transl_rs_mindeg3(s,:), '^', 'Color', [0 0 1], ...
        "DisplayName", "SOM-RS mindeg 3", 'markersize', 11);
    plot (ns, y_transl_rs_mindeg2(s,:), '^', 'Color',[0 0 0.6], ...
        "DisplayName", "SOM-RS mindeg 2", 'markersize', 11);
    lgd = legend;
    fontsize(lgd,6,'points')
    transl_fig_name = convertStringsToChars(strcat("transl_err_progression_sigma", string(sigma)));
    transl_fig_name(transl_fig_name=='.')=[];
    savefigure(transl_fig_name,'epsc',sz)
    hold off

end





end %file function
