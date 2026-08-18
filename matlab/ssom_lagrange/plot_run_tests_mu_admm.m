function plot_run_tests_mu_admm

clc;
clear;
close all;

% sigmas = [0, ...
% 0.0100000000000000, ...
% 0.0200000000000000, ...
% 0.0500000000000000, ...
% 0.100000000000000, ...
% 0.200000000000000, ...
% 0.500000000000000, ...
% 1, ...
% 2];

% num_tests_per_sigma = 50;

% mus = zeros(size(sigmas));

fileList = dir('*.mat');

% mus_admm = [0.100000000000000	0.200000000000000	0.500000000000000	0.800000000000000	1	2	3	4	5	10];


for jj = 1:size(fileList, 1)

    close all;

    f_name = fileList(jj).name;
    
    load(f_name, "mus_admm")
    load(f_name, "sigmas")
    load(f_name, "mus")
    load(f_name, "num_tests_per_sigma")
    load(f_name, "ssom_rot_errs")
    load(f_name, "ssom_transl_errs")
    load(f_name, "ssom_exec_times")
    load(f_name, "ssom_scale_ratios")
    load(f_name, "ssom_transl_errs_norm")
    load(f_name, "results")

    
    disp("jj")
    disp(jj)
    
    disp("filename")
    disp(f_name)
    mu_admm = mus_admm(jj);
    disp("mu_admm")
    disp(mu_admm)    

    disp("sigmas");
    disp(sigmas);
    disp("mus");
    disp(mus);
    disp("num_tests_per_sigma");
    disp(num_tests_per_sigma);
    
    % plot 
    % on x sigmas
    % on y the mean errors and the mean execution times (across all iterations) 
    
    % disp("manopt_sep_rot_errs");
    % disp(manopt_sep_rot_errs);
    % 
    % disp("manopt_sep_transl_errs");
    % disp(manopt_sep_transl_errs);
    % 
    % disp("procrustes_rot_errs");
    % disp(procrustes_rot_errs);
    % 
    % disp("procrustes_transl_errs");
    % disp(procrustes_transl_errs);
    
    disp("ssom_rot_errs");
    disp(ssom_rot_errs);
    
    disp("ssom_transl_errs");
    disp(ssom_transl_errs);
    
    % disp("manopt_sep_exec_times");
    % disp(manopt_sep_exec_times);
    % 
    % disp("procrustes_exec_times");
    % disp(procrustes_exec_times);
    
    disp("ssom_exec_times");
    disp(ssom_exec_times);
    
    disp("ssom_scale_ratios");
    disp(ssom_scale_ratios);
    
    disp("ssom_transl_errs_norm");
    disp(ssom_transl_errs_norm);
    
    
    % results = struct("manopt_sep_rot_errs", manopt_sep_rot_errs, ...
    %     "manopt_sep_transl_errs", manopt_sep_transl_errs, ...
    %     "manopt_sep_exec_times", manopt_sep_exec_times, ...
    %     "procrustes_rot_errs", procrustes_rot_errs, ...
    %     "procrustes_transl_errs", procrustes_transl_errs, ...
    %     "procrustes_exec_times", procrustes_exec_times, ...
    %     "ssom_rot_errs", ssom_rot_errs, ...
    %     "ssom_transl_errs", ssom_transl_errs, ...
    %     "ssom_exec_times", ssom_exec_times, ...
    %     "ssom_scale_ratios", ssom_scale_ratios, ...
    %     "ssom_transl_errs_norm", ssom_transl_errs_norm);
    
    %plot results
    plot_results_ssom(sigmas, results, "ssom");

    disp("jj")
    disp(jj)
    f_name = fileList(jj).name;
    disp("filename")
    disp(f_name)
    mu_admm = mus_admm(jj);
    disp("mu_admm")
    disp(mu_admm)
end

end % file function