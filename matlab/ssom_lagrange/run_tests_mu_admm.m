clear; 
clc;

close all;

mus_admm = [0.1, 0.2, 0.5, 0.8, 1, 2, 3, 4, 5, 10];

for mu_admm = mus_admm

    mu = mu_admm;

    run_lsom_sigma_tests

end

