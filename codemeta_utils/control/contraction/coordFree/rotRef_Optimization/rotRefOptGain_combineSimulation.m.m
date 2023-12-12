% rotRefOptGain_combineSimulation.m
% Combine the multiple simuation data into one data file and correct time

%% Clean up workspace
close all; clear all; clc;

SimTo_3pi_over_4 = load('GainSchedule_3pi_over_4.mat');
SimTo_pi_over_4 = load('GainSchedule_pi_over_4.mat');
SimTo_I = load('GainSchedule_I.mat');

x = [SimTo_3pi_over_4.x, SimTo_pi_over_4.x, SimTo_I.x];
t = [SimTo_3pi_over_4.t; SimTo_pi_over_4.t+SimTo_3pi_over_4.t(end);...
    SimTo_I.t+SimTo_pi_over_4.t(end)+SimTo_3pi_over_4.t(end)];

control = SimTo_3pi_over_4.control;

save('GainSchedule_Combined.mat');