clear all 
close all
% arbitrary input to the network
x = [7.5;3.8];
% x = [5;-10];
L = 3;
% The networkis defined as:
Ab_set = ReluPlexNetwork(L);
flagAddLayer = 0;
flagPlot = 1;
[V_total,Z_total]= main(x,L,Ab_set,flagPlot);