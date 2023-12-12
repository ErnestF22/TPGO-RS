% This function finds a minimum of a network by using DFS
% function [V_total,Z_total] = main(x,L,Ab_set,flagPlot)
clear all 
close all
% arbitrary input to the network
% x = [7.5;3.8];
x = [5;-10];
% x = [-4;-0.9];
flagPlot = 1;
flagAdded = 0;
L = 3+flagAdded;

% The networkis defined as:
Ab_set = ReluPlexNetwork(L);

% 0) Plot the network:
if flagPlot
    ReluPlexPlotNetwork(Ab_set,L-flagAdded);
end

% 1) x is the input of the network, given x find which region x belongs to:
[y,Ab_set] = ReluPlexNetworkOutput(x,Ab_set,L);
[yout,~] = ReluPlexNetworkOutput(x,Ab_set,L-flagAdded);

z=[];
for i=1:L
    z = [z;Ab_set(i).z];
end

% 2) find all the vertices of the region that is defined by z from step2:
[V] = ReluPlexFindAllVerticesOfRegion(z,Ab_set,L);


% 3) find all vertices of the network by DFS method:
[V_total,Z_total] = ReluPlexFindAllVerticesOfNetwork(V,Ab_set,L,flagAdded);

% 4) minimum vertex belongs to this region:
[~,idx] = min(V_total(:,3));
z_min = Z_total(idx,:);

% 5) plot the 
if flagPlot
    plot_vertices(V_total,Ab_set,L-flagAdded,'b')
    hold on
    plot_vertices(V_total,Ab_set,L-flagAdded,'k')
    hold on
    plot3(x(1),x(2),yout,'r*')
    hold off
end


% % 2*) find min robust distance:
% dRobust = ReluPlexMinRobustDis(As,[x;y]);
% x_test1 = x-[dRobust;dRobust];
% x_test2 = x+[dRobust;dRobust];
% x_test3 = x+[0;-dRobust];
% [y1,Ab_set1] = ReluPlexNetworkOutput(x_test1,Ab_set,L);
% [y2,Ab_set2] = ReluPlexNetworkOutput(x_test2,Ab_set,L);
% [y3,Ab_set3] = ReluPlexNetworkOutput(x_test3,Ab_set,L);

% hold on
% plot3(x(1),x(2),y,'r*')
% hold on
% plot3(x_test1(1),x_test1(2),y1,'b*')
% hold on
% plot3(x_test2(1),x_test2(2),y2,'g*')
% hold on
% plot3(x_test3(1),x_test3(2),y3,'k*')


