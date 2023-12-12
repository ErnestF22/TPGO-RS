clear all 
close all
% arbitrary input to the network
x = randi([2,10],[2,1]);

L = 3;
d = 2;
% The networkis defined as:
A1 = randi([-4,15],[2,2]);
b1 = randi([0,10],[2,1]);

A2 = randi([-5,13],[2,2]);
b2 = randi([2,10],[2,1]);

A3 = randi([5,18],[1,2]);
b3 = randi([0,10],1);

Ab_set = struct;
for i=1:L
    Ab_set(i).A = eval(['A' num2str(i)]);
    Ab_set(i).b = eval(['b' num2str(i)]);
end
ReluPlexPlotNetwork(Ab_set,L);
% 1) x is the input of the network, given x find which region x belongs to:
[y,Ab_set] = ReluPlexNetworkOutput(x,Ab_set,L);

z=[];
for i=1:L
    z = [z;Ab_set(i).z];
end

% 2) find all the vertices of the region that is defined by z from step2:
[V] = ReluPlexFindAllVerticesOfRegion(z,Ab_set,L);
