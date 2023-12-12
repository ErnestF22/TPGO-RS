%function [pQueue]=priority_insert(pQueue,key,cost)
%Add an element to the queue.
function [pQueue]=priority_insert(pQueue,key,z)
cost = key(end);
pQueue=[pQueue;struct('V',key,'Z',z,'cost',cost)];
%'ACT',act