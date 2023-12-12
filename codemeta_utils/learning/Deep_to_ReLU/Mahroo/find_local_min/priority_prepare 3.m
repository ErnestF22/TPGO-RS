%function [pQueue]=priority_prepare()
%Create an empty structure array for the queue.
function [pQueue]=priority_prepare()
%'ACT',[]
pQueue=repmat(struct('V',[],'Z',[],'cost',[]),0,1);