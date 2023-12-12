%function cycle1=reverse_cycle(cycle)
%Reverse the order of a single cycle
function cycle1=reverse_cycle(cycle)
cycle1=abs(cycle);
cycle1=max(cycle1)-cycle1+1;
cycle1=cycle1.*-sign(cycle);
