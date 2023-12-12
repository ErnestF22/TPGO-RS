%function gshow_cycle(E,C,cycle,style)
%Show a single cycle on the graph having edges E
function gshow_cycle(E,C,cycle,style)
E=E(cycle~=0,:);
cycle=cycle(cycle~=0);
M=length(cycle);
E=[E(sub2ind([M 2],[1:M]',2-(cycle>0))) E(sub2ind([M 2],[1:M]',2-(cycle<0)))];
gshow(E,'coords',C,'style',style);

