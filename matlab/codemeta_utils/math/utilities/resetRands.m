%function resetRands(seed)
%Resets the states of rand and randn to SEED (default, SEED=0)

%%AUTORIGHTS%%

function resetRands(seed)
if nargin<1
    seed=0;
end

randn('state',seed)
rand('state',seed)
