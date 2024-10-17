function [flag,evals]=isDirectedNetworkStable(E)
%flag to remove known zero root
flagRemoveZeroRoot=true;
%threshold for deciding if real parts of evals are negative
evalsThreshold=0;

G=directedNetworkStabilityMatrix(E);
evals=eig(G);

if flagRemoveZeroRoot
    [~,idxZeroRoot]=min(abs(evals));
    evals(idxZeroRoot(1))=[];
end

flag=all(real(evals)<evalsThreshold);

