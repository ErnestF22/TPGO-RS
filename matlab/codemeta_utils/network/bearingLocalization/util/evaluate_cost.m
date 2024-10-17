function cost = evaluate_cost(Adj,F,M,R)

symAdj = (.5*(Adj + Adj') + 10^(-5)) > 1;   

nE = sum(symAdj(:))/2;
  
nF = zeros(numel(F),1);

for iF=1:numel(F)
    nF(iF) = numel(F{iF});
end

sum_nF = sum(nF);

cost =.5*trace(myreshape3d2dcolumnwise(multitransp(R))'*M*myreshape3d2dcolumnwise(multitransp(R)))+ nE + sum_nF;


end

