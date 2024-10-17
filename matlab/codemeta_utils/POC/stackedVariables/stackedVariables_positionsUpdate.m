function variables=stackedVariables_positionsUpdate(variables)
cnt=0;
for iVariable=1:length(variables)
    d=variables(iVariable).size;
    n=variables(iVariable).cardinality;
    variables(iVariable).position=cnt+reshape(1:(d*n),d,n);
    cnt=cnt+d*n;
end