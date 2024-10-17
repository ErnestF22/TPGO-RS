function stackedVariables_positionsGetPair_test
variables=stackedVariables_set([],{{'x',3,2},{'y',3},{'c1'}});
d=stackedVariables_stackSize(variables);
M=reshape(1:d^2,d,d);
allDescription={{'x'},{'y'},{'c1'}};
[idxStacked,idxStackedTranspose]=stackedVariables_positionsGetPair(variables,allDescription,{'x',2});

disp(M)
disp(M(idxStacked))
disp(M(idxStackedTranspose))
