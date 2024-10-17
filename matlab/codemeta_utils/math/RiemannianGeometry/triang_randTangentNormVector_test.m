function triang_randTangentNormVector_test
x=rand(2,1)*rand(1,2);
f=@() triang_tangentProjVerticalScalar(x,triang_randTangentNormVector(x));
plotfuntrials(f)
