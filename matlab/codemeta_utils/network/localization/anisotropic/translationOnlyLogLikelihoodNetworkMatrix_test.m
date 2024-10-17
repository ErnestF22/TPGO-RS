function translationOnlyLogLikelihoodNetworkMatrix_test
%test by using the matrix representation to compute the derivative
t_node=buildTestTagNetwork();
[Ri0,Ti0]=testNetworkGetRotTransl(t_node);
[Rij,Tij]=testNetworkGetRelativeRotTranslScales(t_node);
[GammaijR,GammaijT]=testNetworkGetDispersionMatricesRotationTranslation(t_node);
E=testNetworkGetEdges(t_node);

[Tit,dTit,Ti0,dTi]=real_randGeodFun(randn(size(Ti0)));

f=@(t) translationLogLikelihoodNetwork(Ri0,Tit(t),Tij,GammaijT,E);

[A,b]=translationOnlyLogLikelihoodNetworkMatrix(Ri0,Tij,GammaijT,E);
df=@(t) dTi(:)'*(A*reshape(Tit(t),[],1)+b);

check_der(f,df)
