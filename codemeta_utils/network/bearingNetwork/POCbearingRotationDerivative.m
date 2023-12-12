function POCbearingRotationDerivative
[Ri,dRi,~,~,wi]=rot_randGeodFun(eye(3));

pi=randn(3,1);
pj=randn(3,1);
em=[0;0;1];
embij=@(t) em'*Ri(t)'*cnormalize(pj-pi);

skew=@(A) (A-A')/2;
%dembij=@(t) trace(hat3(wi)'*);
%funCheckDer(embij,dembij)
gradvRi_embij=@(t) 2*vee3(skew(Ri(t)'*cnormalize(pj-pi)*em'));
funCheckDer(embij,@(t) wi'*gradvRi_embij(t));


