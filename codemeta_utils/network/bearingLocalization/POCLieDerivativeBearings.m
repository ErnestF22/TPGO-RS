%Test the Lie derivatives of the bearing measurements
function POCLieDerivativeBearings
[Ri,dRi,Ri0,dRi0,wi]=rot_randGeodFun();
[Ti,dTi,Ti0,vi]=real_randGeodFun(randn(3,1));
[Tj,dTj,Tj0,vj]=real_randGeodFun(randn(3,1));

E=eye(3);

tij=@(t) cnormalize(Ri(t)'*(Tj(t)-Ti(t)));
dij=@(t) norm(Ri(t)'*(Tj(t)-Ti(t)));
hij=@(k,t) E(:,k)'*tij(t);

%dijSq=@(t) dij(t)^2;
%ddijSq=@(t) 2*(Tj(t)-Ti(t))'*(vj-vi);

%ddij=@(t) 1/dij(t)*(Tj(t)-Ti(t))'*(vj-vi);
ddij=@(t) tij(t)'*Ri(t)'*(vj-vi);

funCheckDer(dij,ddij)



