function [c,gradc]=rotLocRiemannianCostPair(Ri,Rj,Rij,funs)
flagComputeGrad=nargout>1;

e=rot_dist(Rj,Ri*Rij);
if e==0
    c=0;
    gradc=zeros([size(Ri) 2]);
else
    c=funs.f(e);

    if flagComputeGrad
        gradc=rotLocRiemannianCostPair_grad(Ri,Rj,Rij,funs,e);
    end
end
