function gradc=rotLocRiemannianCostPair_grad(Ri,Rj,Rij,funs,e)
if ~exist('e','var') || isempty(e)
    e=rot_dist(Rj,Ri*Rij);
end

if e==0
    gradc=zeros([size(Ri) 2]);
else
    gei=rot_log(Ri,Rj*Rij');
    gej=rot_log(Rj,Ri*Rij);
    dfc=funs.df(e);
    gradc=-dfc*cat(3,gei/e,gej/e);
end
