function Rij=bnFranchi_RzFromBearings(bij,bji)
d=size(bij,1);
bij=[cnormalize(bij(1:2,:));0];
bji=[cnormalize(bji(1:2,:));0];
cij=bearingComputeCosine(-bij,bji);
vij=cnormalize(cross(bij,bji));
Rij=rot3_exp(acos(cij)*hat3(vij));
if d==2
    Rij=Rij(1:2,1:2);
end

