function [zijk,cz]=admmMedoids_auxiliary(muijk,lambdaijk,rho)
d=size(muijk,1);
zijk=zeros(d,1);
cz=0;
for iDim=1:d
    c=admmMedoids_auxiliaryCost(muijk,lambdaijk,rho,'dim',iDim);
    [cMin,idxMin]=min(c);
    cz=cz+cMin;
    zijk(iDim)=muijk(iDim,idxMin);
end