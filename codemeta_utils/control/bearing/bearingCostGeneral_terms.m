function [gradi,Hi]=bearingCostGeneral_terms(yi,ygi,fci,dfci,ci,ddfci)
flagComputeDGrad=false;
if nargout>1
    flagComputeDGrad=true;
end

Pyi=(eye(size(yi,1))-yi*yi');
Pyiygi=Pyi*ygi;
if fci==0
    gradi=zeros(size(yi));
else
    %gradi=-fci*yi-dfci*Pyiygi;
    gradi=(-fci+dfci*yi'*ygi)*yi-dfci*ygi;
end
if flagComputeDGrad
    Hi=ddfci*(Pyiygi*Pyiygi')-(dfci*ci-fci)*Pyi;
end
