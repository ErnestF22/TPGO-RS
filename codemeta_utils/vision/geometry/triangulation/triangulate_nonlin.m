function [XEst,resnorm]=triangulate_nonlin(x,P,X)
d=size(x,1);
NX=size(x,2);

XEst=zeros(d+1,NX);
resnorm=zeros(1,NX);
opts=optimset('Display','Off','Jacobian','On');

if ~exist('X','var')
    X=triangulate_lin(x,P);
end

for iX=1:NX
    [XEst(:,iX),resnorm(iX)]=lsqnonlin(@(z) errorVect(x(:,iX,:),P,z), X(:,iX), [], [], opts);
end

function [errVect,JErrVect]=errorVect(xTruth,P,X)
[err,JxProjected]=reprojectionError(xTruth,P,X);
errVect=err(:);
JErrVect=reshape(permute(squeeze(JxProjected),[1 3 2]),...
    size(JxProjected,1)*size(JxProjected,4),[]);
