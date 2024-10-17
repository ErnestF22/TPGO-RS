%Remove duplicate points
%function p=uniqueApprox(p,tol)
%Produce a version of p where points with other points at less than tol
%distance are removed and replaced with their average.
function p=uniqueApprox(p,tol)
if ~exist('tol','var') || isempty(tol)
    tol=1e-14;
end

d=sqrt(euclideanDistMatrix(p));

idxCurrent=1;
while idxCurrent<=size(p,2)
    idxRemove=idxCurrent+find(d(idxCurrent,idxCurrent+1:end)<tol);
    d(idxRemove,:)=[];
    d(:,idxRemove)=[];
    p(:,idxCurrent)=mean(p(:,[idxCurrent idxRemove]),2);
    p(:,idxRemove)=[];
    idxCurrent=idxCurrent+1;
end
