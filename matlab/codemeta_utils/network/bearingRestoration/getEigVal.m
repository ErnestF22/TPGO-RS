function [eigval]=getEigVal(Erigid,u)
CUnsigned=grCycleBasisMultiComp(Erigid);
C=grOrientCycleBasis(CUnsigned,Erigid)';
M=bearingCluster_measurementMatrixFromC(C,u);
s2=svd(M);
dim1=max(size(M));
dim2=min(size(M));
% if dim1>dim2
%     eigval=[s2;zeros(dim1-dim2,1)];
% else
%     eigval=s2;
% end
%eigval=sort(eigval);
eigval=sort(s2);
end