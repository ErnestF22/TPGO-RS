%Solve flip ambiguity in the QREM by using chirality constraint
%function [QSigned,l,lMask]=essential_solveFlipAmbiguity(Q,x)
function [QSigned,l,lMask]=essential_solveFlipAmbiguity(Q,x)
Nx=size(x,2);
Ns=zeros(1,4);
QSigned=zeros(6,3,4);
l=zeros(2,Nx,4);
lMask=zeros(4,Nx);
for k=1:4
    QSigned(:,:,k)=essential_flipAmbiguity(Q,k);
    l(:,:,k)=essential_triangulateDepths(QSigned(:,:,k),x);
    lMask(k,:)=(l(1,:,k)>0) & (l(2,:,k)>0);
    Ns(:,k)=sum(lMask(k,:));
end

maxNs=max(Ns);
idxMax=find(Ns==maxNs);
switch length(idxMax)
    case 0
        QSigned=[];
    case 1
        QSigned=QSigned(:,:,idxMax);
        l=l(:,:,idxMax);
        lMask=lMask(idxMax,:);
    otherwise
        warning('Multiple solutions with the same number of positive dephts (Ns=%d)',Ns(idxMax))
        QSigned=QSigned(:,:,idxMax);
        l=l(:,:,idxMax);
        lMask=lMask(idxMax,:);
end
