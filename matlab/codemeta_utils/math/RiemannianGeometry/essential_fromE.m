%Compute the point in the QREM corresponding to a given essential matrix E
%function Q=essential_fromE(E)
function Q=essential_fromE(E)
NE=size(E,3);
Q=zeros(6,3,NE);
Rz=blkdiag([0 1;-1 0],1);
for iE=1:NE
    [U,~,V]=svd(E(:,:,iE));
    if det(U)<0
        U=U*diag([1;1;-1]);
    end
    if det(V)<0
        V=V*diag([1;1;-1]);
    end
    Q(1:3,:,iE)=U';
    Q(4:6,:,iE)=Rz*V';
end