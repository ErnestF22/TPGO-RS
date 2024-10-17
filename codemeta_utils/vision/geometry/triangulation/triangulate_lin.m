function X=triangulate_lin(x,P)
method='hom';
NX=size(x,2);
NP=size(P,3);
if NP<2
    error('Triangulation needs at least two views')
end
d=size(P,1)-1;
D=size(P,2)-1;

X=zeros(D,NX);
A=zeros(d*NP,D+1);
idxA=reshape(1:d*NP,d,NP);
for iX=1:NX
    xOrth=xOrthogonalHom(squeeze(x(:,iX,:)));

    for iP=1:NP
        A(idxA(:,iP),:)=xOrth(:,:,iP)'*P(:,:,iP);
    end
    switch method
        case 'hom'
            [U,S,V]=svd(-A,'econ');
            X(:,iX)=V(1:end-1,end)/V(end,end);
        case 'inhom'
            X(:,iX)=-A(:,1:end-1)\A(:,end);
    end
end
