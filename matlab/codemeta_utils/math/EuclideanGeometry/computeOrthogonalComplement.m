function xOrth=computeOrthogonalComplement(x)
[dimData,NX]=size(x);
d=cnormalize(x);
d(end,:)=d(end,:)-1;
d=cnormalize(d);
I=eye(dimData,dimData-1);
xOrth=zeros(dimData,dimData-1,NX);
for iX=1:NX
    xOrth(:,:,iX)=I-2*d(:,iX)*d(1:end-1,iX)';
end

