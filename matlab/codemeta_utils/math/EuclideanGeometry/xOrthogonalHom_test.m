function xOrthogonalHom_test
NX=5;
x=randn(2,NX);

xOrth=xOrthogonalHom(x);
for iX=1:NX
    disp([x(:,iX);1]'*xOrth(:,:,iX))
end

