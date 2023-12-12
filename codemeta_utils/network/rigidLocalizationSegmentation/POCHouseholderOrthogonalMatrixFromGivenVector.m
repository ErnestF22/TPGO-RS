function POCHouseholderOrthogonalMatrixFromGivenVector
resetRands();
v=randn(3,1);
vn=cnormalize(v);

N=size(v,1);

e=zeros(N,1);
e(end)=1;

d=cnormalize(vn-e);

H=eye(N)-2*(d*d');
Horth=eye(N,N-1)-2*d*d(1:end-1)';
disp([vn H*e])
disp(Horth'*v)
