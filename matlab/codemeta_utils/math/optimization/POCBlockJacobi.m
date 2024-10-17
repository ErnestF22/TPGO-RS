function POCBlockJacobi
resetRands();
load POCBlockJacobi_inputData

T=H;
I=eye(size(T));
D=blkZeroOffDiag(H,1);
R=D-H;
DinvH=D\H;
x=randn(size(H,1),1);
g=H*x;
w=1/max(eig(DinvH));
f=w*(D\g);
J=(1-w)*I+w*(D\R);

NIt=300;
u=zeros(size(x,1),NIt);
u(:,1)=g;
for it=2:NIt
    u(:,it)=J*u(:,it-1)+f;
end

[~,e]=cnormalize(H*u-repmat(g,1,NIt));
subplot(2,1,1)
semilogy(e)
subplot(2,1,2)
plot(g'*u)

keyboard
