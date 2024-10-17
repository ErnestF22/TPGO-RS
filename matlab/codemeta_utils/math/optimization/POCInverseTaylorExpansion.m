function POCInverseTaylorExpansion
load POCBlockJacobi_inputData

switch 'H'
    case 'X'
        X=randn(3);
        X=X+X';
        X=1*X/max(abs(eig(X)));
        eig(X)

        I=eye(size(X));
        XSum=matrixSeries(X,100);
        disp([pinv(I-X) XSum])
        
    case 'Xv'
        X=randn(3);
        X=X+X';
        X=X/max(abs(eig(X)));
        eig(X)

        I=eye(size(X));
        v=orthComplementProjector(null(I-X))*rand(3,1);

        XvSum=matrixVectorSeries(X,v,10000);
        disp([pinv(I-X)*v XvSum])

    case 'XH'
        w=2;
        d=size(H,1);
        I=eye(d);

        D=w*blkZeroOffDiag(H,2);
        B=D-H;
        Dhalf=sqrtPSDMatrix(D);

        X=(Dhalf\B)/Dhalf;
        X=(X+X')/2; %this just for numerical stability
        PX=orthComplementProjector(null(I-X));
        X=PX*X*PX';
        XSum=matrixSeries(X,10000);

        v=PX*rand(d,1);
        
        disp([pinv(I-X)*v XSum*v]')
        
        
    case 'XHv'

        w=2;
        d=size(H,1);
        I=eye(d);

        D=w*blkZeroOffDiag(H,1);
        B=D-H;
        Dhalf=sqrtPSDMatrix(D);

        X=(Dhalf\B)/Dhalf;
        X=(X+X')/2; %this just for numerical stability
        eig(X)
        v=orthComplementProjector(null(I-X))*rand(d,1);
        XvSum=matrixVectorSeries(X,v,10000);
        disp([pinv(I-X)*v XvSum]')
        
    case 'H'
        w=2;
        d=size(H,2);
        I=eye(d);

        D=w*blkZeroOffDiag(H,1);
        B=D-H;
        Dhalf=sqrtPSDMatrix(D);
        DhalfInv=inv(Dhalf);

        X=(Dhalf\B)/Dhalf;
        X=(X+X')/2; %this just for numerical stability
        PX=orthComplementProjector(null(I-X));
        X1=X;%PX*X*PX';
        XSum=matrixSeries(X1,5000);
        PH=orthComplementProjector(null(H));
        HInv=((Dhalf\XSum)/Dhalf);
        v=PH*rand(d,1);
        
        disp([pinv(H)*v PH*HInv*v]')

    case 'Hv'

        w=2;
        d=size(H,1);
        I=eye(d);

        D=w*blkZeroOffDiag(H,1);
        B=D-H;
        Dhalf=sqrtPSDMatrix(D);

        X=(Dhalf\B)/Dhalf;
        X=(X+X')/2; %this just for numerical stability
        PX=orthComplementProjector(null(I-X));
        X1=PX*X*PX';
        v=PX*rand(d,1);
        v=Dhalf\v;
        XvSum=matrixVectorSeries(X,v,10000);
        HInvv=Dhalf\XvSum;
        
        disp([pinv(H)*v HInvv]')
        
end

function Ainv=truncatedSPDInverse(A,n)
[U,S,V]=svd(A);
Ainv=V(:,1:end-n)*diag(1./diag(S(1:end-n,1:end-n)))*U(:,1:end-n)';


% 
% x=randn(d,1);
% g=H*x;
% 
% NIt=300;
% u=zeros(size(x,1),NIt);
% u(:,1)=Dhalf\g;
% for it=2:NIt
%     u(:,it)=(X+I)*u(:,it-1);
% end
% 
% [~,e]=cnormalize(H*Dhalf\u-repmat(g,1,NIt));
% subplot(2,1,1)
% semilogy(e)
% subplot(2,1,2)
% plot(g'*u)

function XSum=matrixSeries(X,N,XSum0)
if ~exist('XSum0','var')
    XSum=eye(size(X));
else
    XSum=XSum0;
end
XProd=eye(size(X));
for iN=1:N
    XProd=X*XProd;
    XSum=XSum+XProd;
end

function XvSum=matrixVectorSeries(X,v,N)
XvSum=v;
XvProd=v;
for iN=1:N
    XvProd=X*XvProd;
    XvSum=XvSum+XvProd;
end



