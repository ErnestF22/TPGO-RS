function [T,H]=bearingNetworkPreconditioner(t_node,funs,type,varargin)
NNodes=t_node.NNodes;
E=t_node.E;
yg=t_node.Yijtruth;
nyg=t_node.nYijtruth;
d=size(yg,1);

[~,DgradPhi]=bearingNetworkCostGradient(E,yg,yg,funs,nyg);
H=(DgradPhi+DgradPhi)/2;
I=eye(size(H));

idxNodes=reshape(1:d*NNodes,d,NNodes);
N=null(H);
U=orthComplement(N);
IReduced=eye(size(U,2));

cvx_quiet(true)
switch lower(type)
    case 'eye'
        T=eye(size(H));
    case 'hop0'
        T=zeros(size(H));
        for iNode=1:NNodes
            Hi=H(idxNodes(:,iNode),:);
            Ii=I(idxNodes(:,iNode),:);

            cvx_begin
                variable Ti(d,d) symmetric
                minimize norm(Ti*Hi-Ii,'fro')
                subject to
                    Ti==semidefinite(d);
            cvx_end

            T(idxNodes(:,iNode),idxNodes(:,iNode))=Ti;
        end
    case 'hop1'
        cvx_begin
            variable T(d*NNodes,d*NNodes) symmetric
            minimize norm(T*H-I,'fro')
            subject to
                T(H==0)==0;
                T==semidefinite(d*NNodes);
        cvx_end
    case 'hop1n'
        cvx_begin
            variable T(d*NNodes,d*NNodes) symmetric
            minimize norm(U'*T*H*U-IReduced,'fro')
            subject to
                T(H==0)==0;
                %T*N==zeros(size(N));
                T==semidefinite(d*NNodes);
        cvx_end
    case 'hop1p1'
        cvx_begin
            variable T1(d*NNodes,d*NNodes) symmetric
            minimize norm(T1*H-I,'fro')
            subject to
                T1(H==0)==0;
                T1==semidefinite(d*NNodes);
        cvx_end
        cvx_begin
            variable T2(d*NNodes,d*NNodes) symmetric
            minimize norm(T2*T1*H-I,'fro')
            subject to
                T2(H==0)==0;
                T2==semidefinite(d*NNodes);
        cvx_end
        T=T2*T1;
    case 'hop1np1'
        cvx_begin
            variable T1(d*NNodes,d*NNodes) symmetric
            minimize norm(U'*T1*H*U-IReduced,'fro')
            subject to
                T1(H==0)==0;
                T1==semidefinite(d*NNodes);
        cvx_end
        cvx_begin
            variable T2(d*NNodes,d*NNodes) symmetric
            minimize norm(U'*T2*T1*H*U-IReduced,'fro')
            subject to
                T2(H==0)==0;
                T2==semidefinite(d*NNodes);
        cvx_end
        T=T2*T1;
    case 'ideal'
        cvx_begin
            variable T(d*NNodes,d*NNodes) symmetric
            minimize norm(T*H-I,'fro')
            subject to
                T==semidefinite(d*NNodes);
        cvx_end
    case 'idealn'
        cvx_begin
            variable T(d*NNodes,d*NNodes) symmetric
            minimize norm(U'*T*H*U-IReduced,'fro')
            subject to
                T==semidefinite(d*NNodes);
        cvx_end
    case 'idealpinv'
        T=pinv(H);
    case 'neumann'
        order=varargin{1};
        w=2;

%         D=w*blkZeroOffDiag(H,d);
%         B=D-H;
%         Dhalf=sqrtPSDMatrix(D);
% 
%         X=(Dhalf\B)/Dhalf;
        [X,Dhalf]=makeX(H,d,w);
        X=(X+X')/2; %this just for numerical stability
        XSum=matrixSeries(X,order);
        T=((Dhalf\XSum)/Dhalf);
    case 'neumannopt'
        %get order as the next parameter
        order=varargin{1};
        L=makeL(H,d);
        sL=eig(L);
        sLThresholded=sL(sL>1e-4);
        w=(max(sLThresholded)+min(sLThresholded))/2;

        [X,Dhalf]=makeX(H,d,w);
        X=(X+X')/2; %this just for numerical stability
        XSum=matrixSeries(X,order);
        T=((Dhalf\XSum)/Dhalf);
    otherwise
        error('Preconditioner type not found')
end

function L=makeL(H,d)
D=blkZeroOffDiag(H,d);
Dhalf=sqrtPSDMatrix(D);
L=(Dhalf\H)/Dhalf;

function [X,Dhalf]=makeX(H,d,w)
D=w*blkZeroOffDiag(H,d);
B=D-H;
Dhalf=sqrtPSDMatrix(D);

X=(Dhalf\B)/Dhalf;


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
        
