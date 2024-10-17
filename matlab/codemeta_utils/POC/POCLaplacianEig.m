function POCLaplacianEig
A=adjgallery(7,'locrand',0.7);
A=(A+A')>0;
D=randn(size(A));
D(~A)=0;
idxNewEdge=find(~(A+eye(size(A))));
idxNewEdge=idxNewEdge(randperm(length(idxNewEdge),2));
D(idxNewEdge)=rand(1,length(idxNewEdge));
D=(D+D')/2;
if min(D(:))<0
    D=D/max(abs(min(D(:))));
end

disp([A D])
f=@(t) secondEigLaplacian(A+t*D);

plotfun(f,linspace(0,1))



function lambda2=secondEigLaplacian(A)
L=adj2laplacianmatrix(A);
lambdas=eig(L);
lambda2=lambdas(2);
