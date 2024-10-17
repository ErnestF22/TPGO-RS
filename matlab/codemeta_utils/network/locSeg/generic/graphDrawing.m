function x=graphDrawing(A,dim)
if size(A,1)<dim+1
    x=randn(dim,size(A,1));
    disp('Note: graph drawing requires at least one additional node')
    return
end
L=adj2laplacianmatrix(A);
[V,E]=eig(L);
x=V(:,end-dim:end-1)';
