function acc15_cycleExample

x=[1 0 0; 0 0 0; 0 1 0; 1 1 0; 1 1 1]';
[d,N]=size(x);

A=adjgallery(N,'kneigh',1);
[B,E]=adj2incmatrix(A);
[u,lambda]=bearingCluster_getBearingsScales(x,A);
U=bearingCluster_augmentedBearingMatrix(u);
C=ones(1,size(u,2));
Cd=kron(C,eye(d));

M=Cd*U;
L=bearingCluster_nullSpaceBasis(M);


Nt=11;
pt=linspace(-1,1,Nt);
dLambda=
figure(1)
for it=1:Nt
    

