function acc15_example_bowtie

[A,x]=bearingCluster_generateTest('bowtie');
E=adj2edges(A,'oriented');
u=bearingCluster_getBearingsScalesFromE(x,E);
M=bearingCluster_measurementMatrixFromE(E,u);
L=null(M);

Lij=L(7,:);
N=null(Lij);

l1=cnormalize(L*N*[1;-1]);
l2=cnormalize(l1-0.5*L*Lij');
xRec=bearingCluster_scales2nodes(E,u,l2);

gshow(E,'coords',xRec')

