function quickshift_density_test_lowMemory
resetRands()
N=100;
d=2;

data=randn(d,N);
scales=1+rand(1,N);
phi=@(x) exp(x.^2/2);
f=@(x) log(1+x);

D=sqrt(euclideanDistMatrix(data));
density=quickshift_density(phi,D,'scales',scales,'amplify',f);

densityLowMemory=quickshift_density(phi,data,'lowMemory','numelDataChunk',10,...
    'scales',scales,'amplify',f);

disp(norm(density(:)-densityLowMemory(:)))
