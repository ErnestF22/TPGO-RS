function quickshift_density_testStress
global phi
resetRands()
N=20000;
d=128;

phi=@(x) exp(x.^2/2);
f=@(x) log(1+x);

NTrials=1;
for iTrial=1:NTrials
    disp(['# Trial ' num2str(iTrial)])
    data=randn(d,N);
    scales=1+rand(1,N);
    opts={'scales',scales,'amplify',f};
    disp('Low memory')
    tic
    lowMemoryTest(data,opts)
    toc
    disp('Full distance matrix')
    tic
    %distanceTest(data,opts)
    toc
end

function distanceTest(data,opts)
global phi
D=sqrt(euclideanDistMatrix(data));
density=quickshift_density(phi,D,opts{:});

function lowMemoryTest(data,opts)
global phi
density=quickshift_density(phi,data,'lowMemory',opts{:});


