%Timing QuickMatch on a large synthetic dataset
function quickshift_matching_testTiming
resetRands(1)

%number of points for image
NPointsImage=4000;

%parameters for 
ratioDensity=0.25;
ratioInterCluster=0.67;
thresholdBreakTree=Inf;

%generate points in 4 images
[X,membershipPrior]=quickshift_test_datasets('matching',...
    'NPointsClass',NPointsImage,'nOutliersClass',10);

tic
%compute pairwise distances
D=sqrt(euclideanDistMatrix(X,X));
fprintf('Time after computing distances: ')
toc

%Do the clustering
[membershipMatches,info]=quickshift_matching(D,membershipPrior,...
    'gaussian','ratioDensity',ratioDensity,...
    'ratioInterCluster',ratioInterCluster,...
    'threshold',thresholdBreakTree,...
    'densitySqrtAmplify');
fprintf('Time for entire matching: ')
toc
