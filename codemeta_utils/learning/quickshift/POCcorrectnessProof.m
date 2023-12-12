function POCcorrectnessProof
distIntra=4;
distInter=distIntra/16;

ratioDensity=0.25;
ratioInterCluster=0.67;

X=[0 distIntra];
X=[X X(1)+distInter     X(2)-3*distInter];
X=[X X(1)+2*distInter   X(2)-2*distInter];
X=[X X(1)+3*distInter   X(2)-distInter];
X=[X;zeros(1,size(X,2))];
X=X+0.01*randn(size(X));

membershipPrior=[1 1 2 2 3 3 4 4];
%plotGroups(X,membershipPrior)

D=sqrt(euclideanDistMatrix(X,X));
phi=@(x) max(0,1-(x/2).^2);
scales=quickshift_scalesMembershipPrior(D,membershipPrior);
%scales=distIntra/2;

figure(1)
quickshift_plotDensity(X,phi,'optsDensity',{'scales',ratioDensity*scales},...
    'limits', [-1,-distIntra/2,distIntra+1,distIntra/2]);
axis equal

paramsMatching={'phi',phi,...
    'ratioDensity',ratioDensity,...
    'ratioInterCluster',ratioInterCluster};

[membershipMatches,info]=quickshift_matching(D,membershipPrior,...
    paramsMatching{:});

figure(2)
quickshift_plotMatch(X,membershipMatches,info)