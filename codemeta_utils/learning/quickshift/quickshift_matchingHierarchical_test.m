function quickshift_matchingHierarchical_test
resetRands()
NPointsCluster=20;
[X,membershipPrior]=quickshift_test_datasets('matching',...
    'NPointsClass',NPointsCluster,'nClasses',5);

optsMatching={'gaussian',...
    'ratioDensity',0.25,...
    'ratioInterCluster',0.67,...
    };


[membershipCorrespondences]=quickshift_matchingHierarchical(X,membershipPrior,'optsMatching',optsMatching);

plotGroups(X,membershipCorrespondences)
