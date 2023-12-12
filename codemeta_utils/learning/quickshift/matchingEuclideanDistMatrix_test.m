function matchingEuclideanDistMatrix_test
X=randn(3,10);
membershipPrior=[ones(1,5) 2*ones(1,5)];

[D,DIndicator]=matchingEuclideanDistMatrix(X,'NBest',2,'membershipPrior',membershipPrior);
disp(euclideanDistMatrix(X))
disp(full(D))
disp(full(DIndicator)) 
disp(issparse(DIndicator))
disp(issparse(D))
DIndicator=full(DIndicator);
disp(all(sum(DIndicator)==4))
disp(all(DIndicator(sub2ind([10 10],1:10,1:10))==0))
disp(all(full(D(~DIndicator))'==0))
