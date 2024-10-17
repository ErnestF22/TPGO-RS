function computeAllEuclideanDistancesSq_test
x=[-1  0  0  1  1; 
    1  1  0  0 -1];

allDists=computeAllEuclideanDistancesSq(x);

display(x)
disp('Squared distances')
disp(allDists)
