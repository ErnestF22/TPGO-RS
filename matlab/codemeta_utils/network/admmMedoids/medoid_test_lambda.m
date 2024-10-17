function medoid_test_lambda
%Test that with lambda=0 we get the same result

%resetRands()
nbPoints=100;
nbClusters=4;
dimPoints=2;
allOptsMedoids={{}};

x=rand(dimPoints,nbPoints);
muInit=zeros(dimPoints,nbClusters);

lambda=repmat({[0;0]},1,nbClusters);
allOptsMedoids{end+1}={'bias',lambda};

for iTest=1:length(allOptsMedoids)
    [~,output]=medoids(x,nbClusters,'muInit',muInit,'debug',allOptsMedoids{iTest}{:});
    disp(output.c)
end
