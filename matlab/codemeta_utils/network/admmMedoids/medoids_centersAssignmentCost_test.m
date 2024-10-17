function medoids_centersAssignmentCost_test
resetRands()
nbPoints=20;
dimPoints=2;
kCenter=1;
switch dimPoints
    case 1
        error('Not implemented')
    case 2
        x=rand(2,nbPoints)*10;
        mu=repmat([5 2 3],[2 1]);
        xEval=linspace(0,10,100);
        muEval=@(x) [mu(:,1:(kCenter-1)) x mu(:,(kCenter+1):end)];
        f=@(xEval) medoids_centersAssignmentCost(x,muEval(xEval));
        funImage(xEval,xEval,f,'method','surf')
end        
