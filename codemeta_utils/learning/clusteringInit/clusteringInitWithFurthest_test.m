function clusteringInitWithFurthest_test
switch 2
    case 1
        resetRands()
        K=10;
        A=randn(2,100);

        X=clusteringInitWithFurthest(A,K);

        plotPoints(A);
        hold on
        plotPoints(X,'r')

        for iK=1:K
            text(X(1,iK)+0.05,X(2,iK),num2str(iK))
        end
        hold off
    case 2
        I=eye(4);
        A=I(:,[1 1 1 1 1 2 3 3 3 3 4 4 4 4]);
        XInit=clusteringInitWithFurthest(A,4);
        [membership,X]=kmeans(A',4,'start',XInit');
        X=X';
        membership=membership';
        disp(membership)
        disp(X)
        
end