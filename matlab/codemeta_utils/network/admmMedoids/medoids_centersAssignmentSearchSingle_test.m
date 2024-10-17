function medoids_centersAssignmentSearchSingle_test
resetRands()
nbPoints=20;
dimPoints=2;
kCenter=1;
flagSurfCost=true;

testOpts=3;

switch testOpts
    case 1
        opts={};
    case 2
        opts={'bias',[20;0]};
    case 3
        opts={'priorCenters',[7;2],1};
end

switch dimPoints
    case 1
        x=rand(1,nbPoints)*10;
        mu=[1 2 3];
        muNew=medoids_centersAssignmentSearchSingle(x,mu,kCenter);
        plot(x,zeros(1,nbPoints),'bx')
        hold on
        plot(mu,zeros(1,3),'rs')
        plot(muNew,ones(1,3),'go')
        hold off
        axis([0 10 -5 5])
    case 2
        x=rand(2,nbPoints)*10;
        mu=[
            1 8 3
            1 5 3
            ];
        muNew=medoids_centersAssignmentSearchSingle(x,mu,kCenter,opts{:});
        plotPoints(x,'bx')
        hold on
        plotPoints(mu,'rs')
        plotPoints(muNew,'go')
        if flagSurfCost
            medoids_centersAssignmentSearchSingle_plot(x,mu,kCenter,opts{:})
            xlabel('x')
            ylabel('y')
            %view(0,0)
        end
        hold off
end
disp('Old centers')
disp(mu)
disp('New centers')
disp(muNew)

