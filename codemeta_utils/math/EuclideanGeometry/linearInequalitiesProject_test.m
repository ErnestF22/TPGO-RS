function linearInequalitiesProject_test
nbTrials=20;
C=[1 0; -1 0; 0 1; 0 -1];
d=[0;-1;0;-1];

figure(1)
hold on
for iTrial=1:nbTrials
    x=2*randn(2,1);
    y=linearInequalitiesProject(x,C,d);
    plotPoints(y,'bx')
    plotPoints(x,'rx')
    plotLines(x,y,'k')
end
plot([0 1 1 0 0],[0 0 1 1 0],'k:')
hold off
axis equal

    
