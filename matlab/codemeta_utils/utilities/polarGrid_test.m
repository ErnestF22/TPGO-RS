function polarGrid_test
x=polarGrid(20,37,'thetaMin',-pi+deg2rad(30),'thetaMax',pi-deg2rad(30));

plotPoints(x)
axis equal