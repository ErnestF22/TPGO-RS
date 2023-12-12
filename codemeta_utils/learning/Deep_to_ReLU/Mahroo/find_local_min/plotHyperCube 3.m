function plotHyperCube
[X,Z] = meshgrid(3:0.5:7,3:0.5:6);
Y = 3*ones(size(X));
surf(X,Y,Z)
hold on
[X,Z] = meshgrid(3:0.5:7,5:0.5:7);
Y = 7*ones(size(X));
surf(X,Y,Z)
hold on
[X,Y] = meshgrid(3:0.1:7,3.5:0.1:7);
Z = 7*ones(size(X));
surf(X,Y,Z)
hold on
[X,Y] = meshgrid(3:0.5:7,3:0.5:6);
Z = 3*ones(size(X));
surf(X,Y,Z)
hold on
[Y,Z] = meshgrid(3:0.5:7,3:0.5:7);
X = 7*ones(size(Y));
surf(X,Y,Z)
hold on
[Y,Z] = meshgrid(3:0.5:7,3:0.5:7);
S = ones(size(Y));
X = 3*ones(size(Y));
surf(X,Y,Z)
hold on
[Y,X] = meshgrid(3:0.1:3.5,3:0.1:7);
Z = 2*Y;
surf(X,Y,Z)
hold on
[Y,X] = meshgrid(6:0.5:7,3:0.5:7);
Z = 2*Y-9;
surf(X,Y,Z)