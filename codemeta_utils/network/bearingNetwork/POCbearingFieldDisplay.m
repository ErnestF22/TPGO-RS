function POCbearingFieldDisplay
L=10;
NGrid=20;
x0=[0;0];
xg=[1;0];

%funs=bearingCostFunctions('squared');
funs=consensus_rot3_almostGlobal_functions('type','huber','b',0.1);
y=@(x) bearingCompute(x,x0);
yg=y(xg);
dx=@(x) bearingField(y(x),yg,funs);

x=linspace(-L/2,L/2,NGrid);
[xGrid1,xGrid2]=meshgrid(x,x);
xGrid=cat(3,xGrid1,xGrid2);

yGrid=evalfunVec(y,xGrid);
dxGrid=evalfunVec(dx,xGrid);

quiver(xGrid(:,:,1),xGrid(:,:,2),dxGrid(:,:,1),dxGrid(:,:,2))
axis equal
