function POCControlCost
L=10;
NGrid=100;
NGridSub=20;

optsX={'MarkerSize',10};

switch 2
    case 1
        NX=10;
        X=L*rand(2,NX);
        %XGoal=[L/2;L/2];
        XGoal=L*rand(2,1);
    case 2
        NX=1;
        X=[5;5];
        XGoal=[6;5];
end
uVec=ones(1,NX);


x=linspace(0,L,NGrid);
[gridX,gridY]=meshgrid(x,x);
xSub=linspace(0,L,NGridSub);
[gridXSub,gridYSub]=meshgrid(xSub,xSub);

C=zeros(NGrid);
D=zeros([NGridSub NGridSub 2]);

for iX=1:NGrid
    for iY=1:NGrid
        XEval=[gridX(iX,iY);gridY(iX,iY)];
        C(iX,iY)=cartBearingCost(XEval,X,XGoal);
    end
end
for iX=1:NGridSub
    for iY=1:NGridSub
        XEval=[gridXSub(iX,iY);gridYSub(iX,iY)];
        [~,D(iX,iY,:)]=cartBearingCost(XEval,X,XGoal);
    end
end
CMax=max(C(:));
switch 2
    case 2
        contour(gridX,gridY,C,20)
        hold on
        plot(X(1,:),X(2,:),'r*',optsX{:})
        plot(XGoal(1),XGoal(2),'g*',optsX{:})
        quiver(gridXSub,gridYSub,-D(:,:,1),-D(:,:,2))
        hold off
    case 3
        surf(gridX,gridY,C,'EdgeColor','none')
        hold on
        plot3(X(1,:),X(2,:),CMax*uVec,'r*',optsX{:})
        plot3(XGoal(1),XGoal(2),CMax,'g*',optsX{:})
        hold off
        view(2)
end
axis equal
