function POCCartControl
%state: x=[x;y;theta]
%input: u=[v;w]

fieldNumber=2;
controlNumber=4;

x0=[1;0;-pi/4];
TFinal=50;
D=0.1;

ctheta=@(x) cos(x(3));
stheta=@(x) sin(x(3));
f=@(x,u) [ctheta(x) 0; stheta(x) 0; 0 1]*u;

%[t,x]=ode45(@(t,x) f(x,[1;1]), [0 5], [0;0;0]);
%x=x';

switch fieldNumber
    case 1
        fRef1=@(x) 1;
        fRef2=@(x) 0;
    case 2
        d=@(x) norm([x(1); x(2)]);
        
        fRef1=@(x) d(x)*cos(atan2(x(2),x(1))+pi/2);
        fRef2=@(x) d(x)*sin(atan2(x(2),x(1))+pi/2);
       
        dfRef1=@(x) -d(x)*sin(atan2(x(2),x(1))+pi/2);
        dfRef2=@(x) d(x)*cos(atan2(x(2),x(1))+pi/2);

end

switch controlNumber
    case 1
        %hack
        u=@(x) [[ctheta(x);stheta(x)]'*[fRef1(x);fRef2(x)]; -40*angleDiff(x(3),atan2(fRef2(x),fRef1(x)))];
    case 2
        %input-output linearization around trajectory (not working)
        A=@(x) [ctheta(x) -D*stheta(x); stheta(x) D*ctheta(x)];
        r=@(x) [dfRef1(x); dfRef2(x)];
        u=@(x) A(x)\r(x);
    case 3
        A=@(x) [ctheta(x) -D*stheta(x); stheta(x) D*ctheta(x)];
        u=@(x) A(x)\[fRef1(x);fRef2(x)];
    case 4
        k=2;
        fRef12Sq=@(x) fRef1(x)^2+fRef2(x)^2;
        thetaRef=@(x) atan2(fRef2(x),fRef1(x));
        u=@(x) [sqrt(fRef12Sq(x)); (fRef1(x)*dfRef2(x)-dfRef1(x)*fRef2(x))/fRef12Sq(x)-k*angleDiff(x(3),thetaRef(x))];
end

[t,x]=ode45(@(t,x) f(x,u(x)), [0 TFinal], x0);
x=x';

% A=@(x,u) [0 0 -u(1)*stheta(x); 0 0 u(1)*ctheta(x); 0 0 0];
% B=@(x) [ctheta(x) 0; stheta(x) 0; 0 1];
% 
% 
% 
% NSteps=1000;
% dt=0.01;
% t=dt*(0:NSteps-1);
% x=zeros(3,NSteps);
% u=0.4*[ones(1,NSteps); ones(1,NSteps)];
% 
% x(:,1)=[0;0;0];
% for it=2:NSteps
%     %use midpoint method for the integration
%     xprev=x(:,it-1);
%     uprev=u(:,it-1);
% %     xhalf=xprev+dt/2*(A(xprev,uprev)*xprev+B(xprev)*uprev);
% %     x(:,it)=xprev+dt*(A(xhalf,uprev)*xhalf+B(xhalf)*uprev);
%     x(:,it)=xprev+dt*(A(xprev,uprev)*xprev+B(xprev)*uprev);
% end

xMin=-5;
xMax=5;
xStep=0.7;
xGrid=[xMin:xStep:xMax];
[X,Y]=meshgrid(xGrid,xGrid);
Z=zeros([size(X) 2]);
for iX=1:length(xGrid)
    for iY=1:length(xGrid);
        Z(iX,iY,1)=fRef1([X(iX,iY);Y(iX,iY);0]);
        Z(iX,iY,2)=fRef2([X(iX,iY);Y(iX,iY);0]);
    end
end

subplot(2,1,1)
quiver(X,Y,Z(:,:,1),Z(:,:,2),'k')
hold on
plot(x(1,:),x(2,:))
axis equal
axis([-5 5 -5 5])
hold off

subplot(2,1,2)
plot(t,x(3,:));
