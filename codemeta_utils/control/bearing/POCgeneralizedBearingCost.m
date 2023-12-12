function POCgeneralizedBearingCost
X=zeros(2,1);
d=1;
XGoal=[d;0];
Xe0=XGoal;
%v=cnormalize(randn(2,1));
%v=cnormalize([-1;0.00001]);
%v=[0;1];
%v=cnormalize([-0.5;1]);
%v=cnormalize([-0.78;0.62]);
%v=cnormalize([-0.85;0.2]);
v=cnormalize([0.78;0.62]);
Xe=@(t) Xe0+t*v;

t=linspace(0,1);

Yg=-cnormalize(XGoal);
Ye=@(t) -cnormalize(Xe(t));
c=@(t) Yg'*Ye(t);
nXe=@(t) norm(Xe(t));

if v'*Yg<0
    disp('g_g is negative')
else
    disp('g_g is positive')
end

switch 2
    case 1
        f=@(c) (c-1)^2/2;
        df=@(c) c-1;
    case 2
        f=@(c) acos(c)^2;
        df=@(c) -(2*acos(c))/(1 - c^2)^(1/2);
end
b1=@(c) f(c)-df(c)*c;
b2=@(c) df(c);

tc=linspace(-1,1);
figure(1)
subplot(2,1,1)
plotfun(b2,tc)
ylabel('b2')
xlabel('This should be negative')
subplot(2,1,2)
plotfun(@(t) b1(t)+b2(t),tc)
ylabel('b1+b2')
xlabel('This shold be negative')
% subplot(4,1,1)
% plotfun(b1,tc)
% ylabel('b1')
% xlabel('This should be positive?')
% title('Check conditions on f and df')
% subplot(4,1,4)
% plotfun(@(t) b1(t)-b2(t),tc)
% ylabel('b1-b2')
% %xlabel('This shold be negative')

phi=@(t) nXe(t)*f(c(t));
gradphi=@(t) -f(c(t))*Ye(t)-df(c(t))*(eye(2)-Ye(t)*Ye(t)')*Yg;
dphi=@(t) v'*gradphi(t);
dphib=@(t) -(f(c(t))-df(c(t))*c(t))*v'*Ye(t)-df(c(t))*v'*Yg;

dphic=@(t) -b1(c(t))*v'*Ye(t)-b2(c(t))*v'*Yg;

figure(2)
subplot(4,1,1)
check_der(phi,dphic,t)
title('Check that derivative is correct')
subplot(4,1,2)
plotfun(@(t) v'*Ye(t),t)
hold on
plotfun(@(t) v'*Yg,t,'r')
plotfun(c,t,'c')
plotfun(@(t) v(1),t,'c:')
hold off
legend('v^TY_e','v^TY_g','c','v_1')
subplot(4,1,3)
plotfun(@(t) b1(c(t)),t)
hold on
plotfun(@(t) b2(c(t)),t,'r')
hold off
legend('b1(c)','b2(c)')
subplot(4,1,4)
plotfun(@(t) abs(b1(c(t))*v'*Ye(t)),t)
hold on
plotfun(@(t) abs(b2(c(t))*v'*Yg),t,'r')
plotfun(@(t) abs(b1(c(t))*v'*Yg),t,'c')
plotfun(@(t) abs(b2(c(t))*v'*Ye(t)),t,'m')
hold off
legend('|b1(c) v^TY_e|','|b2(c) v^TY_g|','|b1(c) v^TY_g|','|b2(c) v^TY_e|')

ndphi=@(t) -dphic(t);

ndphiUB=@(t) (v'*Yg> 0)*(b1(c(t))+b2(c(t)))*v'*Yg...
            +(v'*Yg<=0)*(b1(c(t))*v'*Ye(t)-b1(c(t))*v'*Yg);

figure(3)
plotfun(@(t) ndphiUB(t)-ndphi(t),t)
title('Check upper bound on negative derivative (i.e., it is an u.b.)')
ylabel('This shoudl be positive')

figure(4)
plotfun( ndphiUB,t)
title('Check that the upper bound on negative derivative is negative (i.e., we can use the u.b.)')
ylabel('This should be negative')

figure(5)
check_fun(@(t) v'*Yg-c(t)*v'*Ye(t), @(t) v'*(eye(2)-Ye(t)*Ye(t)')*Yg,t)
