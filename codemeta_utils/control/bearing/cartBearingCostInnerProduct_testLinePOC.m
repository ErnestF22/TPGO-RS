function cartBearingCostInnerProduct_testLinePOC
%resetRands();
X=zeros(2,1);
%d=randn;
d=1;
XGoal=[d;0];
XEval0=XGoal;
%v=cnormalize(randn(2,1));
%v=cnormalize([-1;0.00001]);
%v=[0;1];
%v=cnormalize([-0.5;1]);
%v=cnormalize([-0.78;0.62]);
v=cnormalize([0.78;0.62]);
XEval=@(t) XEval0+t*v;
disp(v)

%check_der(@(t) cost(XEval(t),X,XGoal,v),'function',linspace(0,10))

Yg=-cnormalize(XGoal);
Ye=@(t) -cnormalize(XEval(t));
c=@(t) Yg'*Ye(t);
f=@(t) cost(XEval(t),X,XGoal,v);
dfa=@(t) -0.5*(1-c(t))^2*v'*Ye(t);
dfb=@(t) (1-c(t))*v'*(eye(2)-(Ye(t)*Ye(t)'))*Yg;
%df=@(t) -0.5*(1-c(t))^2*v'*Ye(t)+(1-c(t))*v'*(eye(2)-(Ye(t)*Ye(t)'))*Yg;
%df=@(t) (1-c(t))*(-0.5*(1-c(t))*v'*Ye(t)+v'*(eye(2)-(Ye(t)*Ye(t)'))*Yg);
%df=@(t) (1-c(t))*(0.5*(c(t)-1)*v'*Ye(t)+v'*Yg-c(t)*v'*Ye(t));
df=@(t) (1-c(t))*(-0.5*(1+c(t))*v'*Ye(t)+v'*Yg);

t0=-v(1);
tMax=100;

dfLB=@(t) (1-c(t))*(...
    (v'*Yg>=0)*(...
        (t<=t0)*(-v'*Ye(t)+v'*Yg) +...
        (t> t0)*(v'*Yg))+...
    (v'*Yg< 0)*(...
        (0.5*(1+c(t))*(-v'*Ye(t))+v'*Yg))...
    );

t=linspace(0,tMax);
nplots=[1 3 7];

subplot(1+length(nplots),1,1)
check_der(f,df,t)
for iplot=1:length(nplots)
    subplot(1+length(nplots),1,1+iplot)
    switch nplots(iplot)
        case 1
            plotfun(@(t) v'*Ye(t),t)
            hold on
            plotfun(@(t) Ye(t)'*Yg,t,'c')
            plotfun(@(t) v'*Yg,t,'r')
            plotfun(@(t) -v'*Yg,t,'m')
            if t0>0 && t0<tMax
                plot(t0,0,'bo');
            end
            hold off
            legend('v^TYe','Ye^TYg','v^TYg','-v^TYg')
        case 2
            dfa=@(t) -0.5*(1-c(t))^2*v'*Ye(t);
            dfb=@(t) (1-c(t))*v'*(eye(2)-(Ye(t)*Ye(t)'))*Yg;
            plotfun(dfa,t)
            hold on
            plotfun(dfb,t,'r')
            hold off
            legend('dfa','dfb')
        case 3
            plotfun(@(t) 1-c(t),t)
            hold on
            plotfun(@(t) 1+c(t),t,'r')
            hold off
            legend('1-c(t)','1+c(t)')
        case 4
            dfa=@(t) -0.5*(1+c(t))*v'*Ye(t);
            dfb=@(t) v'*Yg;
            plotfun(dfa,t)
            hold on
            plotfun(dfb,t,'r')
            hold off
            legend('dfa','dfb')
        case 5
            plotfun(df,t)
            hold on
            plotfun(dfLB,t,'r')
            hold off
            legend('df','LB on df')
        case 6
            plotfun(@(t) df(t)-dfLB(t),t)
            legend('df - LB on df')
        case 7
            plotfun(dfa,t)
            hold on
            plotfun(dfb,t,'r')
            hold off
            legend('dfa','dfb')
    end
end

function [c,dc]=cost(XEval,X,XGoal,v)
[c,gradc]=cartBearingCostInnerProduct(XEval,X,XGoal);
dc=gradc'*v;
