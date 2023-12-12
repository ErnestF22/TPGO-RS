function cartBearingCost_testLinePOC
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
v=cnormalize([-0.78;0.62]);
XEval=@(t) XEval0+t*v;
disp(v)

%check_der(@(t) cost(XEval(t),X,XGoal,v),'function',linspace(0,10))

Yg=-cnormalize(XGoal);
Ye=@(t) -cnormalize(XEval(t));
c=@(t) Yg'*Ye(t);
a=@(t) acos(c(t));
f=@(t) cost(XEval(t),X,XGoal,v);
df=@(t) -0.5*a(t)^2*v'*Ye(t)+a(t)/sqrt(1-c(t)^2)*v'*(eye(2)-(Ye(t)*Ye(t)'))*Yg;
%df=@(t) -0.5*a(t)^2*v'*Ye(t)+a(t)/sqrt(1-c(t)^2)*v'*Yg-a(t)/sqrt(1-c(t)^2)*v'*Ye(t)*c(t);
%df=@(t) (-c(t)*a(t)/sqrt(1-c(t)^2)-a(t)^2/2)*(v'*Ye(t))+a(t)/sqrt(1-c(t)^2)*(v'*Yg);
cosinca=@(t) a(t)/sin(a(t));
%df=@(t) (-c(t)*cosinca(t)-a(t)^2/2)*(v'*Ye(t))+cosinca(t)*(v'*Yg);
cota=@(t) cot(a(t));
%df=@(t) (-a(t)*cota(t)-a(t)^2/2)*(v'*Ye(t))+cosinca(t)*(v'*Yg);

t=linspace(0,0.7);
nplots=[1 12 13];

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
            hold off
            legend('v^TYe','Ye^TYg','v^TYg')
        case 2
            check_fun(@(t) a(t)/sqrt(1-c(t)^2), @(t) a(t)/sin(a(t)),t)
        case 3
            plotfun(@(t) (-a(t)*cota(t)-a(t)^2/2)*(v'*Ye(t)),t)
            hold on
            plotfun(@(t) cosinca(t)*(v'*Yg),t,'r')
            hold off
            legend('First term','Second term')
        case 4
            plotfun(@(t) (-a(t)*cota(t)-a(t)^2/2),t)
            hold on
            plotfun(@(t) cosinca(t),t,'r')
            hold off
            legend('First term coeff','Second term coeff')
        case 5
            dfLB=@(t) (-a(t)*cota(t)-a(t)^2/2)*(v'*Ye(t))+cosinca(t)*(v'*Yg);
            plotfun(df,t)
            hold on
            plotfun(dfLB,t,'r')
            hold off
            legend('Derivative','Lower bound on derivative')
        case 6
            plotfun(a,t)
            legend('a(t)')
        case 7
            plotfun(c,t)
            hold on
            plotfun(@(t) v(1),t,'r')
            hold off
            legend('c(t)','v1')

        case 8
            plotfun(@(t) (-a(t)*cota(t)-a(t)^2/2)-cosinca(t),t,'r')
            legend('Difference of coefficients')
        case 9
            plotfun(@(t) -0.5*a(t)^2*v'*Ye(t),t)
            hold on
            plotfun(@(t) +cosinca(t)*v'*(eye(2)-(Ye(t)*Ye(t)'))*Yg,t,'r')
            plot(-v(1),-0.5*a(-v(1))^2*v'*Ye(-v(1)),'o')
            hold off
            legend('df1 f2','f1 df2')
        case 10
            plotfun(@(t) -v'*Ye(t),t)
            hold on
            plotfun(@(t) +a(t)/sqrt(1-c(t)^2)*v'*(eye(2)-(Ye(t)*Ye(t)'))*Yg,t,'r')
            plot(-v(1),-v'*Ye(-v(1)),'o')
            hold off
            legend('df1','f1 . df2')
        case 11
            plotfun(@(t) (-a(t)*cota(t))*(v'*Ye(t)),t)
            hold on
            plotfun(@(t) (-a(t)^2/2)*(v'*Ye(t)),t,'c')
            plotfun(@(t) cosinca(t)*(v'*Yg),t,'r')
            hold off
            legend('First term, firs half','First term, second half','Second term')
        case 12
            plotfun(@(t) 0.5*a(t)^2*v'*Ye(t),t)
            hold on
            plotfun(@(t) cosinca(t)*(v'*Yg-(Ye(t)'*Yg)*v'*Ye(t)),t,'r')
            plotfun(@(t) cosinca(t)*(v'*Yg-v'*Ye(t)),t,'m')
            plot(-v(1),-v'*Ye(-v(1)),'o')
            hold off
            legend('-df1 f2','f1 df2', 'LB f1 df2')
        case 13
            plotfun(@(t) (v'*Yg-v'*Ye(t))/(v'*Ye(t))-0.5*a(t)*sin(a(t)),t,'m')
            legend('LB f1 df2- (-df1 f2)')
        case 14
            check_fun(@(t) (v'*Yg-v'*Ye(t))/(v'*Ye(t))-0.5*a(t)^2/cosinca(t),...
                @(t) (v'*Yg-v'*Ye(t))/(v'*Ye(t))-0.5*a(t)*sin(a(t)),t)
        case 15
            plotfun(@(t) 0.5*a(t)^2/cosinca(t),t)
            hold on
            plotfun(@(t) (v'*Yg-v'*Ye(t))/(v'*Ye(t)),t,'m')
            plot(-v(1),-v'*Ye(-v(1)),'o')
            hold off
            legend('smaller','larger')
    end
end

function [c,dc]=cost(XEval,X,XGoal,v)
[c,gradc]=cartBearingCost(XEval,X,XGoal);
dc=gradc'*v;

function ignore
df1f2=-0.5*a(t)^2*v'*Ye(t)
f1df2=+cosinca(t)*(v'*Yg-(Ye(t)'*Yg)*v'*Ye(t))
