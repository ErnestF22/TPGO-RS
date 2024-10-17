function POCbearingNetworkCostPair_test
global funs
%resetRands(4)

funs=bearingCostFunctions('cosine');
f=funs.f;
df=funs.df;
ddf=funs.ddf;

switch 1
    case 1
        [Ti,~,Ti0,vi]=real_randGeodFun([0;0]);
        [Tj,~,Tj0,vj]=real_randGeodFun([1;0],'speed',rand);
    case 2
        Ti0=[0;0];
        Tj0=[1;0];
        vi=[0;1];
        vj=[0;1]/2;
        Ti=@(t) Ti0+t*vi;
        Tj=@(t) Tj0+t*vj;
end
Yijtruth=bearingCompute(Ti0,Tj0);
z=zeros(2,1);
t=linspace(0,5,100);

switch 4
    case 1
        %check derivative
        check_der(@(t) costAndDer(Ti(t),Tj(t),vi,vj,Yijtruth), 'function',t)
    case 2
        %first split of derivative in 4 components
        I=eye(2);
        Yij=@(t) bearingCompute(Ti(t),Tj(t));
        c=@(t) bearingComputeCosine(Yij(t),Yijtruth);
        fc=@(t) f(c(t));
        dfc=@(t) df(c(t));

        dci1=@(t)  -2*fc(t)*(vi'*Yij(t));
        dci2=@(t)  -2*dfc(t)*vi'*(I-Yij(t)*Yij(t)')*Yijtruth;
        dcj1=@(t)   2*fc(t)*(vj'*Yij(t));
        dcj2=@(t)   2*dfc(t)*vj'*(I-Yij(t)*Yij(t)')*Yijtruth;
        dc=@(t) dci1(t)+dci2(t)+dcj1(t)+dcj2(t);

        plotfun(@(t) der(Ti(t),Tj(t),vi,vj,Yijtruth),t)
        hold on
        plotfun(@(t) dc(t),t,'rx')
        plotfun(@(t) dci1(t),t,'c--')
        plotfun(@(t) dci2(t),t,'m--')
        plotfun(@(t) dcj1(t),t,'k--')
        plotfun(@(t) dcj2(t),t,'g--')
        plotfun(@(t) dci1(t)+dci2(t),t,'c')
        plotfun(@(t) dcj1(t)+dcj2(t),t,'k')
        hold off
        legend('der','dc','dci1','dci2','dcj1','dcj2','dci','dcj')
    case 3
        %second split of derivative in 4 components
        I=eye(2);
        Yij=@(t) bearingCompute(Ti(t),Tj(t));
        c=@(t) bearingComputeCosine(Yij(t),Yijtruth);
        fc=@(t) f(c(t));
        dfc=@(t) df(c(t));
        gc=@(t) (fc(t)-dfc(t)*c(t));

        dci1=@(t)  -2*gc(t)*(vi'*Yij(t));
        dci2=@(t)  -2*dfc(t)*(vi'*Yijtruth);
        dcj1=@(t)   2*gc(t)*(vj'*Yij(t));
        dcj2=@(t)   2*dfc(t)*(vj'*Yijtruth);
        dc=@(t) dci1(t)+dci2(t)+dcj1(t)+dcj2(t);

        legendText={};
        plotfun(@(t) der(Ti(t),Tj(t),vi,vj,Yijtruth),t);
        legendText=[legendText 'der'];
        hold on
        plotfun(@(t) dc(t),t,'rx')
        legendText=[legendText 'dc'];

        plotfun(@(t) dci1(t)+dcj1(t),t,'c')
        legendText=[legendText 'dc1'];

        plotfun(@(t) dci2(t)+dcj2(t),t,'k')
        legendText=[legendText 'dc2'];

        plotfun(@(t) dci1(t),t,'c--')
        legendText=[legendText 'dci1'];
% 
%         plotfun(@(t) dci2(t),t,'m--')
%         legendText=[legendText 'dci2'];
% 
        plotfun(@(t) dcj1(t),t,'k--')
        legendText=[legendText 'dcj1'];
% 
%         plotfun(@(t) dcj2(t),t,'g--')
%         legendText=[legendText 'dcj2'];

        hold off
        legend(legendText)
    case 4
        %check equivalence under time-dependent translation
        XEval=@(t) Tj(t)-Ti(t);
        dXEval=vj-vi;
        X=[0;0];
        YGoal=bearingCompute(XEval(0),X);
        
        %first check cost and derivative of single cost
        f1=@(t) costAndDerSingle(XEval(t),X,dXEval,YGoal);
        f2=@(t) costAndDer(Ti(t),Tj(t),vi,vj,Yijtruth);
        figure(1)
        check_der(f1,'function',t)
        title('Check single term function and derivative')
        
        %the compare the two versions of the cost
        figure(2)
        plotfun(f1,t,'b');
        hold on
        plotfun(f2,t,'rx');
        hold off
        title('Check single versus pair interpretation')
end


function [c,dc]=costAndDer(Ti,Tj,vi,vj,Yijtruth)
global funs
[Yij,nYij]=bearingCompute(Ti,Tj);
[c,gradc]=bearingNetworkCostPair(Yij,Yijtruth,nYij,funs);
c=c;
dc=[vi;vj]'*gradc(:);

function [c,dc]=costAndDerSingle(XEval,X,dXEval,YGoal)
global funs
[YEval,nYEval]=bearingCompute(XEval,X);
[c,gradc]=bearingCostGeneral(YEval,YGoal,nYEval,funs);
dc=gradc'*dXEval;

function dc=der(Ti,Tj,vi,vj,Yijtruth)
global funs
Yij=bearingCompute(Ti,Tj);
gradc=bearingNetworkCostPair_grad(Yij,Yijtruth,funs);
dc=[vi;vj]'*gradc(:)*2;

