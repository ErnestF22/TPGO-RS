function POCessentialDistanceLog
resetRands(3)
t=linspace(-pi,pi);
switch 8
    case 1
        %explicit form for thetai
        e3=[0;0;1];
        R0=rot_randn();
        Rzt=@(t) rot(t*e3);
        
        Rt=@(t) R0*Rzt(t);
        DLog=@(t) rot3_logDiff(eye(3),Rt(t));
        
        thetat=@(t) rot_dist(eye(3),Rt(t));

        ai=R0(1,1)+R0(2,2);
        bi=R0(1,2)-R0(2,1);
        ci=R0(3,3);
        
        thetabt=@(t) acos(((ai*cos(t)+bi*sin(t)+ci)-1)/2);
        
        m=norm([ai;bi]);
        p=sign(ai)*acos(bi/m);
        
        thetact=@(t) acos(((m*sin(t+p)+ci)-1)/2);

        tBreak=3/2*pi-p;
        
        disp(-1-ci+m)
        
        subplot(2,1,1)
        plotfun(thetat,t,'bx')
        hold on
        plotfun(thetabt,t,'ro')
        plotfun(thetact,t,'g')
        plot(tBreak,thetact(tBreak),'+','MarkerSize',10)
        hold off
        subplot(2,1,2)
        plotfun(@(t) ai*cos(t)+bi*sin(t),t)
        hold on
        plotfun(@(t) -1-ci,t,'--')
        hold off

    case 2
        e3=[0;0;1];
        %R0=eye(3);
        R0=rot_randn();
        %R0=rot_randn(eye(3),0.3);
        Rzt=@(t) rot(t*e3);
        
        Rt=@(t) R0*Rzt(t);
        Logt=@(t) logrot(Rt(t));
        DLogt=@(t) rot3_logDiff(eye(3),Rt(t));
        thetat=@(t) rot_dist(eye(3),Rt(t));
        ut=@(t) Logt(t)/thetat(t);

        ft=@(t) thetat(t)^2/2;
        dft=@(t) Logt(t)'*e3;
        ddft=@(t) e3'*DLogt(t)*e3;
        ct=@(t) thetat(t)/2*cot(thetat(t)/2);
        ddfbt=@(t) (1-ct(t))*(e3'*ut(t))^2+ct(t);

%         figure(1)
%         check_der(ft,dft,'angle')
%         figure(2)
%         check_der(dft,ddfbt,'angle')
%         plotfun(@(t) (e3'*ut(t))^2,'angle')
        plotfun(ct,'angle')
%         feval=evalfun(ft,t);
%         ceval=evalfun(ct,t);
%         plot(ceval,feval)
    case 3
        e3=[0;0;1];
        Q10=rot_randn([],[],2);
        Q110=Q10(:,:,1);
        Q120=Q10(:,:,2);
        Q20=rot_randn([],[],2);
        Q210=Q20(:,:,1);
        Q220=Q20(:,:,2);
        Rzt=@(t) rot(t*e3);
        
        Q21t=@(t) Q20(:,:,1)*Rzt(t);
        Q22t=@(t) Q20(:,:,2)*Rzt(t);
        
        ft=@(t) 0.5*(rot_dist(Q110,Q21t(t))^2+rot_dist(Q120,Q22t(t))^2);
        %dft=@(t) e3'*(logrot(Q110'*Q21t(t))+logrot(Q120'*Q22t(t)));
        
        tBreak1=discontinuityDistance(Q110'*Q210);
        tBreak2=discontinuityDistance(Q120'*Q220);
        
        %first min
        tSearch1=tBreak1;
        tSearch2=tBreak2;
        if tSearch1>tSearch2
            tSearch1=tSearch1-2*pi;
        end
        [tMin1,fMin1]=fminbnd(ft,tSearch1,tSearch2);
        tMin1=modAngle(tMin1);
        tSearch1=tSearch1+2*pi;
        [tMin2,fMin2]=fminbnd(ft,tSearch2,tSearch1);
        tMin2=modAngle(tMin2);
        
        if fMin1<fMin2
            tMin=tMin1;
            fMin=fMin1;
        else
            tMin=tMin2;
            fMin=fMin2;
        end
        plotfun(ft,'angle')
        hold on
        plot([tBreak1 tBreak2],[ft(tBreak1) ft(tBreak2)],'k+','MarkerSize',10)
        plot([tMin1 tMin2],[fMin1 fMin2],'rx','MarkerSize',10)
        plot(tMin,fMin,'ro','MarkerSize',10)
        hold off
        
        Q21Comp=Q21t(tMin);
        Q22Comp=Q22t(tMin);
        v1=logrot(Q110'*Q21Comp);
        v2=logrot(Q120'*Q22Comp);
        disp(v1(3)+v2(3))
        
    case 4
        %same as 3 but with left rotation multiplication
        e3=[0;0;1];
        Q10=rot_randn([],[],2);
        Q110=Q10(:,:,1);
        Q120=Q10(:,:,2);
        Q20=rot_randn([],[],2);
        Q210=Q20(:,:,1);
        Q220=Q20(:,:,2);
        Rzt=@(t) rot(t*e3);
        
        Q21t=@(t) Rzt(t)*Q20(:,:,1);
        Q22t=@(t) Rzt(t)*Q20(:,:,2);
        
        ft=@(t) 0.5*(rot_dist(Q110,Q21t(t))^2+rot_dist(Q120,Q22t(t))^2);
        %dft=@(t) e3'*(logrot(Q110'*Q21t(t))+logrot(Q120'*Q22t(t)));
        
        tBreak1=discontinuityDistance(Q210*Q110');
        tBreak2=discontinuityDistance(Q220*Q120');
        
        %first min
        tSearch1=tBreak1;
        tSearch2=tBreak2;
        if tSearch1>tSearch2
            tSearch1=tSearch1-2*pi;
        end
        [tMin1,fMin1]=fminbnd(ft,tSearch1,tSearch2);
        tMin1=modAngle(tMin1);
        tSearch1=tSearch1+2*pi;
        [tMin2,fMin2]=fminbnd(ft,tSearch2,tSearch1);
        tMin2=modAngle(tMin2);
        
        if fMin1<fMin2
            tMin=tMin1;
            fMin=fMin1;
        else
            tMin=tMin2;
            fMin=fMin2;
        end
        plotfun(ft,'angle')
        hold on
        plot([tBreak1 tBreak2],[ft(tBreak1) ft(tBreak2)],'k+','MarkerSize',10)
        plot([tMin1 tMin2],[fMin1 fMin2],'rx','MarkerSize',10)
        plot(tMin,fMin,'ro','MarkerSize',10)
        hold off
        
        Q21Comp=Q21t(tMin);
        Q22Comp=Q22t(tMin);
        v1=logrot(Q110'*Q21Comp);
        v2=logrot(Q120'*Q22Comp);
        disp(v1(3)+v2(3))
        
    case 5
        %explicit form for function and derivatives
        e3=[0;0;1];
        Q10=rot_randn([],[],1);
        Q20=rot_randn([],[],1);
        Rzt=@(t) rot(t*e3);
        
        Q2t=@(t) Rzt(t)*Q20;
        
        ft=@(t) 0.5*rot_dist(Q10,Q2t(t))^2;
        dft=@(t) e3'*Q10*logrot(Q10'*Q2t(t));
        ddft=@(t) e3'*Q10*rot3_expDiffInvMat(eye(3),Q10'*Q2t(t))*Q10'*e3;
        
        disp('Check analytic, non-closed form derivatives')
        check_der(ft,dft,'angle')
        check_der(dft,ddft,'angle')
        
        R0=Q20*Q10';
        [tBreak,a,b,c,m,p]=discontinuityDistance(R0);
        
        %thetat=@(t) acos((a*cos(t)+b*sin(t)+c-1)/2);
        thetat=@(t) acos((m*sin(t+p)+c-1)/2);
        dthetat=@(t) -(m*cos(t+p))/(2*sin(thetat(t)));

        fbt=@(t) 0.5*thetat(t)^2;
        dfbt=@(t) thetat(t)*dthetat(t);
        ddfbt=@(t) -m/(2*sin(thetat(t))^2)*((dthetat(t)*cos(t+p)-thetat(t)*sin(t+p)) *sin(thetat(t))...
            -thetat(t)*cos(t+p)*cos(thetat(t))*dthetat(t));
        
        disp('Comparing form based on theta with analytic form')
        check_fun(ft,fbt,'angle')
        check_fun(dft,dfbt,'angle')
        check_fun(ddft,ddfbt,'angle')
        
        disp('Comparing closed form for R-R'', u and ezai''*u')
        skew=@(A) (A-A')/2;
        RRtveet=@(t) vee(skew(Q10'*Q2t(t)));
        %I=eye(3);
        %trERRR=@(k,t) trace(hat(I(:,k))*Q10'*Q2t(t));
        %RRtveebt=@(t) -[trERRR(1,t) trERRR(2,t) trERRR(3,t)]/2;
        av1=Q10(1,2)*Q20(2,3) + Q10(2,3)*Q20(1,2) - Q10(1,3)*Q20(2,2) - Q10(2,2)*Q20(1,3) ;
        bv1=Q10(1,3)*Q20(1,2) - Q10(1,2)*Q20(1,3) + Q10(2,3)*Q20(2,2) - Q10(2,2)*Q20(2,3) ;
        cv1=Q10(3,3)*Q20(3,2) - Q10(3,2)*Q20(3,3);
        av2=Q10(1,3)*Q20(2,1) + Q10(2,1)*Q20(1,3) - Q10(1,1)*Q20(2,3) - Q10(2,3)*Q20(1,1);
        bv2=Q10(1,1)*Q20(1,3) - Q10(1,3)*Q20(1,1) + Q10(2,1)*Q20(2,3) - Q10(2,3)*Q20(2,1);
        cv2=Q10(3,1)*Q20(3,3) - Q10(3,3)*Q20(3,1);
        av3=Q10(1,1)*Q20(2,2) + Q10(2,2)*Q20(1,1) - Q10(1,2)*Q20(2,1) - Q10(2,1)*Q20(1,2);
        bv3=Q10(1,2)*Q20(1,1) - Q10(1,1)*Q20(1,2) - Q10(2,1)*Q20(2,2) + Q10(2,2)*Q20(2,1);
        cv3=Q10(3,2)*Q20(3,1) - Q10(3,1)*Q20(3,2);
        RRtveebt=@(t) ([av1;av2;av3]*sin(t)+[bv1;bv2;bv3]*cos(t)+[cv1;cv2;cv3])/2;
        check_fun(RRtveet,RRtveebt,'angle')
        ut=@(t) logrot(Q10'*Q2t(t))/thetat(t);
        ubt=@(t) RRtveet(t)/sin(thetat(t));
        check_fun(ut,ubt,'angle')
        eztut=@(t) e3'*Q10*ut(t);
        eztubt=@(t) -m*cos(t+p)/(2*sin(thetat(t)));
        check_fun(eztut, eztubt,'angle')
        
        disp('Comparing closed form for df and ddf')
        check_fun(dft,@(t) dfc(m,p,c,t),'angle')
        check_fun(ddft,@(t) ddfc(m,p,c,t),'angle')
        
        disp('Comparing alternative way to compute sin(theta), cot(theta)')
        st=@(t) sin(thetat(t));
        sbt=@(t) sqrt(1-((m*sin(t+p)+c-1)/2)^2);
        check_fun(st,sbt,'angle')
        cott=@(t) cot(thetat(t));
        cotbt=@(t) (m*sin(t+p)+c-1)/2/sqrt(1-((m*sin(t+p)+c-1)/2)^2);
        check_fun(cott,cotbt,'angle')
     case 6
        %check analytic value for df around discontinuity
        e3=[0;0;1];
        Q10=rot_randn([],[],1);
        Q20=rot_randn([],[],1);
        Rzt=@(t) rot(t*e3);
        
        Q2t=@(t) Rzt(t)*Q20;
        
        ft=@(t) 0.5*rot_dist(Q10,Q2t(t))^2;
        dft=@(t) e3'*Q10*logrot(Q10'*Q2t(t),'quaternion');
        ddft=@(t) e3'*Q10*rot3_expDiffInvMat(eye(3),Q10'*Q2t(t))*Q10'*e3;
        
        R0=Q20*Q10';
        [tBreak,a,b,c,m,p]=discontinuityDistance(R0);
        
        
        Q2Break=Q2t(tBreak);
        v0=rot3_LogPi(Q10'*Q2Break);
        dftBreakp=e3'*Q10*v0;
        dftBreakm=-dftBreakp;

%         figure(1)
%         plotfun(@(t) logrot(Q10'*Q2t(t),'quaternion'),'angle')
%         hold on
%         plot(tBreak, v0,'x')
%         plot(tBreak, -v0,'x')
%         hold off
%         
%         figure(2)
%         plotfun(@(t) Q10*logrot(Q10'*Q2t(t),'quaternion'),'angle')
%         hold on
%         plot(tBreak, Q10*v0,'x')
%         plot(tBreak, -Q10*v0,'x')
%         hold off
%         
%         figure(3)
        plotfun(dft,'angle')
        hold on
        plot([tBreak tBreak],[dftBreakp dftBreakm],'x')
        hold off
    case 7
        %test Newton method (1 rotation)
        e3=[0;0;1];
        Q10=rot_randn([],[],1);
        Q20=rot_randn([],[],1);
        Rzt=@(t) rot(t*e3);
        
        Q2t=@(t) Rzt(t)*Q20;
        
        R0=Q20*Q10';
        [tBreak,a,b,c,m,p]=discontinuityDistance(R0);
        tMin=dfNewtonSingle(m,p,c,0);
        
        ft=@(t) 0.5*rot_dist(Q10,Q2t(t))^2;
        plotfun(ft,'angle')
        hold on
        plot(tMin,ft(tMin),'x')
        hold off
    case 8
        %test Newton method (2 rotations)
        e3=[0;0;1];
        Q10=rot_randn([],[],2);
        Q110=Q10(:,:,1);
        Q120=Q10(:,:,2);
        Q20=rot_randn([],[],2);
        Q210=Q20(:,:,1);
        Q220=Q20(:,:,2);
        Rzt=@(t) rot(t*e3);
        
        Q21t=@(t) Q20(:,:,1)*Rzt(t);
        Q22t=@(t) Q20(:,:,2)*Rzt(t);
        
        ft=@(t) 0.5*(rot_dist(Q110,Q21t(t))^2+rot_dist(Q120,Q22t(t))^2);
        %dft=@(t) e3'*(logrot(Q110'*Q21t(t))+logrot(Q120'*Q22t(t)));
        
        [tBreak1,a1,b1,c1,m1,p1]=discontinuityDistance(Q110'*Q210);
        [tBreak2,a2,b2,c2,m2,p2]=discontinuityDistance(Q120'*Q220);
        tMin=dfNewton(m1,p1,c1,m2,p2,c2,0);
        
        figure(1)
        check_der(ft,@(t) dfc(m1,p1,c1,t)+dfc(m2,p2,c2,t),'angle')
        hold on
        plot(tMin,ft(tMin),'mx')
        plot(tMin,0,'mx')
        hold off
        
end

function d=dfc(m,p,c,t)
theta=acos((m*sin(t+p)+c-1)/2);
dtheta= -(m*cos(t+p))/(2*sin(theta));
d=theta*dtheta;

function dd=ddfc(m,p,c,t)
theta=acos((m*sin(t+p)+c-1)/2);
eztu=-m*cos(t+p)/(2*sin(theta));
dd= eztu^2+theta/2*cot(theta/2)*(1-eztu^2);

function tMin=dfNewton(m1,p1,c1,m2,p2,c2,tMin)
tolDist=1e-12;
while true
    d=dfc(m1,p1,c1,tMin)+dfc(m2,p2,c2,tMin);
    dd=ddfc(m1,p1,c1,tMin)+ddfc(m2,p2,c2,tMin);
    tOld=tMin;
    tMin=tOld-d/dd;
    if abs(tMin-tOld)<tolDist
        break
    end
end

function tMin=dfNewtonSingle(m,p,c,tMin)
tolDist=1e-12;
while true
    d=dfc(m,p,c,tMin);
    dd=ddfc(m,p,c,tMin);
    tOld=tMin;
    tMin=tOld-d/dd;
    if abs(tMin-tOld)<tolDist
        break
    end
end

function [tBreak,a,b,c,m,p]=discontinuityDistance(R0)
a=R0(1,1)+R0(2,2);
b=R0(1,2)-R0(2,1);
c=R0(3,3);

m=norm([a;b]);
p=sign(a)*acos(b/m);

tBreak=modAngle(3/2*pi-p);

function a=modAngle(a)
a=mod(a+pi,2*pi)-pi;
