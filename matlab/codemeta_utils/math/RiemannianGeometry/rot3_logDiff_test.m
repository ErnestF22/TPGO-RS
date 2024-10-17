function rot3_logDiff_test
RBase=rot_randn(eye(3));
testNum=2;

switch testNum
    case 1
        [Rt,dRt,R0,dR0,v]=rot_randGeodFun(rot_randn(eye(3)));

        LogRt=@(t) logrot(RBase'*Rt(t));

        u=@(t) LogRt(t);

        Du=@(t) rot3_logDiff(RBase,Rt(t));

        du=@(t) Du(t)*v;

        check_der(u,du)
    case 2
        [R1t,dR1t,R10,dR10,v1]=rot_randGeodFun(rot_randn(eye(3)));
        [R2t,dR2t,R20,dR20,v2]=rot_randGeodFun(rot_randn(eye(3)),'speed',0);

        %R1ttR2t=@(t) R1t(t)'*R2t(t);
        %dR1ttR2t=@(t) -R1t(t)'*R2t(t)*hat(R2t(t)'*R1t(t)*v1);
        %check_der(R1ttR2t,dR1ttR2t,'angle')
        
        LogRt=@(t) logrot(RBase'*R1t(t)'*R2t(t));

        u=@(t) LogRt(t);

        Du=@(t) rot3_logDiff(RBase,R1t(t)'*R2t(t));

        du=@(t) Du(t)*(v2-R2t(t)'*R1t(t)*v1);

        check_der(u,du)
    case 3
        [R1t,dR1t,R10,dR10,v1]=rot_randGeodFun(rot_randn(eye(3)));
        [R2t,dR2t,R20,dR20,v2]=rot_randGeodFun(rot_randn(eye(3)),'speed',0);

        %R1ttR2t=@(t) R1t(t)'*R2t(t);
        %dR1ttR2t=@(t) -R1t(t)'*R2t(t)*hat(R2t(t)'*R1t(t)*v1);
        %check_der(R1ttR2t,dR1ttR2t,'angle')
        
        LogRt=@(t) logrot(RBase'*R1t(t)'*R2t(t));

        u=@(t) 0.5*LogRt(t)'*LogRt(t);

        DLogRt=@(t) rot3_logDiff(RBase,R1t(t)'*R2t(t));

        du=@(t) LogRt(t)'*DLogRt(t)*(v2-R2t(t)'*R1t(t)*v1);

        check_der(u,du)
        
end
