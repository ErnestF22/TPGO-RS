function essential_fromRT_Diff_test
switch 2
    case 1
        %test directional derivative of Q and E
        e3hat=[0 -1 0; 1 0 0; 0 0 0];

        [T1,~,~,vT1]=real_randGeodFun(randn(3,1));
        [T2,~,~,vT2]=real_randGeodFun(randn(3,1));
        [R1,~,~,~,vR1]=rot_randGeodFun();
        [R2,~,~,~,vR2]=rot_randGeodFun();

        Q=@(t) essential_fromRT(R1(t),T1(t),R2(t),T2(t));
        D=@(t) essential_fromRT_Diff(R1(t),T1(t),R2(t),T2(t));

        v=@(t) D(t)*[vR1;vT1;vR2;vT2];

        dQ=@(t) essential_hat(Q(t),v(t));

%         check_der(Q,dQ)
        
        E=@(t) essential_getE(Q(t));
        dE=@(t) essential_getR1(dQ(t))'*e3hat*essential_getR2(Q(t))...
            +essential_getR1(Q(t))'*e3hat*essential_getR2(dQ(t));
        
%         check_der(E,dE)
        
        dQProj=@(t) essential_tangentProj(Q(t),dQ(t));
        dEProj=@(t) essential_getR1(dQProj(t))'*e3hat*essential_getR2(Q(t))...
            +essential_getR1(Q(t))'*e3hat*essential_getR2(dQProj(t));

        check_der(E,dEProj)
        
%         Q1=@(t) essential_getR1(Q(t));
%         dQ1=@(t) essential_getR1(dQ(t));
%         Q2=@(t) essential_getR2(Q(t));
%         dQ2=@(t) essential_getR2(dQ(t));
        
    case 2
        %check horizontality of vectors
        T1=randn(3,1);
        T2=randn(3,1);
        R1=rot_randn();
        R2=rot_randn();
        Q=essential_fromRT(R1,T1,R2,T2);
        D=essential_fromRT_Diff(R1,T1,R2,T2);
        
        dQrand=@() essential_tangentProjVerticalScalar(Q,essential_hat(Q,D*randn(12,1)));
        plotfuntrials(dQrand,100)
end
        